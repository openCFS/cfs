#include "FeFunctions.hh"

#include <boost/bind.hpp>
#include <boost/tr1/type_traits.hpp>

#include "PDE/SinglePDE.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/Assemble.hh"
#include "BaseFE.hh"
#include "H1/H1Elems.hh"
#include "HCurl/HCurlElems.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"

namespace CoupledField {
DECLARE_LOG(fefunc)
 DEFINE_LOG(fefunc, "feFunction")
 
  BaseFeFunction::BaseFeFunction(){
    
    fctId_ = NO_FCT_ID;
    mHandle_ = domain->GetMathParser()->GetNewHandle();
    algsys_ = NULL;
    
    // initialize members of coefficient function
    dependType_ = CoefFunction::GENERAL;
    isAnalytic_ = false;
    dimType_ = NO_DIM;

  }
  
  BaseFeFunction::~BaseFeFunction(){
    
    // delete all mapping contexts
    std::map<std::string, MapContext*>::iterator it = entityCtx_.begin(); 
    for( ; it != entityCtx_.end(); ++it ) {
      delete it->second;
    }
    entityCtx_.clear();
  }
  
  
  shared_ptr<ResultInfo> BaseFeFunction::GetResultInfo(){
    return result_;
  }
  
  void BaseFeFunction::SetFeSpace( shared_ptr<FeSpace> space ){
    feSpace_ = space;
    
    // consistency check: if we have a HCurl space,
    // even a scalar unknown results in a vector-valued
    // function
    if( space->GetSpaceType() == FeSpace::HCURL ) {
      dimType_ = CoefFunction::VECTOR;
    }
        
  }
  
  shared_ptr<FeSpace> BaseFeFunction::GetFeSpace(){
    return feSpace_;
  }
  
  shared_ptr<BaseTimeScheme> BaseFeFunction::GetTimeScheme(){
    return timeScheme_;
  }

  void BaseFeFunction::SetTimeScheme(shared_ptr<BaseTimeScheme> scheme){
    timeScheme_  = scheme;
  }

  void BaseFeFunction::AddEntityList( shared_ptr<EntityList> list ){
    // Security check: If entity list was already added, leave
    
    // Note: As the shared_ptr to an Entitylist is not
    // unique within CFS, we have to ensure, that the names of 
    // the entity lists rather than the pointers match!
    
    std::string entName = list->GetName();
    bool found = false;
    for( UInt i = 0; i < entities_.GetSize(); ++i ) {
      if( entities_[i]->GetName() == entName )
        found = true;
    }
    if( found ) {
      // nothing to be done
    } else {
      entities_.Push_back(list);
      regions_.insert( list->GetRegion() );
      
    }
  }
  
  StdVector< shared_ptr<EntityList> > BaseFeFunction::GetEntityList(){
    return entities_;
  }
  
  const std::set<RegionIdType>& BaseFeFunction::GetRegions() const {
    return regions_;
  }
   
  void BaseFeFunction::SetFctId(FeFctIdType id ) {
    fctId_ = id;
  }
  
  //! Set the PDE pointer of the function
  void BaseFeFunction::SetPDE(SinglePDE* pde){
    pde_ = pde;
  }

  //! Get the PDE Pointer
  SinglePDE* BaseFeFunction::GetPDE(){
    return pde_;
  }

  void BaseFeFunction::SetGrid(Grid* grid){
    grid_ = grid;
  }

  Grid* BaseFeFunction::GetGrid(){
    return grid_;
  }

  void BaseFeFunction::SetSystem( AlgebraicSys* sys){
    algsys_ = sys;
  }

  AlgebraicSys* BaseFeFunction::GetSystem(){
    return algsys_;
  }

  void BaseFeFunction::AddHomDirichletBc( shared_ptr<HomDirichletBc> bc ){
    hdBcs_.Push_back(bc);
    entities_.Push_back(bc->entities);
  }

  void BaseFeFunction::AddInhomDirichletBc( shared_ptr<InhomDirichletBc> bc ){
    idBcs_.Push_back(bc);
    entities_.Push_back(bc->entities);
  }

  void BaseFeFunction::AddConstraint( shared_ptr<Constraint> bc ){
    constraints_.Push_back(bc);
    entities_.Push_back(bc->masterEntities);
    entities_.Push_back(bc->slaveEntities);
  }
  
  UInt BaseFeFunction::GetVecSize() const {
    assert( result_ );
    assert( dimType_ == CoefFunction::VECTOR );
    return result_->dofNames.GetSize();
  }
  
  void BaseFeFunction::GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION( "GetTensorSize not implemented");
  }
  
  std::string BaseFeFunction::ToString() const {
    std::string ret;
    ret += "FeFunction for result '" + 
        SolutionTypeEnum.ToString(result_->resultType) + "'";
    return ret;
  }
  
  
  // -----------------------------------
  BaseFeFunction::MapContext::MapContext() {
    algSys = NULL;
    assemble = NULL;
    sol = NULL;
  }
  
  BaseFeFunction::MapContext::~MapContext() {
    delete algSys;
    algSys = NULL;
    
    delete assemble;
    assemble = NULL;
    
    delete sol;
    sol = NULL;
    
    entityEqns.Clear();
    
  }
  
  // ========================================================================
  //  T E M P L A T I Z E D    V E R S I O N 
  // ========================================================================

  template<typename T>
  FeFunction<T>::FeFunction(){
    coeffs_ = NULL;
    factor_ = 1.0;
    timeDerivOrder_ = 0;
    idOp_ = NULL;
    isComplex_ = std::tr1::is_same<T,Complex>::value;
    MathParser * mp = domain->GetMathParser();
    
    // Add expression for calculating the time derivative in the
    // harmonic case
    if( IsComplex( )) {
      mp->SetExpr(mHandle_, "2*pi*f");
      mp->AddExpChangeCallBack(
          boost::bind(&FeFunction<T>::UpdateTimeDeriv, this ),
          mHandle_ );
    }
  }
  
  

  template<typename T>
  FeFunction<T>::~FeFunction(){
    if( domain) {
      MathParser * mp = domain->GetMathParser();
      mp->ReleaseHandle( mHandle_ );
    }
    if( idOp_ )
      delete idOp_;
    idOp_ = NULL;
    
    // only delete internal pointer, if no 
    // time-derivative was set
    if( timeDerivOrder_ == 0 && coeffs_)
      delete coeffs_;
    coeffs_ = NULL;
  }
  
  template<typename T>
  void FeFunction<T>::SetTimeDerivOrder( UInt i, 
                                         shared_ptr<FeFunction<T> >) {
    EXCEPTION("Only meaningful in harmonic case");
  }
  
  template<>
  void FeFunction<Complex>::SetTimeDerivOrder( UInt i, 
                                               shared_ptr<FeFunction<Complex> > feFct ) {
    if( i > 2 )  {
      EXCEPTION("Only time derivatives up to order 2 defined!");
    }
    timeDerivOrder_ = i;
    
    // delete old coefficient fector and 
    if( coeffs_ )
      delete coeffs_;
    coeffs_ = feFct->coeffs_;
    
    // not sure if other data members have to be copied as well ....
    
  }
  
  template<typename T>
  void FeFunction<T>::UpdateTimeDeriv() {
    EXCEPTION("Only meaningful in harmonic case");
  }
  
  
  template<>
    void FeFunction<Complex>::UpdateTimeDeriv() {
    MathParser * mp = domain->GetMathParser();
    
    Double omega = mp->Eval( mHandle_ );
    switch( timeDerivOrder_ ) {
      case 0:
        factor_ = Complex(1.0,0);
        break;
      case 1:
        factor_ = Complex(0,omega);
        break;
      case 2:
        factor_ = Complex(-(omega*omega),0);
        break;
    }
  }
  
  template<typename T>
  void FeFunction<T>::Finalize(){
    
    // assert that functionId was set
    static bool warn = true;
    if(warn) 
    {
      WARN("Add some more consistency checks here");
      warn = false;
    }
    
    if (fctId_ == NO_FCT_ID ) {
      EXCEPTION("No fctId was set!");
    }
    
    // only create new vector, if we are not a time
    // derivative fe function
    if( timeDerivOrder_ == 0 ) {
      coeffs_ = new Vector<T>(feSpace_->GetNumEquations());
    }
    
   
  }
  
  template<typename T>
  void FeFunction<T>::SetResultInfo( shared_ptr<ResultInfo> info ){
     result_ = info;

     // now initialize the members correctly
     switch( info->entryType ) {
       case ResultInfo::SCALAR:
         dimType_ = CoefFunction::SCALAR;
         break;
       case ResultInfo::VECTOR:
         dimType_ = CoefFunction::VECTOR;
         break;
       default:
         EXCEPTION("Only entry types SCALAR, VECTOR and TENSOR "
                    << "are handled" );
         break;
     }
     if( feSpace_ ) {
     if( feSpace_->GetSpaceType() == FeSpace::HCURL ) {
          dimType_ = CoefFunction::VECTOR;
        }
     }
     
     // Create interpolation operator
     UInt dim = grid_->GetDim();
     UInt numDofs = feSpace_->GetNumDofs();
     idOp_ = GenerateInterpolationOperator( dim, numDofs ); 
   }

  template<typename T>
  void FeFunction<T>::ExtractResult( shared_ptr<BaseResult> res ) {

    ResultInfo& resInfo = *(res->GetResultInfo() );
    UInt numDofs = resInfo.dofNames.GetSize();

    shared_ptr<EntityList> list = res->GetEntityList();
    Vector<T> & actSol = dynamic_cast<Result<T>&>(*res).GetVector();
    actSol.Resize( list->GetSize() * numDofs );

    EntityIterator it = list->GetIterator();
    actSol.Init();

    StdVector<Integer> eqnNums;
    UInt pos = 0;
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // get equation numbers
      feSpace_->GetEqns( eqnNums, it );

      // In case no equation was found, this indicates that the nodes, for which
      // the results should be calculated could not be found. Thus we have
      // to use interpolation to interpolate the continuous result to the
      // nodal locations in the entity list
      if( eqnNums.GetSize() == 0){
        //ok so the space does not know about this particular entity
        //we try to determine its value via interpolation
        Vector<Double> globCoord;
        Vector<T> elemSolution;
        Vector<T> dofSol;
        if(it.GetType()== EntityList::NODE_LIST){
          //now we obtain the global coords of the
          //node assuming that everything is the same grid. if not, we are in trouble anyway
          UInt curNodeNum = it.GetNode();
          grid_->GetNodeCoordinate(globCoord,curNodeNum,true);
        }else if(it.GetType() == EntityList::ELEM_LIST ||
                 it.GetType() == EntityList::SURF_ELEM_LIST){
          //determine global coord of element midpoint
          EXCEPTION("Interpolation for extract result not implemented for the Element case");
        }
        //Obtain intersecting element
        //const Elem* myElem  = grid_->GetElemAtGlobalCoord(globCoord,locCoord);
        
        
        // try to find the correct element, being one belonging to the regionlist of
        // this fefunction
        LocPoint lp;
        const Elem* myElem = 
            grid_->GetElemAtGlobalCoord(globCoord, lp, regions_ ); 
        if( !myElem ) {
          WARN("Some elements were skipped during the interpolation");
        }

        shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap( myElem, true );

        LocPointMapped lpm;
        lpm.Set(lp,esm);
        

        this->GetElemSolution(elemSolution,myElem);
        BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
        idOp_->ApplyOp(dofSol, lpm, ptFe, elemSolution );
        for(UInt iDim = 0; iDim < numDofs; iDim++ ) {
          actSol[pos++] = dofSol[iDim];

        }
      }else{
        Vector<T> & vals = *coeffs_;
        for ( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ){
          // check for homogeneous Dirichlet boundary condition
          if ( eqnNums[iDof] != 0 ) {
            actSol[pos++] = factor_ * vals[abs(eqnNums[iDof])-1];
          } else {
            actSol[pos++] = 0.0;
          }
        }
      }
    }
  }
  
  template<typename T>
  BaseBOperator* FeFunction<T>::GenerateInterpolationOperator(UInt dim, UInt dofDim){
    BaseBOperator* myOP = NULL;
    FeSpace::SpaceType curType = feSpace_->GetSpaceType();
    switch(curType){
      case FeSpace::H1:
      case FeSpace::L2:
        if(dim==2){
          if(dofDim==1)
            myOP = new IdentityOperator<FeH1,2,1,T>();
          else if(dofDim==2)
            myOP = new IdentityOperator<FeH1,2,2,T>();
          else if(dofDim==3)
            myOP = new IdentityOperator<FeH1,2,3,T>();
        }else{
          if(dofDim==1)
            myOP = new IdentityOperator<FeH1,3,1,T>();
          else if(dofDim==2)
            myOP = new IdentityOperator<FeH1,3,2,T>();
          else if(dofDim==3)
            myOP = new IdentityOperator<FeH1,3,3,T>();
        }
        break;
      case FeSpace::HCURL:
        if(dim==2){
          if(dofDim==1)
            myOP = new IdentityOperator<FeHCurl,2,1,T>();
        }else{
          if(dofDim==1)
            myOP = new IdentityOperator<FeHCurl,3,1,T>();
        }
        break;
      default:
        EXCEPTION("FeSpace type not suited for interpolation");
        break;
    }
    return myOP;
  }

  
  template<typename T> BiLinearForm* FeFunction<T>::
  GenerateInterpolBilinForm( UInt spaceDim, UInt dofDim ) {
    BiLinearForm *massInt = NULL;
    FeSpace::SpaceType curType = feSpace_->GetSpaceType();
    PtrCoefFct unity;
    if ( std::tr1::is_same<T,Complex>::value ) { 
      unity = CoefFunction::Generate(Global::COMPLEX, "1.0");
    } else {
      unity = CoefFunction::Generate(Global::REAL, "1.0");
    }
    switch(curType){
      case FeSpace::H1:
      case FeSpace::L2:
        if(spaceDim==1){
          // =============
          //  1D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<IdentityOperator<FeH1,1,1,T>,T,T>(unity, 1.0);
          } else if(dofDim==2) {
            massInt = new BBInt<IdentityOperator<FeH1,1,2,T>,T,T>(unity, 1.0);
          }else if(dofDim==3){
            massInt = new BBInt<IdentityOperator<FeH1,1,3,T>,T,T>(unity, 1.0);
          }
        } else if(spaceDim==2){
          // =============
          //  2D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<IdentityOperator<FeH1,2,1,T>,T,T>(unity, 1.0);
          } else if(dofDim==2) {
            massInt = new BBInt<IdentityOperator<FeH1,2,2,T>,T,T>(unity, 1.0);
          }else if(dofDim==3){
            massInt = new BBInt<IdentityOperator<FeH1,2,3,T>,T,T>(unity, 1.0);
          }
        } else if(spaceDim==3){
          // =============
          //  3D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<IdentityOperator<FeH1,3,1,T>,T,T>(unity, 1.0);
          } else if(dofDim==2) {
            massInt = new BBInt<IdentityOperator<FeH1,3,2,T>,T,T>(unity, 1.0);
          }else if(dofDim==3){
            massInt = new BBInt<IdentityOperator<FeH1,3,3,T>,T,T>(unity, 1.0);
          }
        }
        break;
      case FeSpace::HCURL:
        break;
      default:
        EXCEPTION("FeSpace type not suited for interpolation");
        break;
    }
    return massInt;
  }
  
  template<typename T> LinearForm* FeFunction<T>::
  GenerateInterpolLinForm( UInt spaceDim, UInt dofDim, PtrCoefFct coefFct ) {
    LinearForm * rhsInt = NULL;
    FeSpace::SpaceType curType = feSpace_->GetSpaceType();
    switch(curType){
      case FeSpace::H1:
      case FeSpace::L2:
        if(spaceDim==1){
          // =============
          //  1D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,1,T>,T>(1.0, coefFct);
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,2,T>,T>(1.0, coefFct);
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,3,T>,T>(1.0, coefFct);
          }
        } else if(spaceDim==2){
          // =============
          //  2D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,1,T>,T>(1.0, coefFct);
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,2,T>,T>(1.0, coefFct);
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,3,T>,T>(1.0, coefFct);
          }
        } else if(spaceDim==3){
          // =============
          //  3D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,1,T>,T>(1.0, coefFct);
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,2,T>,T>(1.0, coefFct);
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,3,T>,T>(1.0, coefFct);
          }
        }
        break;
      case FeSpace::HCURL:
        break;
      default:
        EXCEPTION("FeSpace type not suited for interpolation");
        break;
    }
    return rhsInt;
  }

  template<typename T>
  void FeFunction<T>::GetEntitySolution( SingleVector& elemSol, 
                                         const EntityIterator& it ){
    Vector<T> & temp = dynamic_cast<Vector<T>&>(elemSol);
    StdVector<Integer> eqns;
    Vector<T> & vals = *coeffs_;
    feSpace_->GetEqns(eqns, it);
    temp.Resize(eqns.GetSize());
    for(UInt iDof = 0 ; iDof < eqns.GetSize(); iDof++){
      if( eqns[iDof] != 0 ) {
        temp[iDof] = factor_ * vals[eqns[iDof]-1];
      } else {
        temp[iDof] = 0.0;
      }
    }
  }

  //! Get solution as matrix for specific entity
   template<typename T>
   void FeFunction<T>::GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                   const EntityIterator& it ){
     //for now we put the unkowns in the columns
     //and the dof entrys in rows
     Matrix<T> & temp = dynamic_cast<Matrix<T>&>(elemSol);
     UInt dofsPerUnknown = result_->dofNames.GetSize();
     StdVector<Integer> eqns;
     Vector<T> & vals = *coeffs_;
     temp.Resize(feSpace_->GetNumFunctions(it), dofsPerUnknown);
     for(UInt iDof = 0 ; iDof < dofsPerUnknown ; iDof++){
       feSpace_->GetEqns(eqns, it,iDof);
       for(UInt iEqn = 0;iEqn < eqns.GetSize() ; iEqn++){
         if( eqns[iEqn] > 0 ) {
         temp[iDof][iEqn] = factor_ * vals[eqns[iEqn]-1];
       } else {
         temp[iDof][iEqn] = 0.0;
       }
     }
   }
   }
   
  template<typename T>
  void FeFunction<T>::GetElemSolution( Vector<T>& elemSol,
                                         const Elem* elem ) {
    StdVector<Integer> eqns;
    Vector<T> & vals = *coeffs_;
    feSpace_->GetElemEqns(eqns, elem);
    elemSol.Resize(eqns.GetSize());
    for(UInt i= 0 ; i< eqns.GetSize(); i++){
      if( eqns[i] != 0 ) {
        elemSol[i] = factor_ * vals[std::abs(eqns[i])-1]; 
      } else {
        elemSol[i] = 0.0;
      }
    }
  }
  
  template<typename T>
  void FeFunction<T>::MapCoefFctToSpace(shared_ptr<EntityList> entityList, 
                                        shared_ptr<CoefFunction> coefFct,
                                        std::map <Integer, T>& vals,
                                        bool cache,
                                        const std::set<UInt>& comp ) {

    
    // Perform some checks at first:
    // a) if coefficient function is constant -> perform "easy mapping", without
    //    any caching. Here we simply take the constant value and 
    // b) if coefficient function has a general (spatial) dependency,
    //    we create solve FE problem on the subset of the FeSpace, as defined
    //    by the boundary condition, i.e. a mass-bilinearform is created and an
    //    according RHS integrator. The auxilliary data needed can be cached for
    //    repeated access / changing values (e.g. time / frequency dependend boundary+
    //    conditions.
    
    if (IS_LOG_ENABLED(fefunc, trace)) {
      StdVector<UInt> compVec;
      std::set<UInt>::const_iterator it = comp.begin();
      for (; it != comp.end(); ++it ) {
        compVec.Push_back(*it);
      }
      LOG_TRACE(fefunc) << "Mapping coeffct " << coefFct->ToString() 
                          << " on " << entityList->GetName() << " for dofs "
                          << compVec.ToString() << " for FeFunction "
                          << SolutionTypeEnum.ToString(result_->resultType);
    }
    
    
    if( coefFct->GetDependency() == CoefFunction::CONST ||
        coefFct->GetDependency() == CoefFunction::TIMEFREQ ) {
      // --------------------------
      //  SIMPLE MAPPING MECHANISM
      // --------------------------
      LocPointMapped lpm;
      StdVector<Integer> eqns, vEqns;
      Vector<T> dummyVec;
      if( coefFct->GetDimType() == CoefFunction::SCALAR ) {
        dummyVec.Resize(1);
        coefFct->GetScalar(dummyVec[0], lpm);
      } else {
        coefFct->GetVector( dummyVec, lpm);
      }
              
      // loop over all dofs
      std::set<UInt>::const_iterator it = comp.begin();
      for( ; it != comp.end(); ++it) {
      
        T val = dummyVec[*it];
        feSpace_->GetEntityListEqns( eqns, entityList, *it );
        UInt numEqns = eqns.GetSize();
        
        // switch depending on space type:
        if( !feSpace_->IsHierarchical() ) {
          // --- non-hierachical case ---
          for( UInt i = 0; i < numEqns; ++i ) {
            vals[eqns[i]] = val; 
          }
        } else {
          // -- hierachical case ---
          for( UInt i = 0; i < numEqns; ++i ) {
            vals[eqns[i]] = 0.0; 
          }
          
          // get nodal vertex nodes
          feSpace_->GetEntityListEqns(vEqns, entityList, *it, BaseFE::VERTEX );
          UInt numVEqns = vEqns.GetSize();
          for( UInt i = 0; i < numVEqns; ++i ) {
            vals[vEqns[i]] = val; 
          }
        }
        
      } // loop: dofs
    } else 
    {
      // ---------------------------
      //  GENERAL MAPPING MECHANISM
      // ---------------------------
      std::string name = entityList->GetName();
      MapContext* ctx = NULL; 
      
      // check, if caching is activated and if entry can be found in cache
      if( cache == true 
          && entityCtx_.find(name) != entityCtx_.end() ) {
        ctx = entityCtx_[name];
      } else { 
        // generate new context
        ctx = new MapContext();
        
        bool isComplex = std::tr1::is_same<T,Complex>::value;
        
        // generate new algebraic system
        ctx->olasNode.reset(new ParamNode());
        ctx->olasNode->SetName(std::string("IDBC-") + name);
        ctx->infoNode.reset(new ParamNode(ParamNode::INSERT));
        ctx->algSys = new AlgebraicSys( ctx->olasNode, ctx->infoNode, isComplex );
        
        // generate new assemble class    
        BasePDE::AnalysisType aType = 
            isComplex ? BasePDE::HARMONIC : BasePDE::STATIC;
        ctx->assemble = new Assemble( ctx->algSys, aType, 0 );
        
        // --------------------------------------------------------------------
        // generate (dimensional and space-dependent) interpolation (bi)linear
        // forms for the auxiliary problem
        // --------------------------------------------------------------------
        const Elem* el = entityList->GetIterator().GetElem();
        UInt dim = Elem::shapes[el->type].dim;
        UInt dofDim = feSpace_->GetNumDofs();
        BiLinearForm *massInt = GenerateInterpolBilinForm(dim, dofDim);
        LinearForm * rhsInt = GenerateInterpolLinForm(dim, dofDim, coefFct);

        BiLinFormContext * massCtx = new BiLinFormContext( massInt, STIFFNESS);
        massInt->SetName("Interpolator");
        massCtx->SetEntities( entityList,entityList );
        massCtx->SetFeFunctions(shared_from_this(), shared_from_this());
        ctx->assemble->AddBiLinearForm(massCtx);

        LinearFormContext * rhsCtx = new LinearFormContext( rhsInt );
        rhsInt->SetName("RhsInterpolator");
        rhsCtx->SetEntities( entityList );
        rhsCtx->SetFeFunction(shared_from_this());
        ctx->assemble->AddLinearForm(rhsCtx);

        // generate a mapping global eqnNr -> entityList local one for all dofs
        std::map<Integer,Integer> eqnMap;
        std::set<Integer> registeredEqns;
        
        // loop over all dofs
        if( comp.size() == 0 ) {
          // scalar problem
          feSpace_->GetEntityListEqns( ctx->entityEqns, entityList );
        } else {
          StdVector<Integer> tmp;
          std::set<UInt>::const_iterator it = comp.begin();
          for( ; it != comp.end(); it++ ) {
            feSpace_->GetEntityListEqns( tmp, entityList, *it );
            for( UInt k = 0; k < tmp.GetSize(); ++k ) {
              ctx->entityEqns.Push_back( tmp[k] ) ;
            }
          }
        }
        // fill map and pass it to assemble class
        UInt numEqns = ctx->entityEqns.GetSize();
        for( UInt i = 0; i < numEqns; ++i ) {
          eqnMap[ctx->entityEqns[i]] = i+1;
          registeredEqns.insert(i+1);
        }
        
        // ------------------------
        //  setup algebraic system
        // ------------------------
        // 1) setup matrix graph
        std::map<FeFctIdType, FeFctIdType> feFctIdMap;
        ctx->algSys->GraphSetupInit( 1,  false );
        FeFctIdType newFctId = ctx->algSys->ObtainFctId("interpolation");
        ctx->algSys->RegisterFct(newFctId, eqnMap.size(), eqnMap.size() );
        feFctIdMap[newFctId] = fctId_;
        ctx->assemble->SetEqnCustomMap(eqnMap, feFctIdMap);
        
        // define matrix graph and SBM blocks
        AlgebraicSys::SBMBlockDef sbmBlock;
        sbmBlock[newFctId] = registeredEqns;
        ctx->algSys->DefineSBMMatrixBlock( sbmBlock, false );
        ctx->algSys->FinishRegistration();
        ctx->assemble->SetupMatrixGraph(newFctId, newFctId);
        ctx->algSys->GraphSetupDone();
        ctx->algSys->CreateLinSys();
        ctx->algSys->InitMatrix();
        
        // create solver and preconditioner
        ctx->algSys->CreateSolver();
        ctx->algSys->CreatePrecond();

        // now reset AlgebraicSystem 
        ctx->algSys->InitRHS();
        ctx->algSys->InitSol();
        
        // assemble mapping matrix
        ctx->assemble->AssembleMatrices();

        // setup the preconditioner and solver
        ctx->algSys->SetupPrecond(ctx->infoNode);
        ctx->algSys->SetupSolver(ctx->infoNode);
        
        // initialize solution SBM vector
        ctx->sol = new SBM_Vector();
        ctx->sol->Resize(1);
        Vector<T> * tmp = new Vector<T>(eqnMap.size());
        ctx->sol->SetSubVector(tmp,0);
        ctx->sol->SetOwnership(true);
        
        
        // Store the context in the map
        if( cache ) {
          entityCtx_[name] = ctx;
        }
      } // end of first time setup of matrix

      // update the RHS due to the new coefficient vector
      ctx->algSys->InitRHS();
      ctx->assemble->AssembleLinRHS(NULL);

      
      // solve system and aquire solution
      ctx->algSys->Solve(ctx->infoNode);
      ctx->algSys->GetSolutionVal(*(ctx->sol));
      
      Vector <T> & sol = dynamic_cast<Vector<T> &>(*(ctx->sol->GetPointer(0)));
      
      // store result values to vals-map
      UInt numEqns = ctx->entityEqns.GetSize();
      
      for( UInt i = 0; i < numEqns; ++i ) {
        if( ctx->entityEqns[i] != 0 ) {
          vals[ctx->entityEqns[i]] = sol[i];
        }
      }
      
      // Print out information about the solution of the system
      // ctx->infoNode->ToXML(std::cerr);

    } // end of general mapping section
  }

 

  template<typename T>
  void FeFunction<T>::ApplyBC(){
    //===================================================
    // InHomogeneous BCs
    // ==================================================
    //loop over all inhomogeneous BCs
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
      InhomDirichletBc const & actBc = *(idBcs_[i]);
      
      // check, if entity list is defined on elements or nodes
      if( actBc.entities->GetType() == EntityList::ELEM_LIST ||
          actBc.entities->GetType() == EntityList::SURF_ELEM_LIST ) {

        // -------------------------
        // 1) Element Based Mapping
        // -------------------------
        
        // Map coefficient function onto the actual FeSpace
        std::map<Integer, T> coefs;
        this->MapCoefFctToSpace( actBc.entities, actBc.value, coefs, 
                                 true, actBc.dofs );

        // Loop over all entries and set them
        typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
        for( ; coefIt != coefs.end(); ++coefIt ) {
          Integer eqnNr = coefIt->first; 
          T val = coefIt->second;
          
          // In case of effective mass-formulation, 
          // the bcs have to be adjusted
          if( this->GetTimeScheme() ) {
            this->GetTimeScheme()->AdaptBC(val,val,0,eqnNr);
          }
          
          algsys_->SetDirichlet(  fctId_, eqnNr, val);
        }  // loop coefs 

      } else {
        // ------------------------------------------
        // 2) Legacy Nodal Based Coefficient Mapping
        // ------------------------------------------

        // Note: The legacy based is implemented only for 
        // coefficient functions not depending on space
        if( actBc.value->GetDependency() == CoefFunction::GENERAL ||
            actBc.value->GetDependency() == CoefFunction::SOLUTION) {
          EXCEPTION("Boundary condition, which are not defined on elements "
              << "are not allowed to be spatially dependent!");
        }
        
        // Evaluate coefficient function (general vector valued case)
        Vector<T> dummyVec;
        LocPointMapped lpm;
        StdVector<Integer> eqns;
        if( actBc.value->GetDimType() == CoefFunction::SCALAR ) {
          dummyVec.Resize(1);
          actBc.value->GetScalar(dummyVec[0], lpm);
        } else {
          actBc.value->GetVector( dummyVec, lpm);
        }
        // loop over all dofs
        std::set<UInt>::const_iterator it = actBc.dofs.begin();
        for( ; it != actBc.dofs.end(); ++it) {

          T val = dummyVec[*it];
          feSpace_->GetEntityListEqns( eqns,  actBc.entities, *it );
          UInt numEqns = eqns.GetSize();

          for( UInt i = 0; i < numEqns; ++i ) {
            // In case of effective mass-formulation, 
            // the bcs have to be adjusted
            if( this->GetTimeScheme() ) {
              this->GetTimeScheme()->AdaptBC(val,val,0,eqns[i]);
            }
            algsys_->SetDirichlet(  fctId_, eqns[i], val);
          } // loop: eqns
        } // loop: dofs
      } // if: not defined on elements
    } // loop over idbcs
  }
  
  template<typename T>
  void FeFunction<T>::GetVector(Vector<T>& vec, 
                                const LocPointMapped& lpm ){
    // get element solution
    Vector<T> elemSol;
    GetElemSolution( elemSol, lpm.ptEl);
    // apply identity operator to it
    BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
    idOp_->ApplyOp(vec, lpm, ptFe, elemSol );
  }
  
  
  template<typename T>
  void FeFunction<T>::GetScalar(T & scal, 
                                const LocPointMapped& lpm ){
    // get element solution
    Vector<T> elemSol, temp;
    GetElemSolution( elemSol, lpm.ptEl);
    // apply identity operator to it
    BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
    idOp_->ApplyOp(temp, lpm, ptFe, elemSol );
    scal = temp[0];
  }
  
  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class FeFunction<Double>;
  template class FeFunction<Complex>;
#endif
}
