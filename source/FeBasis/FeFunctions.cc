#include "FeFunctions.hh"

#include <boost/bind.hpp>
#include <boost/tr1/type_traits.hpp>
#include <def_expl_templ_inst.hh>

#include "PDE/SinglePDE.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "BaseFE.hh"
#include "H1/H1Elems.hh"
#include "HCurl/HCurlElems.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"

namespace CoupledField {
DECLARE_LOG(fefunc)
 DEFINE_LOG(fefunc, "feFunction")

 BaseFeFunction::BaseFeFunction(MathParser* mp){

  fctId_ = NO_FCT_ID;
  if(mp) {
    mp_ = mp;
    mHandle_ = mp_->GetNewHandle();
  } else  {
    mp_ = NULL;
  }
  algsys_ = NULL;

  // initialize members of coefficient function
  dependType_ = CoefFunction::GENERAL;
  isAnalytic_ = false;
  dimType_ = NO_DIM;

  }

  BaseFeFunction::~BaseFeFunction(){
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

  void BaseFeFunction::AddLoadCoefFunction( shared_ptr<CoefFunction> load ){
    this->loadCoefs_.Push_back(load);
    //entities for this case should already be definied....
    //entities_.Push_back(bc->entities);
  }

  void BaseFeFunction::AddExternalDataSource( shared_ptr<CoefFunction> coef ){
    this->externalDataCoefs_.Push_back(coef);
  }

  void BaseFeFunction::AddConstraint( shared_ptr<Constraint> bc ){
    constraints_.Push_back(bc);
    entities_.Push_back(bc->masterEntities);
    entities_.Push_back(bc->slaveEntities);
  }
  
  UInt BaseFeFunction::GetVecSize() const {
    assert( result_ ); assert( dimType_ == CoefFunction::VECTOR );
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
  
  
  
  // ========================================================================
  //  T E M P L A T I Z E D    V E R S I O N 
  // ========================================================================

  template<typename T>
  FeFunction<T>::FeFunction(MathParser* mp) :
     BaseFeFunction(mp)
  {
    coeffs_ = NULL;
    factor_ = 1.0;
    timeDerivOrder_ = 0;
    idOp_ = NULL;
    isComplex_ = std::tr1::is_same<T,Complex>::value;
    
    if( mp_ ) {
    // harmonic case
    if( IsComplex( )) {
      this->mp_->SetExpr(mHandle_, "2*pi*f");
      this->mp_->AddExpChangeCallBack(
          boost::bind(&FeFunction<T>::UpdateTimeDeriv, this ),
          mHandle_ );
      }
    }    
  }
  
  

  template<typename T>
  FeFunction<T>::~FeFunction(){
    if( domain) {
      if (mp_) {
        mp_->ReleaseHandle( mHandle_ );
      }
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
    if( !mp_ ) 
      return;
    Double omega = this->mp_->Eval( mHandle_ );
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
        Vector<T> elemSolution;
        Vector<T> dofSol;
        UInt nodeNum = 0;
        if(it.GetType()== EntityList::NODE_LIST){
          //now we obtain the global coords of the
          //node assuming that everything is the same grid. if not, we are in trouble anyway
          nodeNum = it.GetNode();
        }else if(it.GetType() == EntityList::ELEM_LIST ||
                 it.GetType() == EntityList::SURF_ELEM_LIST){
          //determine global coord of element midpoint
          EXCEPTION("Interpolation for extract result not implemented for the Element case");
        }
        // try to find the correct element, being one belonging to the regionlist of
        // this fefunction
        LocPoint lp;
        const Elem* myElem = 
            grid_->GetElemAtNode(nodeNum, lp, regions_ ); 
        if( !myElem ) {
          WARN("Some elements were skipped during the interpolation");
          for(UInt iDim = 0; iDim < numDofs; iDim++ ) {
            actSol[pos++] = 0.0;
          }
          continue;
        }

        shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap( myElem, true );

        LocPointMapped lpm;
        lpm.Set(lp,esm,0.0);
        

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
  GenerateInterpolBilinForm( UInt spaceDim, UInt dofDim, bool updatedGeo  ) {
    BiLinearForm *massInt = NULL;
    FeSpace::SpaceType curType = feSpace_->GetSpaceType();
    PtrCoefFct unity;
    if ( std::tr1::is_same<T,Complex>::value ) { 
      unity = CoefFunction::Generate(mp_, Global::COMPLEX, "1.0");
    } else {
      unity = CoefFunction::Generate(mp_, Global::REAL, "1.0");
    }
    switch(curType){
      case FeSpace::H1:
      case FeSpace::L2:
        if(spaceDim==1){
          // =============
          //  1D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,1,1,T>(),
                                   unity, 1.0, updatedGeo );
          } else if(dofDim==2) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,1,2,T>(), 
                                   unity, 1.0, updatedGeo );
          }else if(dofDim==3){
            massInt = new BBInt<T>(new IdentityOperator<FeH1,1,3,T>(), 
                                   unity, 1.0, updatedGeo );
          }
        } else if(spaceDim==2){
          // =============
          //  2D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,2,1,T>(), 
                                   unity, 1.0, updatedGeo );
          } else if(dofDim==2) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,2,2,T>(), 
                                   unity, 1.0, updatedGeo );
          }else if(dofDim==3){
            massInt = new BBInt<T>(new IdentityOperator<FeH1,2,3,T>(), 
                                   unity, 1.0, updatedGeo );
          }
        } else if(spaceDim==3){
          // =============
          //  3D Entities
          // =============
          if(dofDim==1) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,3,1,T>(), 
                                   unity, 1.0, updatedGeo );
          } else if(dofDim==2) {
            massInt = new BBInt<T>(new IdentityOperator<FeH1,3,2,T>(), 
                                   unity, 1.0, updatedGeo );
          }else if(dofDim==3){
            massInt = new BBInt<T>(new IdentityOperator<FeH1,3,3,T>(), 
                                   unity, 1.0, updatedGeo );
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
  GenerateInterpolLinForm( UInt spaceDim, UInt dofDim, PtrCoefFct coefFct, 
                           bool updatedGeo ) {
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
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,1,T>,T>(1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,2,T>,T>(1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,3,T>,T>(1.0, coefFct, updatedGeo );
          }
        } else if(spaceDim==2){
          // =============
          //  2D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,1,T>,T>(1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,2,T>,T>(1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,2,3,T>,T>(1.0, coefFct, updatedGeo );
          }
        } else if(spaceDim==3){
          // =============
          //  3D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,1,T>,T>(1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,2,T>,T>(1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<IdentityOperator<FeH1,3,3,T>,T>(1.0, coefFct, updatedGeo );
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
        temp[iDof] = factor_ * vals[std::abs(eqns[iDof])-1];
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
         temp[iDof][iEqn] = factor_ * vals[std::abs(eqns[iEqn])-1];
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
        feSpace_->MapCoefFctToSpace( actBc.entities, actBc.value, coefs, 
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
  void FeFunction<T>::ApplyLoads(){
    //loop over all loads
    for ( UInt i = 0; i < loadCoefs_.GetSize(); i++ ) {
      if(loadCoefs_[i]->IsConservative()){
        loadCoefs_[i]->MapConservative(this->feSpace_,*this->coeffs_);
      }else{
        //ok here we pass again the work to the space
        shared_ptr<CoefFunction> curFnc = loadCoefs_[i];
        if(curFnc->GetEntityList().GetSize() > 1){
          EXCEPTION("There are currently more than one entity lists defined within the function. This implementation can not yet cope with that.")
        }
        shared_ptr<EntityList> curEnt = curFnc->GetEntityList()[0];
        // check, if entity list is defined on elements or nodes
        if( curEnt->GetType() == EntityList::ELEM_LIST ||
            curEnt->GetType() == EntityList::SURF_ELEM_LIST ) {
          // Map coefficient function onto the actual FeSpace
          std::map<Integer, T> coefs;
          feSpace_->MapCoefFctToSpace( curEnt, curFnc, coefs,
                                       false );
          Vector<T> & myVals = *this->coeffs_;
          typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
          for( ; coefIt != coefs.end(); ++coefIt ) {
            Integer eqnNr = coefIt->first;
            T val = coefIt->second;
            myVals[eqnNr-1] += val;
          }
        }
      }
    }
  }

  template<typename T>
  void FeFunction<T>::ApplyExternalData(){
    //loop over all loads
    for ( UInt i = 0; i < externalDataCoefs_.GetSize(); i++ ) {
      shared_ptr<CoefFunction> curFnc = externalDataCoefs_[i];
      if(curFnc->GetEntityList().GetSize() > 1){
        EXCEPTION("There are currently more than one entity lists defined within the function. This implementation can not yet cope with that.")
      }
      shared_ptr<EntityList> curEnt = curFnc->GetEntityList()[0];
      // check, if entity list is defined on elements or nodes
      if( curEnt->GetType() == EntityList::ELEM_LIST ||
          curEnt->GetType() == EntityList::SURF_ELEM_LIST ) {
        // Map coefficient function onto the actual FeSpace
        std::map<Integer, T> coefs;
        feSpace_->MapCoefFctToSpace( curEnt, curFnc, coefs,
                                     false );
        Vector<T> & myVals = *this->coeffs_;
        typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
        for( ; coefIt != coefs.end(); ++coefIt ) {
          Integer eqnNr = coefIt->first;
          T val = coefIt->second;
          myVals[eqnNr-1] = val;
        }

      }else{
        if( curFnc->GetDependency() == CoefFunction::GENERAL ||
            curFnc->GetDependency() == CoefFunction::SOLUTION) {
          EXCEPTION("Boundary condition, which are not defined on elements "
              << "are not allowed to be spatially dependent!");
        }
        EXCEPTION("Node lists are not yet supported");
      }
    }
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
