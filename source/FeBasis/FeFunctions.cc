#include "FeFunctions.hh"

#include <boost/bind/bind.hpp>
#include "Utils/mathParser/mathParser.hh"
#include "Utils/tools.hh"
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
  DEFINE_LOG(fefunc, "feFunction")

  
  BaseFeFunction::BaseFeFunction(MathParser* mp) 
  {
    if(mp) {
      mp_ = mp;
      mHandle_ = mp_->GetNewHandle();
    } 

    // initialize members of coefficient function
    dependType_ = CoefFunction::GENERAL;
    isAnalytic_ = false;
    dimType_ = NO_DIM;
  }

  BaseFeFunction::~BaseFeFunction(){
    idBcs_.Clear();
  }
  
  string BaseFeFunction::ToString() {
    return SolutionTypeEnum.ToString(result_->resultType) + ": ";
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
    
    // Note: As the shared_ptr to an Entitylist is not unique within CFS, we have to ensure,
    // that the names of the entity lists rather than the pointers match!

    // This can be very expensive for periodic B.C. (Constraints) as there
    // are at least each 2-node pair added twice (@see BaseFeFunction::AddConstraint())
    // for that case only check via two_node_entries_cache_
    bool found = false;

    if(list->GetSize() == 2 && list->GetType() == EntityList::NODE_LIST)
    {
      assert(boost::dynamic_pointer_cast<NodeList>(list) != NULL); // we ask for it in if, so make a static_cast
      shared_ptr<NodeList> nl = boost::static_pointer_cast<NodeList>(list);

      assert(nl->GetNodes().GetSize() == 2);
      std::list<unsigned int>& list = two_node_entries_cache_[nl->GetNodes()[0]];
      if(find(list.begin(), list.end(), nl->GetNodes()[1]) != list.end())
        found = true;
      else
        list.push_back(nl->GetNodes()[1]); // not found? Then add to cache and below the entry will be added
    }
    else
    {
      const std::string& entName = list->GetName();
      for(unsigned int i = 0; !found && i < entities_.GetSize(); ++i) {
        if(entities_[i]->GetName() == entName) {
          found = true;
          break;
        }
      }
    }

    if(!found) {
      entities_.Push_back(list);
      regions_.insert(list->GetRegion());
    }
  }
  
  void BaseFeFunction::GetEntitySolution(SingleVector& elemSol, const Elem* elem)
  {
    ElemList elemList(grid_);
    elemList.SetElement(elem);
    const EntityIterator& it = elemList.GetIterator();

    GetEntitySolution(elemSol, it);
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
    LOG_DBG(fefunc) << ToString() << "AddHomDirichletBc()";
    hdBcs_.Push_back( bc );
    entities_.Push_back(bc->entities);
    if( bc->entities->GetDefineType() == EntityList::REGION) {
      regions_.insert( bc->entities->GetRegion());
    }
  }

  void BaseFeFunction::AddInhomDirichletBc( shared_ptr<InhomDirichletBc> bc ){
    LOG_DBG(fefunc) << ToString() << "AddInhomDirichletBc()";
    idBcs_.Push_back(bc);
    entities_.Push_back(bc->entities);
    if( bc->entities->GetDefineType() == EntityList::REGION) {
          regions_.insert( bc->entities->GetRegion());
        }
  }

  void BaseFeFunction::AddLoadCoefFunction( PtrCoefFct load,
                                            const StdVector<shared_ptr<EntityList> >& lists){
    this->loadCoefs_[load] = lists;

    //entities for this case should already be definied....
    //entities_.Push_back(bc->entities);
  }

  void BaseFeFunction::AddLoadCoefFunction( PtrCoefFct coef,
                                            shared_ptr<EntityList >& list) {
    this->loadCoefs_[coef].Push_back(list);

  }

  void BaseFeFunction::AddExternalDataSource( PtrCoefFct coef,
                                              const StdVector<shared_ptr<EntityList> >& lists){
    LOG_DBG(fefunc) << ToString() << "AddExternalDataSource()" << "size of lists is " << lists.GetSize();
    this->externalDataCoefs_[coef] = lists;
  }

  void BaseFeFunction::AddExternalDataSource( PtrCoefFct coef,
                                              shared_ptr<EntityList > list){
    this->externalDataCoefs_[coef].Push_back(list);
  }
  
  void BaseFeFunction::RemoveExternalDataSource() {
    this->externalDataCoefs_.clear();
  }
 
  void BaseFeFunction::AddConstraint( shared_ptr<Constraint> bc ){
    constraints_.Push_back(bc);
    entities_.Push_back(bc->masterEntities);
    entities_.Push_back(bc->slaveEntities);
  }
  
  bool BaseFeFunction::HasPeriodicBC() const  {
    for(unsigned int i = 0; i < constraints_.GetSize(); i++)
      if(constraints_[i]->periodic)
        return true;
    return false;
  }

  bool BaseFeFunction::HasConstraint(std::string& name, unsigned int dof) const {
    for(unsigned int i = 0; i < constraints_.GetSize(); i++) {
      if(constraints_[i]->name == name && constraints_[i]->masterDof == dof)
        return true;
    }
    return false;
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
    ret += "FeFunction for result '" + SolutionTypeEnum.ToString(result_->resultType) + "'";
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
    isComplex_ = std::is_same<T,Complex>::value;
    
    
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

    // Note: Currently we can not release the handle in a safe way,
    // as Fefunctions are only passed as shared_ptr, so we can not
    // know exactly, when the last FeFunction is freed. If this is within
    // a child Domain, the mathParser might already be deleted when
    // we try to attempt this.
    
    
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
    // FIXME
    static bool warn = false;
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
      delete coeffs_;
      coeffs_ = new Vector<T>(feSpace_->GetNumEquations());
      coeffs_->Init();
    }

    /* Check: If boundary conditions are defined on node lists with more than
     * one node and the space has no grid mapping, we issue a warning, as in
     * this case only the vertex-associated dofs would be fixed.
     */
    if( this->feSpace_->GetMapType(ALL_REGIONS) != FeSpace::GRID ) {

      std::string nodeListNames;
      shared_ptr<EntityList> list;
      {
        // Loop over all HDBCs
        HdBcList::iterator hdbcIt = hdBcs_.Begin();
        for(; hdbcIt != hdBcs_.End(); ++hdbcIt ) {
          list = (*hdbcIt)->entities;
          if( list->GetType() == EntityList::NODE_LIST &&
              list->GetSize() > 1 ) {
            nodeListNames += "\t" + list->GetName() + "\n";
          }
        }
      }

      {
        // Loop over all IDBCs
        IdBcList::iterator idbcIt = idBcs_.Begin();
        for(; idbcIt != idBcs_.End(); ++idbcIt ) {
          list = (*idbcIt)->entities;
          if( list->GetType() == EntityList::NODE_LIST &&
              list->GetSize() > 1 ) {
            nodeListNames += "\t" + list->GetName() + "\n";
          }
        }
      }

      if( !nodeListNames.empty()  ) {
        WARN( "In case of general / higher order approximation, boundary "
            << "conditions should be applied on (surface) elements instead "
            << "of node lists. For the quantity '" 
            << SolutionTypeEnum.ToString(result_->resultType)
            << "' the following node lists are used in boundary "
            << "conditions:\n\n" << nodeListNames 
            << "\nPlease consider changing them to (surface) element lists "
            << "or use Lagrangian polynomials with grid order!")
      }
    }
  }
  
  template<typename T>
  void FeFunction<T>::CleanUp(){
    if( domain) {
      if (mp_) {
        mp_->ReleaseHandle( mHandle_ );
        mp_ = NULL;
      }
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
     if(feSpace_ && feSpace_->GetSpaceType() == FeSpace::HCURL)
         dimType_ = CoefFunction::VECTOR;

     // Create interpolation operator
     UInt dim = grid_->GetDim();
     UInt numDofs = feSpace_->GetNumDofs();
     idOp_ = GenerateInterpolationOperator( dim, numDofs ); 
   }

  template<typename T>
  void FeFunction<T>::ExtractResult(shared_ptr<BaseResult> res) 
  {
    ResultInfo* resInfo = res->GetResultInfo().get();
    UInt numDofs = resInfo->dofNames.GetSize();

    shared_ptr<EntityList> list = res->GetEntityList();
    assert(dynamic_cast<Result<T>*>(res.get()) != nullptr);
    Vector<T>& actSol = dynamic_cast<Result<T>&>(*res).GetVector();
    actSol.Resize(list->GetSize() * numDofs);
    LOG_DBG(fefunc) << ToString() << "ExtractResult for size " << actSol.GetSize();

    EntityIterator it = list->GetIterator();
    actSol.Init();

    StdVector<Integer> eqnNums;
    UInt pos = 0;
    for(it.Begin(); !it.IsEnd(); it++) {
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
        }else if(it.GetType() == EntityList::ELEM_LIST || it.GetType() == EntityList::SURF_ELEM_LIST){
          //determine global coord of element midpoint
          throw Exception("Interpolation for extract result not implemented for the Element case");
        }
        // try to find the correct element, being one belonging to the regionlist of
        // this fefunction
        LocPoint lp;
        const Elem* myElem = grid_->GetElemAtNode(nodeNum, lp, regions_ );
        
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
        // the equations case
        assert(coeffs_ != nullptr && coeffs_->GetSize() > 0);
        Vector<T>& vals = *coeffs_;
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
        // currently the HCURL space exists only in 3D
        myOP = new IdentityOperator<FeHCurl,3,1,T>();
        break;
      case FeSpace::CONSTANT:
        // this space does not approximate space
        WARN("The interpolation operator for FeSpaceConst is not initialized.");
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
    if ( std::is_same<T,Complex>::value ) { 
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
        if( spaceDim == 3 ) {
          // =============
          //  3D Entities
          // =============
            massInt = new BBInt<T>(new IdentityOperator<FeHCurl,3,1,T>(), 
                                   unity, 1.0, updatedGeo );
        } else {
          EXCEPTION("HCURL mapping only working in 3D")
        }
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
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,1,T>(), 1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,2,T>(), 1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,3,T>(), 1.0, coefFct, updatedGeo );
          }
        } else if(spaceDim==2){
          // =============
          //  2D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,1,T>, 1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,2,T>, 1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,2,3,T>, 1.0, coefFct, updatedGeo );
          }
        } else if(spaceDim==3){
          // =============
          //  3D Entities
          // =============
          if(dofDim==1) {
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,3,1,T>, 1.0, coefFct, updatedGeo );
          } else if(dofDim==2) {
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,3,2,T>, 1.0, coefFct, updatedGeo );
          }else if(dofDim==3){
            rhsInt = new BUIntegrator<T>( new IdentityOperator<FeH1,3,3,T>, 1.0, coefFct, updatedGeo );
          }
        }
        break;
        
      case FeSpace::HCURL:
        if( spaceDim == 3 ) {
          // =============
          //  3D Entities
          // =============
          rhsInt = new BUIntegrator<T>(new IdentityOperator<FeHCurl,3,1,T>(), 1.0, coefFct, updatedGeo );
        } else {
          EXCEPTION("HCURL mapping only working in 3D")
        }

        break;
        
      default:
        EXCEPTION("FeSpace type not suited for interpolation");
        break;
    }
    return rhsInt;
  }

  template<typename T>
  void FeFunction<T>::GetEntitySolution( SingleVector& elemSol, const EntityIterator& it )
  {
    LOG_DBG(fefunc) << ToString() << "FE:GES: #es=" << elemSol.GetSize() << " it=" << it.ToString();
    Vector<T>& temp = dynamic_cast<Vector<T>&>(elemSol);

    StdVector<Integer> eqns;
    Vector<T>& vals = *coeffs_;
    feSpace_->GetEqns(eqns, it);

    temp.Resize(eqns.GetSize());
    // In case no equation was found, this indicates that the nodes, for which
    // the results should be calculated could not be found. Thus we have
    // to use interpolation to interpolate the continuous result to the
    // nodal locations in the entity list
    if( eqns.GetSize() == 0){

      //ok so the space does not know about this particular entity
      //we try to determine its value via interpolation
      Vector<T> elemSolution;
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
      const Elem* myElem = grid_->GetElemAtNode(nodeNum, lp, regions_);
      if( !myElem ) {
        WARN("Some elements were skipped during the interpolation");
        temp.Init();
      }

      shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap( myElem, true );

      LocPointMapped lpm;
      lpm.Set(lp,esm,0.0);

      this->GetElemSolution(elemSolution,myElem);
      BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);

      idOp_->ApplyOp(temp, lpm, ptFe, elemSolution );
    } else {
      LOG_DBG3(fefunc) << "vals=" << vals.ToString() << "  eqns=" << eqns.ToString();
      for(UInt iDof = 0 ; iDof < eqns.GetSize(); iDof++){
        if( eqns[iDof] != 0 ) {
          temp[iDof] = factor_ * vals[std::abs(eqns[iDof])-1];
        } else {
          temp[iDof] = 0.0;
        }
      }
    }
  }

  //! Get solution as matrix for specific entity
   template<typename T>
   void FeFunction<T>::GetEntitySolutionAsMatrix( DenseMatrix& elemSol, const EntityIterator& it ){
     //for now we put the unkowns in the columns
     //and the dof entrys in rows
     Matrix<T>& temp = dynamic_cast<Matrix<T>&>(elemSol);
     UInt dofsPerUnknown = result_->dofNames.GetSize();
     StdVector<Integer> eqns;
     Vector<T>& vals = *coeffs_;
     temp.Resize(feSpace_->GetNumFunctions(it), dofsPerUnknown);
     for(UInt iDof = 0 ; iDof < dofsPerUnknown ; iDof++) {
       feSpace_->GetEqns(eqns, it,iDof);
       for(UInt iEqn = 0; iEqn < eqns.GetSize() ; iEqn++) {
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
    LOG_DBG(fefunc) << ToString() << "GetElemSolution()";

    // NOTE: Element solution caching was disabled because it doesn't properly
    // track when the solution vector changes during nonlinear iterations.
    // This caused stale values to be returned, leading to incorrect results.
    // TODO: If caching is needed for performance, add a solution version counter
    // or invalidation mechanism.

    // Use thread-local work buffer to avoid per-call allocations
    StdVector<Integer>& eqns = work_eqns_.Mine();
    const Vector<T>& vals = *coeffs_;
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
    LOG_DBG(fefunc) << ToString() << "FE:ApplyBC() (inhomogeneous)";

    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
      InhomDirichletBc const& actBc = *(idBcs_[i]);
      
      // get time derivation order of boundary condition
      // (e.g. 2 for mech acceleration, 0 for mech displacement
      UInt bcOrder = actBc.timeDerivOrder;
            
      // check, if entity list is defined on elements or nodes
      if( actBc.entities->GetType() == EntityList::ELEM_LIST ||
          actBc.entities->GetType() == EntityList::SURF_ELEM_LIST ) {

          // -------------------------
          // 1) Element Based Mapping
          // -------------------------

          // Map coefficient function onto the actual FeSpace
          std::map<Integer, T> coefs;
          StdVector<shared_ptr<EntityList> > list(1);
          list[0] = actBc.entities;
          feSpace_->MapCoefFctToSpace( list, actBc.value, shared_from_this(), coefs,
                                        true, actBc.dofs, actBc.updatedGeo );

          Complex harm;
          UInt harmInt;
          LocPointMapped lpm;
          
          // Loop over all entries and set them
          typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
          for( ; coefIt != coefs.end(); ++coefIt ) {
            Integer eqnNr = coefIt->first;
            T val = coefIt->second;

          // if the solution order and the bc order do not match, adaptBC
          // (e.g. in mechanics solution order = 2 (mass formulation) but bc order = 0 (mech displ) or
          //  solution order = 0 (stiff formulation) but bc order = 2 (mech acc.) )
          if( this->GetTimeScheme() ) {
              this->GetTimeScheme()->AdaptBC(val,val,bcOrder,eqnNr);
          }
          /*
            // In case of effective mass-formulation,
            // the bcs have to be adjusted
            if( this->GetTimeScheme() ) {
              this->GetTimeScheme()->AdaptBC(val,val,0,eqnNr);
            }

              */

          if(!algsys_->IsMultHarm()){
            algsys_->SetDirichlet(  fctId_, eqnNr, val);
          } else {
            actBc.harm->GetScalar(harm, lpm);
            harmInt = harm.real();
            algsys_->SetDirichletMH(  fctId_, eqnNr, val, harmInt);
          }
        }  // loop coefs 

      } else {
        // ------------------------------------------
        // 2) Legacy Nodal Based Coefficient Mapping
        // ------------------------------------------

        // we can handle spatial coef, e.g. expressions
        if(actBc.value->GetDependency() == CoefFunction::SOLUTION)
          EXCEPTION("Boundary condition, which are not defined on elements are not allowed to be solution dependent!");
        
        // Evaluate coefficient function (general vector valued case)
        Complex harm;
        unsigned int harmInt = 0;
        LocPointMapped lpm; // usually blind dummy, for spacial we set lpm.lp.number with the element number
        StdVector<int> eqns;

        // do we have to handle spacial dependency? We don't keepp the list, so no shared pointer
        NodeList* nl = actBc.value->IsSpacialDependent() ? dynamic_cast<NodeList*>(actBc.entities.get()) : NULL;
        // we obtain a realistic dummy value, e.g. from an expression coef, other would ignore it
        if(nl)
          lpm.lp.number = nl->GetNodes()[0];

        // the dummyVec is mostly used for non-spacial values
        Vector<T> dummyVec(actBc.dofs.size());
        if(actBc.value->GetDimType() == CoefFunction::SCALAR)
          actBc.value->GetScalar(dummyVec[0], lpm);
        else
          actBc.value->GetVector(dummyVec, lpm);

        if(algsys_->IsMultHarm()){
          actBc.harm->GetScalar(harm, lpm);
          harmInt = harm.real();
        }
        LOG_DBG(fefunc) << "ABC: bc dofs=" << ToStringCont(actBc.dofs);
        LOG_DBG(fefunc) << "ABC: dummyVec=" << dummyVec.ToString();
        LOG_DBG(fefunc) << "ABC: bc entities=" << actBc.entities->ToString();
        LOG_DBG(fefunc) << "ABC: nl nodes=" << (nl != NULL ? nl->GetNodes().ToString() : "no node list in usage");

        // loop over all dofs
        for(std::set<UInt>::const_iterator it = actBc.dofs.begin(); it != actBc.dofs.end(); ++it)
        {
          // transform the node list to an equations list. Do cfs -M to validate
          feSpace_->GetEntityListEqns( eqns,  actBc.entities, *it );

          // in the spatial case we need to overwrite this value but save the effort for the normal case
          T val = dummyVec[*it];
          LOG_DBG(fefunc) << "ABC dof=" << *it << " val=" << val;

          //LOG_DBG3(fefunc) << "ABC: eqns=" << eqns.ToString();
          for( UInt i = 0; i < eqns.GetSize(); ++i )
          {
            // handle spatial case
            if(nl)
            {
              lpm.lp.number = nl->GetNodes()[i];
              if(actBc.value->GetDimType() == CoefFunction::SCALAR)
                actBc.value->GetScalar(dummyVec[0], lpm);
               else
                actBc.value->GetVector(dummyVec, lpm); // done multiple times for all dofs :(
              val = dummyVec[*it];
              LOG_DBG3(fefunc) << "ABC dof=" << *it << " i=" << i << " node=" << nl->GetNodes()[i] << " eqn=" << eqns[i] << " val=" << val;
            }

            // if the solution order and the bc order do not match, adaptBC
            // (e.g. in mechanics solution order = 2 (mass formulation) but bc order = 0 (mech displ) or
            //  solution order = 0 (stiff formulation) but bc order = 2 (mech acc.) )
            if( this->GetTimeScheme() )
              this->GetTimeScheme()->AdaptBC(val,val,bcOrder,eqns[i]);
            /*
            // In case of effective mass-formulation, 
            // the bcs have to be adjusted
            if( this->GetTimeScheme() )
              this->GetTimeScheme()->AdaptBC(val,val,0,eqns[i]);
            */
            if(!algsys_->IsMultHarm())
              algsys_->SetDirichlet(  fctId_, eqns[i], val);
            else
              algsys_->SetDirichletMH(  fctId_, eqns[i], val, harmInt);
          } // loop: eqns
        } // loop: dofs
      } // if: not defined on elements
    } // loop over idbcs
  }
  
  template<typename T>
  void FeFunction<T>::ApplyLoads(){
	  //LOG_DBG(fefunc) << ToString() << "ApplyLoads()";
	  //loop over all loads
	  LoadCoefList::iterator it = loadCoefs_.begin();
	  // Loop over all coeffunctions
	  for ( ; it != loadCoefs_.end(); ++it  ) {
		  PtrCoefFct ptCoef = it->first;
		  StdVector<shared_ptr<EntityList> >& lists = it->second;

		  CoefFunction::CoefInverseType type = ptCoef->GetInverseType();
		  CoefFunction::CoefInverseSourceApprox typeApprox = ptCoef->GetInverseSourceApproxType();
		  if(ptCoef->IsConservative()){
			  //this is a little circumfencial allocaing and releasing memory
			  //in each step. perhaps it would be better to make a class variable or do it
			  //differently somehow
			  Vector<T> loadVec(this->coeffs_->GetSize());
			  loadVec.Init();
			  ptCoef->MapConservative(this->feSpace_,loadVec);
			  //this is a hack for source localization
			  if ( type != CoefFunction::INVSOURCE || typeApprox == CoefFunction::DELTA)
				  this->algsys_->SetFncRHS(loadVec,this->fctId_);
			  }
			  else{
				  //ok here we pass again the work to the space

				  // Map coefficient function onto the actual FeSpace
				  std::map<Integer, T> coefs;
				  feSpace_->MapCoefFctToSpace( lists, ptCoef, shared_from_this(), coefs, false );

				  typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
				  for( ; coefIt != coefs.end(); ++coefIt ) {
					  this->algsys_->SetNodeRHS(coefIt->second,this->fctId_,(Integer)coefIt->first);
				  }
			  }

		  }// loop: coefs
  }

  template<typename T>
  void FeFunction<T>::ApplyLoads(PtrCoefFct& values) {
	  //loop over all loads
	  LoadCoefList::iterator it = loadCoefs_.begin();

	  // Loop over all coeffunctions
	  for ( ; it != loadCoefs_.end(); ++it  ) {
		  PtrCoefFct ptCoef = it->first;
		  StdVector<shared_ptr<EntityList> >& lists = it->second;

		  std::map<Integer, T> coefs;
		  feSpace_->MapCoefFctToSpace( lists, ptCoef, shared_from_this(), coefs, false );

		  typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
		  for( ; coefIt != coefs.end(); ++coefIt ) {
			  //std::cout << "Val: " << coefIt->second << std::endl;
			  this->algsys_->SetNodeRHS(coefIt->second,this->fctId_,(Integer)coefIt->first);
		  }
	  }// loop: coefs
  }


  template<typename T>
  void FeFunction<T>::ApplyGeomInfo(){

    boost::unordered_map< Integer, EqNodeGeom> eqIndGeomMap;
    UInt maxEqn, dim;
    feSpace_->CreateEquIndGeomMap(eqIndGeomMap, maxEqn, dim);
    //this EqNodeGeom could be used later on, but now we only need
    //index<->coordinate, so extract the coordinate from the map, eqnNum remains the key
    StdVector< Vector<Double> > indGeom;
    bool edge = false;
    if( feSpace_->GetSpaceType() != feSpace_->HCURL ){
      indGeom.Resize(eqIndGeomMap.size());
    }else{
      edge = true;
      indGeom.Resize(2 * eqIndGeomMap.size() + 1);
    }


    boost::unordered_map< Integer , EqNodeGeom >::const_iterator eqIt = eqIndGeomMap.begin();
    if( edge != true){
      while(eqIt != eqIndGeomMap.end() ){
        //used for mech and poisson
        indGeom[eqIt->second.indexNum - 1] = eqIt->second.coord;
       eqIt++;
      }
      this->algsys_->SetGeomIndexMap(indGeom, dim);
    }else{

      boost::unordered_map< Integer, StdVector<Integer> > eNodes;
      boost::unordered_map< Integer, Double > lengths;
      while(eqIt != eqIndGeomMap.end() ){
        eNodes[eqIt->second.indexNum - 1] = eqIt->second.eNodes;
        Vector<Double> n1, n2;
        n1 = eqIt->second.eCoords[0];
        n2 = eqIt->second.eCoords[1];
        lengths[eqIt->second.indexNum - 1] = n1.NormL2(n2);
        eqIt++;
      }
      this->algsys_->SetEdgeIndexMap(lengths, eNodes);
    }
  }

  template<typename T>
  void FeFunction<T>::ApplyExternalData(){

    LOG_DBG(fefunc) << ToString() << "ApplyExternalData()";
    LoadCoefList::iterator it = externalDataCoefs_.begin();
    // Loop over all loads
    for ( ; it != externalDataCoefs_.end(); ++it  ) {
      PtrCoefFct curFnc = it->first;
      StdVector<shared_ptr<EntityList> >& lists = it->second;

      // Map coefficient function onto the actual FeSpace
      std::map<Integer, T> coefs;
      feSpace_->MapCoefFctToSpace( lists, curFnc, shared_from_this(), coefs, false );
      Vector<T>& myVals = *this->coeffs_;
      typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
      for( ; coefIt != coefs.end(); ++coefIt ) {
        Integer eqnNr = coefIt->first;
        T val = coefIt->second;
        myVals[eqnNr-1] = val;
      }
    }
  }
  
  template<typename T>
  void FeFunction<T>::InitFromFeFunction( shared_ptr<BaseFeFunction> feFct ) {
   
    LOG_DBG(fefunc) << ToString() << "Initialize from other FeFunction";
    
    // Get related spaces
    shared_ptr<FeSpace> otherSpace = feFct->GetFeSpace();
    bool sameApprox = true;
    
    
    // Get common entity lists
    //FeFunction<T> & otherFct =
    //      dynamic_cast<FeFunction<T>& >(*feFct);
    //feFct->entities_
    //

    bool convertDoubleToComplex = (std::is_same<T,Complex>::value && !feFct->IsComplex());

    StdVector<shared_ptr<EntityList> > intersect;
    EntityList::Intersect( this->entities_, feFct->GetEntityList(), intersect );

    // Loop over all entities
    LOG_DBG3(fefunc) << ToString() << "Entities to be transferred: ";
    for(UInt iList = 0; iList < intersect.GetSize(); ++iList ) {
      // Check if all entities have the same approximation
      bool isSame = this->feSpace_->IsSameEntityApproximation(intersect[iList],otherSpace);
      LOG_DBG3(fefunc) << ToString() << "\t" << intersect[iList]->GetName() << " (same approximation: " << (isSame ? "true" : "false") << ")";
      sameApprox &= isSame;
    }
    
    if( !sameApprox ) {
      LOG_DBG(fefunc)<< "=> Applying General Mapping Algorithm <= ";
      // -----------------
      //  General Mapping
      // -----------------
      // Map coefficient function onto the actual FeSpace using the general mechanism
      std::map<Integer, T> coefs;
      feSpace_->MapCoefFctToSpace( intersect, feFct, shared_from_this(), coefs, false );
      Vector<T>& myVals = *this->coeffs_;
      typename std::map<Integer, T>::const_iterator coefIt = coefs.begin();
      for( ; coefIt != coefs.end(); ++coefIt ) {
        Integer eqnNr = coefIt->first;
        T val = coefIt->second;
        myVals[eqnNr-1] = val;
      }

    } else {
      // ---------------
      //  Copy entries
      // ---------------
      LOG_DBG(fefunc)<< "=> Applying Copy Based Algorithm <= ";
      Vector<T>& myVec = *(this->coeffs_);
      myVec.Init();
      if(convertDoubleToComplex){
        FeFunction<Double>& otherFct =
              dynamic_cast<FeFunction<Double>& >(*feFct);
        for( UInt iList = 0; iList < intersect.GetSize(); ++iList ) {
          shared_ptr<EntityList> actList = intersect[iList];

          EntityIterator it = actList->GetIterator();
          StdVector<Integer> eqns;
          // loop over all elements
          for( it.Begin(); !it.IsEnd(); it++ ) {
            Vector<Double> elemSol;
            otherFct.GetEntitySolution( elemSol, it);
            feSpace_->GetEqns( eqns, it);
            UInt numEqns = eqns.GetSize();
            // Loop over all equations
            for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
              if( eqns[iEqn] > 0 ) {
                myVec[eqns[iEqn]-1] = elemSol[iEqn];
              }
            } // loop: eqns
          } // loop: elements
        } // loop: lists
      }else{
        FeFunction<T>& otherFct = dynamic_cast<FeFunction<T>& >(*feFct);
        for( UInt iList = 0; iList < intersect.GetSize(); ++iList ) {
          shared_ptr<EntityList> actList = intersect[iList];

          EntityIterator it = actList->GetIterator();
          StdVector<Integer> eqns;
          // loop over all elements
          for( it.Begin(); !it.IsEnd(); it++ ) {
            Vector<T> elemSol;
            otherFct.GetEntitySolution( elemSol, it);
            feSpace_->GetEqns( eqns, it);
            UInt numEqns = eqns.GetSize();
            // Loop over all equations
            for( UInt iEqn = 0; iEqn < numEqns; ++iEqn ) {
              if( eqns[iEqn] > 0 ) {
                myVec[eqns[iEqn]-1] = elemSol[iEqn];
              }
            } // loop: eqns
          } // loop: elements
        } // loop: lists
      }
    }
    LOG_DBG(fefunc) << ToString() << "Finished initialization from other FeFunction";
  }

  template<typename T>
  void FeFunction<T>::GetVector(Vector<T>& vec,
                                const LocPointMapped& lpm ){
    Vector<T> elemSol;
    GetElemSolution( elemSol, lpm.ptEl);
    // apply identity operator to it
    BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
    idOp_->ApplyOp(vec, lpm, ptFe, elemSol );
    LOG_DBG(fefunc) << ToString() << ": Requesting vector solution of point in elem " << lpm.ptEl->elemNum << ". Sol[0] = " << vec[0];
  }

  template<typename T>
  void FeFunction<T>::GetAvgElemValue(T& vec, const Elem* elem) {
    Vector<T> elemSol;
    GetElemSolution( elemSol, elem);
    // average
    vec = 0.;
    for (UInt i = 0; i < elemSol.GetSize(); i++) {
      vec += elemSol[i];
    }
    vec /= elemSol.GetSize();
    LOG_DBG(fefunc) << ToString() << ": Requesting average solution of element " << elem->elemNum << ". Sol = " << vec;
  }
  
  template<typename T>
  void FeFunction<T>::GetScalar(T& scal, const LocPointMapped& lpm ){
    Vector<T> elemSol, temp;
    GetElemSolution( elemSol, lpm.ptEl);
    // apply identity operator to it
    if( feSpace_->GetSpaceType() == FeSpace::SpaceType::CONSTANT){
      temp = elemSol;
    }else{
      BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
      idOp_->ApplyOp(temp, lpm, ptFe, elemSol );
    }
    scal = temp[0];
  }
  
// Explicit template instantiation
template class FeFunction<Double>;
template class FeFunction<Complex>;

} // end of namespace
