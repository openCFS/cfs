#include "FeFunctions.hh"

#include <boost/bind.hpp>

#include "PDE/SinglePDE.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include <boost/tr1/type_traits.hpp>

#include "BaseFE.hh"
#include "H1/H1Elems.hh"
#include "HCurl/HCurlElems.hh"

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
      WARN("Entitylist with name '" << entName 
           << "' was already added" );
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
  
  // ========================================================================
  //  T E M P L A T I Z E D    V E R S I O N 
  // ========================================================================

  template<typename T>
  FeFunction<T>::FeFunction(){
    coeffs_ = NULL;
    factor_ = 1.0;
    timeDerivOrder_ = 0;
    idOp_ = NULL;
    
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

    UInt dim = grid_->GetDim();
    scoped_ptr<BaseBOperator<T> > interOp ( GenerateInterpolationOperator(dim,numDofs) );
    
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
#ifdef USE_INTERPOLATION
        Vector<Double> globCoord;
        StdVector<Vector<Double> > locCoord;
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
        StdVector<const Elem*> elems = grid_->GetElemsAtGlobalCoord(globCoord,locCoord);
        LocPoint myLp;
        const Elem* myElem = NULL;
        for ( UInt i = 0; i < elems.GetSize(); ++i ) {
          if( this->regions_.find( elems[i]->regionId ) != this->regions_.end() ) {
            myElem = elems[i];
            myLp = locCoord[i];
            break;
          }
        }
        
        if( !myElem ) {
          WARN("Some elements were skipped during the interpolation");
        }

        shared_ptr<ElemShapeMap> esm = grid_->GetElemShapeMap( myElem, true );

        LocPointMapped lpm;
        lpm.Set(myLp,esm);
        

        this->GetElemSolution(elemSolution,myElem);
        BaseFE * ptFe = feSpace_->GetFe(lpm.ptEl->elemNum);
        interOp->ApplyOp(dofSol, lpm, ptFe, elemSolution );
        for(UInt iDim = 0; iDim < numDofs; iDim++ ) {
          actSol[pos++] = dofSol[iDim];

        }
#else
        WARN("Interpolation not enabled but needed to extract this result. enable it plz!");
#endif
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
  BaseBOperator<T>* FeFunction<T>::GenerateInterpolationOperator(UInt dim, UInt dofDim){
    BaseBOperator<T>* myOP = NULL;
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
  void FeFunction<T>::ApplyBC(){
    //basic algo idea for nodes;
    //1. Compute Dirichlet Values with the help of math-parser and FeSpace
    //2. Hand the value to the PDE for transformation according to time-Stepping ETC
    //3. Incorporate to OLAS (Penalty or projection)
    //
    // The general Scheme (also for Higher Order Entities) could be:
    // 1. Obtain eqn numbers from feSpace
    // 2. Get BC Coordinates from elemShape Object
    //    the elemSchape obj. asks the FeSpace for the local DOF Coordinates 
    //    and maps them to global
    // 3. Evaluate BC at the given coordinate and timestep
    // 4. Pass the value to OLAS
    //  
    Vector<Double> globCoord;
    
    // get global coordinate system and math parser
    CoordSystem * coosy = domain->GetCoordSystem();
    MathParser * parser = domain->GetMathParser();
    UInt dof;
    Double val;

    //===================================================
    // Homogenious BCs
    // Not treated right now due to old school equation numbering
    // ==================================================

    //loop over all homogenious BCs
    /*for ( UInt i = 0; i < hdBcs_.GetSize(); i++ ) {
      HomDirichletBc const & actBc = *(hdBcs_[i]);

      // Get EntityIterator
      EntityIterator it = actBc.entities->GetIterator();
      dof = actBc.dof;

      //RESULTINFO???
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        StdVector<Integer> eqns;
        feSpace_->GetEqns( eqns, it, dof  );
        //for( UInt iEqn = 0; iEqn < eqns.GetSize(); iEqn++){
        //   algsys_->SetDirichlet(  fctId_, eqns[iEqn], 0);

        //}
      }
    }*/

    //===================================================
    // InHomogenious BCs
    // ==================================================

    //loop over all inhomogeneous BCs
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
      InhomDirichletBc const & actBc = *(idBcs_[i]);

      // Get EntityIterator
      EntityIterator it = actBc.entities->GetIterator();

      //Get the grid Pointer 
      Grid * ptGrid = actBc.entities->GetGrid();
      dof = actBc.dof;

      //RESULTINFO???
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        StdVector<Integer> eqns, vEqns;
        // check type of space: if it is the "normal" H1, we set all of the
        // boundary nodes to the specified value.
        // If it is the H1Hi, we set the EDGE /FACE / INNER nodes to zero,
        // as we assume a linear boundary condition.
        if( feSpace_->IsHierarchical() ) {
         feSpace_->GetEqns(vEqns, it, dof, BaseFE::VERTEX ); 
        }
        feSpace_->GetEqns( eqns, it, dof  );

        for( UInt iEqn = 0; iEqn < eqns.GetSize(); iEqn++){
          UInt eqnNr = eqns[iEqn];
          //Evaluate the BC Value
          if( it.GetType() == EntityList::NODE_LIST ) {
            // Get node coordinate
            ptGrid->GetNodeCoordinate( globCoord, it.GetNode() );
            parser->SetCoordinates( mHandle_, *coosy, globCoord );
            //ResultCache::SetIndex(it.GetPos());
          }
          else {
            // register element midpoint
            const Elem* ptElem = it.GetElem();
            shared_ptr<ElemShapeMap> esm = 
                ptGrid->GetElemShapeMap(ptElem, false );
            esm->GetGlobMidPoint(globCoord);
            parser->SetCoordinates( mHandle_, *coosy, globCoord );
          }

          // Now evaluate value of IDBC
          //ResultCache::SetOutputType(ResultCache::OUT_AMPL);
          parser->SetExpr( mHandle_, actBc.value );
          val = parser->Eval( mHandle_ );

          if( this->GetTimeScheme() ) {
            this->GetTimeScheme()->AdaptBC(val,val,0,eqnNr);
          }
            //val = 1;

          //pde_->TransformBC(val, val, eqnNr);
          
          
          // if we have the higher order H1 space, we just set
          // the dirichelt values for the vertex unknowns.
          // The edge/face/inner ones are assumed to be constant.
          if( feSpace_->IsHierarchical()  
              && iEqn >= vEqns.GetSize() )  {
            algsys_->SetDirichlet(  fctId_, eqnNr, 0.0);
          } else {
            if ( std::tr1::is_same<T,Complex>::value ) {
              // === Harmonic ===
              parser->SetExpr( mHandle_, actBc.phase );
              Double phase = parser->Eval( mHandle_ );
              Complex complexValue( val * std::cos( phase / 180 * PI ),
                                    val * std::sin( phase / 180 * PI ) );
              algsys_->SetDirichlet(  fctId_, eqnNr, complexValue);
            } else {
              // === Transient ===
              algsys_->SetDirichlet(  fctId_, eqnNr, val);
            } // complex / transient
          } 
        } // for eqns
      }
    }
  }
  
  template<typename T>
  void FeFunction<T>::GetVector(Vector<T>& vec, 
                                const LocPointMapped& lpm ){
    assert(dimType_ == CoefFunction::VECTOR);
    
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
    assert(dimType_ == CoefFunction::SCALAR);
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
