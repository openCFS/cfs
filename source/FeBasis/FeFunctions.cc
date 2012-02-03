#include "FeFunctions.hh"

#include "PDE/SinglePDE.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <boost/tr1/type_traits.hpp>

namespace CoupledField {
DECLARE_LOG(fefunc)
 DEFINE_LOG(fefunc, "feFunction")
 
  BaseFeFunction::BaseFeFunction(){
    
    fctId_ = NO_FCT_ID;
    mHandle_ = domain->GetMathParser()->GetNewHandle();
    algsys_ = NULL;

  }
  
  BaseFeFunction::~BaseFeFunction(){
  }
  
  void BaseFeFunction::SetResultInfo( shared_ptr<ResultInfo> info ){
    result_ = info;
  }
  
  shared_ptr<ResultInfo> BaseFeFunction::GetResultInfo(){
    return result_;
  }
  
  void BaseFeFunction::SetFeSpace( shared_ptr<FeSpace> space ){
    feSpace_ = space;
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
    
    // Note: As the shared_ptr to an Entitylist is not a unique
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
    }
  }
  
  StdVector< shared_ptr<EntityList> > BaseFeFunction::GetEntityList(){
    return entities_;
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

  void BaseFeFunction::AddInhomNeumannBC( shared_ptr<InhomNeumannBc> bc ){
    inBcs_.Push_back(bc);
    entities_.Push_back(bc->entities);
  }

  void BaseFeFunction::AddConstraint( shared_ptr<Constraint> bc ){
    constraints_.Push_back(bc);
  }
  
  // ========================================================================
  //  T E M P L A T I Z E D    V E R S I O N 
  // ========================================================================

  template<typename T>
  FeFunction<T>::FeFunction(){
  }

  template<typename T>
  FeFunction<T>::~FeFunction(){
  }
  
  template<typename T>
  void FeFunction<T>::Finalize(){
    
    // assert that functionId was set
    WARN("Add some more consistency checks here");
    if (fctId_ == NO_FCT_ID ) {
      EXCEPTION("No fctId was set!");
    }
        
    coeffs_.Resize( feSpace_->GetNumEquations());
    coeffs_.Init();
  }

  template<typename T>
  void FeFunction<T>::ExtractResult( shared_ptr<BaseResult> baseRes ) {

    ResultInfo& resInfo = *(baseRes->GetResultInfo() );
    //UInt numDofs = resInfo.dofNames.GetSize();
    shared_ptr<EntityList> list = baseRes->GetEntityList();
    Vector<T> & actSol = dynamic_cast<Result<T>&>(*baseRes).GetVector();
    actSol.Resize( list->GetSize() * resInfo.dofNames.GetSize() );
    EntityIterator it = list->GetIterator();
    actSol.Init();

    StdVector<Integer> eqnNums;
    UInt pos = 0;
    for ( it.Begin(); !it.IsEnd(); it++ ) {

      // get equation numbers
      feSpace_->GetEqns( eqnNums, it );
      for ( UInt iDof = 0; iDof < eqnNums.GetSize(); iDof++ ){
   
        // check for homogeneous Dirichlet boundary condition
        if ( eqnNums[iDof] != 0 ) {
          actSol[pos++] = coeffs_[abs(eqnNums[iDof])-1];
        } else {
          actSol[pos++] = 0.0;
        }
      }
    }
  }
  
  template<typename T>
  void FeFunction<T>::GetEntitySolution( SingleVector& elemSol, 
                                         const EntityIterator& it ){
    Vector<T> & temp = dynamic_cast<Vector<T>&>(elemSol);
    StdVector<Integer> eqns;
    feSpace_->GetEqns(eqns, it);
    temp.Resize(eqns.GetSize());
    for(UInt iDof = 0 ; iDof < eqns.GetSize(); iDof++){
      if( eqns[iDof] != 0 ) {
        temp[iDof] = coeffs_[eqns[iDof]-1];
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
     temp.Resize(feSpace_->GetNumFunctions(it), dofsPerUnknown);
     for(UInt iDof = 0 ; iDof < dofsPerUnknown ; iDof++){
       feSpace_->GetEqns(eqns, it,iDof);
       for(UInt iEqn = 0;iEqn < eqns.GetSize() ; iEqn++){
         if( eqns[iEqn] > 0 ) {
         temp[iDof][iEqn] = coeffs_[eqns[iEqn]-1];
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
    
    feSpace_->GetElemEqns(eqns, elem);
    elemSol.Resize(eqns.GetSize());
    for(UInt i= 0 ; i< eqns.GetSize(); i++){
      if( eqns[i] != 0 ) {
        elemSol[i] = coeffs_[std::abs(eqns[i])-1]; 
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
            // this case needs to be implemented ...
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
  // Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class FeFunction<Double>;
  template class FeFunction<Complex>;
#endif
}
