// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "fefunction.hh"
#include "PDE/SinglePDE.hh"
#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  BaseFeFunction::BaseFeFunction(){
    mHandle_ = domain->GetMathParser()->GetNewHandle();
    isoOrder_ = 0;
    anIsoOrder_.Resize(0);

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
  
  void BaseFeFunction::AddEntityList( shared_ptr<EntityList> list ){
    entities_.Push_back(list);
  }
  
  StdVector< shared_ptr<EntityList> > BaseFeFunction::GetEntityList(){
    return entities_;
  }
  
  //! Set the PDE pointer of the function
  void BaseFeFunction::SetPDE(shared_ptr<SinglePDE> pde){
    pde_ = pde;
  }

  //! Get the PDE Pointer
  shared_ptr<SinglePDE> BaseFeFunction::GetPDE(){
    return pde_;
  }

  void BaseFeFunction::SetGrid(shared_ptr<Grid> grid){
    grid_ = grid;
  }

  shared_ptr<Grid> BaseFeFunction::GetGrid(){
    return grid_;
  }

  void BaseFeFunction::SetSystem( shared_ptr<BaseSystem> sys){
    algsys_ = sys;
  }

  shared_ptr<BaseSystem> BaseFeFunction::GetSystem(){
    return algsys_;
  }

  void BaseFeFunction::SetIsoOrder(UInt order){
    isoOrder_ = order;
  }

  UInt BaseFeFunction::GetIsoOrder(){
    return isoOrder_;
  }

  void BaseFeFunction::SetAnIsoOrder(StdVector<UInt> order){
    anIsoOrder_ = order;
  }

  StdVector<UInt> BaseFeFunction::GetAnIsoOrder(){
    return anIsoOrder_;
  }

  void BaseFeFunction::AddHomDirichletBc( shared_ptr<HomDirichletBc> bc ){
    hdBcs_.Push_back(bc);
  }

  void BaseFeFunction::AddInhomDirichletBc( shared_ptr<InhomDirichletBc> bc ){
    idBcs_.Push_back(bc);
  }

  void BaseFeFunction::AddInhomNeumannBC( shared_ptr<InhomNeumannBc> bc ){
    inBcs_.Push_back(bc);
  }

  void BaseFeFunction::AddConstraint( shared_ptr<Constraint> bc ){
    constraints_.Push_back(bc);
  }

  //! Get solution for specific entity
  void BaseFeFunction::GetEntitySolution( SingleVector& elemSol, 
                        const EntityIterator& it ){
  }
                        
  //! Get solution as matrix for specific entity
  void BaseFeFunction::GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                  const EntityIterator& it ){
  }


  template<typename T>
  FeFunction<T>::FeFunction(){
  }

  template<typename T>
  FeFunction<T>::~FeFunction(){
  }

  template<typename T>
  void FeFunction<T>::GetEntitySolution( SingleVector& elemSol, 
                        const EntityIterator& it ){
    Vector<T> & temp = dynamic_cast<Vector<T>&>(elemSol);
    StdVector<Integer> eqns;
    feSpace_->GetEqns(eqns, it);
    temp.Resize(eqns.GetSize());
    for(UInt iDof = 0 ; iDof < eqns.GetSize(); iDof++){
      temp[iDof] = coeffs_[eqns[iDof]]; 
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
    elemSol.Resize(feSpace_->GetNumFunctions(it), dofsPerUnknown);
    for(UInt iDof = 0 ; iDof < dofsPerUnknown ; iDof++){
      feSpace_->GetEqns(eqns, it,iDof);
      for(UInt iEqn = 0;iEqn < eqns.GetSize() ; iEqn++){
        temp[iDof][iEqn] = coeffs_[eqns[iEqn]];
      }
    }
  }


  void FeFunction<Complex>::ApplyBC(){
    EXCEPTION("Complex Valued BCs are not supported yet");
  }

  void FeFunction<Double>::ApplyBC(){
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

    //loop over all homogenious BCs
    for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
      InhomDirichletBc const & actBc = *(idBcs_[i]);

      // Get EntityIterator
      EntityIterator it = actBc.entities->GetIterator();

      //Get the grid Pointer 
      Grid * ptGrid = actBc.entities->GetGrid();
      dof = actBc.dof;

      //RESULTINFO???
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        StdVector<Integer> eqns;
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
          pde_->TransformBC(val, val, eqnNr);

          algsys_->SetDirichlet(  fctId_, eqnNr, val);
        }
      }
    }
  }


}
