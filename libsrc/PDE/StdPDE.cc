#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"


namespace CoupledField {


  StdPDE::StdPDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
                 WriteResults *aOutFile, TimeFunc *aTimeFunc)
    : BasePDE(),
      nonLin_(FALSE),
      incStopCrit_(1e-2), 
      residualStopCrit_(1e-3),
      firstTimeStepStatic_(TRUE),
      isaxi_(FALSE),
      isComplex_(FALSE),
      numPDENodes_(0),
      numElems_(0),
      sol_(NULL),
      isAlwaysStatic_(FALSE),
      dampingType_(NONE),
      needsDampingMatrix_(FALSE),
      isIncrFormulation_(FALSE),
      TS_alg_(NULL),
      materialData_(NULL),
      effectiveMass_(FALSE)
    
  {
    ENTER_FCN( "StdPDE::StdPDE");

   
    // =====================================================================
    // set file pointers
    // =====================================================================
    inFile_     = aInFile;
    outFile_    = aOutFile;
    ptTimeFunc_ = aTimeFunc;
    ptgrid_     = aptgrid;
    ptBCs_      = aptBCs;
    assemble_   = NULL;
    algsys_     = NULL;
    eqnData_    = NULL;
    
    // =====================================================================
    // set analysis parameters
    // =====================================================================
    actlevel_ = 0;
    actFrequency_ = 0;
    couplingBCsCounter_ = 0;
    numDirichletBCs_ = 0;
    pdeIsCoupled_ = FALSE;
    updateCouplingBCs_ = FALSE;
    dim_ = ptgrid_->GetDim();
    geoUpdate_ = FALSE;
    iterCoupledCounter_ = 0;
    effectiveMass_ = FALSE;

  }


  const Vector<Double>& StdPDE::getS1() const {
  
    ENTER_FCN( "StdPDE::getS1" );
  
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv1();
    }
    else {
      (*error) << pdename_ << ":getS1: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );
    }
  }
  
  
  const Vector<Double>& StdPDE::getS2() const {
    
    ENTER_FCN( "StdPDE::getS2" );
    
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv2();
    }
    else {
      (*error) << pdename_ << ":getS2: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );
    }
  }
  

  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  void StdPDE::GetElemCoords( const StdVector<Integer> connect, 
                              Matrix<Double> &coordMat, 
                              const Integer level ) {

    ENTER_FCN( "StdPDE::GetElemCoords" );

    Integer pdeNode;
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    if (deltCoords_.GetSizeRow() != 0 && geoUpdate_ == TRUE) {
      for (Integer i=0; i<coordMat.GetSizeRow(); i++) {
        for (Integer j=0; j<coordMat.GetSizeCol(); j++) {
          pdeNode = eqnData_->Mesh2PDENode(connect[j]);
          //coordMat(i,j) += deltCoords_(i,mesh2PDENode_ [connect[j] - 1] - 1);
          coordMat(i,j) += deltCoords_(i, pdeNode - 1);
        }
      }
    }
  }

  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double StdPDE::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "StdSolveStep::RhsL2Norm" );

    Integer node, dof, eqnNr, eqnDof;
  
    std::list<Integer> nodes;
  
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    // Eliminate dirichlet node from RHS (due to penalty formulation)
    for (Integer i=0; i< bcs_hd_.GetSize(); i++)
      {
        dof = 1;
        if ( dofspernode_ > 1 ) {
          dof = GetBCDof( homDirichDof_[i] );
        }

        nodes=ptBCs_->GetNodesLevel(bcs_hd_[i]);
      
        for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
          {
            node=*p;
            eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
            if (eqnNr != 0){
              actRHS[(eqnNr-1)*dofsPerEQN + eqnDof-1] = 0;
            }
          }
      }
    return actRHS.NormL2();
  }


  

  Integer StdPDE::GetBCDof( const std::string dofStartString ) {
    
    ENTER_FCN( "StdPDE::GetBCDof" );
    
    Integer nrActDof = 0;
    
    if ( dofStartString == "ux" )
      nrActDof = 1;
    if ( dofStartString == "uy" )
      nrActDof = 2;
    if ( dofStartString == "uz" )
      nrActDof = 3;
    // HARD coded for piezo case
    if ( dofStartString == "ep" )
      nrActDof = dofspernode_;
    if ( nrActDof == 0 )
      Error( "Unknown dof-type in homog. BC; substring must start with ux, uy,"
             "uz or ep!!", __FILE__, __LINE__);

    return nrActDof;
  }

  void StdPDE::GetMemento(PDEMemento & memento) {

    ENTER_FCN( "StdPDE::GetMemento" );

    // first get memento of coupling object
    if (pdeIsCoupled_) {
      ptCoupling_->GetMemento(memento.couplingMemento_);
    }

    // then write own data to PDEMemento
    memento.analysisType_ = analysistype_;

    if (memento.sol_)
      delete memento.sol_;

    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // --- Real values --
      // Set solution
      memento.sol_ = new Vector<Double>;
      dynamic_cast<Vector<Double>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).GetAlgSysVector();

      if (analysistype_ == TRANSIENT) {
        Warning( "Currently only the first derivative is stored in the "
                 "memento object!", __FILE__, __LINE__ );

        // Set first derivative
        memento.solDeriv1_ = getS1();   
        
        // Set second derivative
        memento.solDeriv2_ = getS2();
      }
    }
    else {

      // --- Complex values --      
      // Set solution
      memento.sol_ = new Vector<Complex>;
      dynamic_cast<Vector<Complex>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).GetAlgSysVector();
    } 

    // now memento is initialized
    memento.isSet_ = TRUE;
  }


  void StdPDE::SetMemento( PDEMemento &memento ) {

    ENTER_FCN( "StdPDE::SetMemento" );
  
    // if there is no information in the menmento just leave
    if ( memento.isSet_ == FALSE ) {
      return;
    }
  
    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // --- Real values --
      // Set solution
      dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).SetAlgSysVector
        (dynamic_cast<Vector<Double>&>(*(memento.sol_)));
      
      // if previous step was transient and the current step is also
      // then give the time derivative to the timestepping algorithm
      if (analysistype_ == TRANSIENT
          && memento.analysisType_ == TRANSIENT) {
        TS_alg_->SetDeriv1(memento.solDeriv1_);
        TS_alg_->SetDeriv2(memento.solDeriv2_);
      }
    }

    else {

      // --- Complex values --      
      // Set solution
      dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).SetAlgSysVector
        (dynamic_cast<Vector<Complex>&>(*(memento.sol_)));
    }

    // If PDE is coupled, set the according coupling memento
    if (pdeIsCoupled_) {
      ptCoupling_->SetMemento(memento.couplingMemento_);
    }
  }

  void StdPDE::CreateMatrices_Solver() {
    
    ENTER_FCN( "StdPDE::CreateMatrices_Solver" );
    
    // create and initialize matrices 
    assemble_->CreateMatrices();
    assemble_->InitMatrices();
    
    // create solver and preconditioner
    algsys_->CreateSolver();
    
    algsys_->CreatePrecond();

    // now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();
  }

  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================
  void StdPDE::SetAlgSys() {

    ENTER_FCN( "StdPDE::SetAlgSys" );

    // Set parameter for solver and preconditioner

    (*cla) <<  "--- PDE: " << pdename_ << " ---" << std::endl;

    // Set parameters for OLAS
    std::string amExpert;
    params->Get( "override", amExpert, "expert" );
    CFSOLASParams::SetParams( pdename_, params, olasParams_,(amExpert=="yes"));

    // Set the graph type used for the system matrices
    assemble_->SetupMatrixGraph(eqnData_->GetNumEQNs());

    // Allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }

  Double StdPDE::GetFracDampMatrixCoeff(Integer actSD) {
    
    ENTER_FCN( "StdPDE::GetFracDampMatrixCoeff" );

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    Double y  = materialData_[actSD].GetDampingBeta();
    Double dt = TS_alg_->GetTimeStep();

    coeff = exp(-(y-1.0)*log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

    return coeff;
  }
  
  void StdPDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                   StdVector<Integer>& connecth ) {

    ENTER_FCN( "StdPDE::GetSolVecOfElement" );

    // displacement of element nodes
    elemSol.Resize(dofspernode_ * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();
  
    for(Integer actNode=0; actNode<connecth.GetSize(); actNode++) {
      for(Integer actDof=0; actDof < dofspernode_; actDof++) {
        eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
        if (eqnNr!= 0) {
          elemSol[actDof + actNode*dofspernode_] =
            sol[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
        }
        else {
          elemSol[actDof + actNode*dofspernode_] = 0.0;
        }
      }
    }
  }

  
  void StdPDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                       StdVector<Integer>& connecth) {

    ENTER_FCN( "StdPDE::GetDerivSolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();
  
    if (analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der1 = getS1();
    
      for( Integer actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( Integer actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if ( eqnNr != 0 ) {
            sol[actDof + actNode*dofspernode_] =
              sol_der1[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }

  void StdPDE::GetDeriv2SolVecOfElement( Vector<Double>& sol,
                                         StdVector<Integer>& connecth ) {

    ENTER_FCN( "StdPDE::GetDeriv2SolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    if (analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der2 = getS2();
    
      for( Integer actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( Integer actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if (eqnNr!= 0) {
            sol[actDof + actNode*dofspernode_] =
              sol_der2[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }


  void StdPDE::GetDerivSolOfElement( Matrix<Double>& sol,
                                     StdVector<Integer>& connect_PDE ) {

    ENTER_FCN( "StdPDE::GetDerivSolOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_, connect_PDE.GetSize());

    const Vector<Double>& sol_der1 = getS1();  

    for( Integer actNode = 0; actNode < connect_PDE.GetSize(); actNode++ ) {
      for( Integer actDof = 0; actDof < dofspernode_; actDof++) {
        sol[actDof][actNode] =
          sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];
      }
    }
  }

  void StdPDE::sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted){

    ENTER_FCN( "StdPDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voit: Txx Tyy Tzz Tyz Txz Txy

    if (subType_ == "3d") {
      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }
    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
    else {
      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  void StdPDE::sortStresses( Vector<Complex> &unsorted,
                             Vector<Complex> &sorted ) {
    ENTER_FCN( "StdPDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voit: Txx Tyy Tzz Tyz Txz Txy

    if (subType_ == "3d") {

      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }

    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }

    else {

      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void StdPDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "StdPDE::StoreAlgsysToVec" );

    //const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();

    vec.Resize(numElems);
    for (Integer i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }


} // end of namespace
