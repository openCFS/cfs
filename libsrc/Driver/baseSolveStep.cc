#include <fstream>
#include <iostream>
#include <string>

#include "baseSolveStep.hh"


namespace CoupledField {

  BaseSolveStep::BaseSolveStep(BasePDE& apde) : PDE_(apde)
  {

    ENTER_FCN( "BaseSolveStep::BaseSolveStep" );
 
    pdename_      = PDE_.pdename_;
    numPDENodes_  = PDE_.numPDENodes_;
    isaxi_        = PDE_.isaxi_;
    subdoms_      = PDE_.subdoms_;
    materialData_ = PDE_.getPDEMaterialData();
    ptgrid_       = PDE_.ptgrid_;
    algsys_       = PDE_.algsys_;
    sol_          = PDE_.sol_;
    eqnData_      = PDE_.eqnData_;
    assemble_     = PDE_.assemble_;
    actlevel_     = PDE_.actlevel_;

    matrix_factor_ = PDE_.matrix_factor_;
    lasttimecalc_ = PDE_.lasttimecalc_;
    laststepcalc_ = PDE_.laststepcalc_;
    TS_alg_       = PDE_.TS_alg_;

    lineSearch_       = PDE_.lineSearch_;
    geoUpdate_        = PDE_.geoUpdate_;
    incStopCrit_      = PDE_.incStopCrit_;
    residualStopCrit_ = PDE_.residualStopCrit_;
    nonLinMaxIter_    = PDE_.nonLinMaxIter_;
    nonLinMethod_     = PDE_.nonLinMethod_;
    nonLinLogging_    = PDE_.nonLinLogging_ ;
	nonLinPDEName_    = PDE_.nonLinPDEName_;


    pdeIsCoupled_       = PDE_.pdeIsCoupled_;
    iterCoupledCounter_ = &PDE_.iterCoupledCounter_;

}

  
  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void BaseSolveStep::SolveStepStatic(const Integer kstep, const Double asteptime,
				      const Integer level, const Boolean reset) {

    ENTER_FCN( "BaseSolveStep::SolveStepStatic" );

    PDE_.lasttimecalc_ = asteptime;
    PDE_.laststepcalc_ = kstep;

    if (PDE_.nonLin_) {
      StepStaticNonLin(kstep,asteptime,level,reset);
    }
    else {
      StepStaticLin(kstep,asteptime,level,reset);
    }
  }


  void BaseSolveStep::StepStaticLin( const Integer kstep, const Double aTime,
                               const Integer level, const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::StepStaticLin" );

    Integer job = 3; // only update BCs
    Double * ptsol;
    PDE_.lasttimecalc_ = aTime; // for correct output in unv-file


    // If the geometry has changed or the system matrix
    // is calculated for the first time,
    // the matrices have to be reassembled and therfore
    // the preconditioner has to be recalculated

    if ( PDE_.geoUpdate_ == TRUE ||  PDE_.firstTimeStepStatic_ == TRUE) {
      PDE_.assemble_->AssembleMatrices(level);
      job = 1; // calc new preconditioner
    }
  
    // The RHS-sources have to be reassembled each time
    PDE_.assemble_->AssembleSrcRHS(level, aTime);

    PDE_.SetBCs(level, aTime);

    // Incorporate Boundary conitions and
    // recalc the prconditioner eventually
    PDE_.algsys_->BuildInDirichlet();
    PDE_.algsys_->SetupPrecond( job );
    PDE_.algsys_->SetupSolver( job );

    // Solve problem
    PDE_.algsys_->Solve();

    // Get the solution and store it
    ptsol = PDE_.algsys_->GetSolutionVal();
    PDE_.sol_->CopyFromAlgSysDataPointer(ptsol);

    PDE_.firstTimeStepStatic_ = FALSE;
  }


  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  // initialize PDEs before iteration (done for each time step)
  void BaseSolveStep::InitStepTransCoupled( Double asteptime ) {
    ENTER_FCN( "BaseSolveStep::InitStepTransCoupled" );
    PDE_.lasttimecalc_ = asteptime;
    PDE_.iterCoupledCounter_ = 0;
  }


  void BaseSolveStep::PreStepTrans( const Integer kstep, const Double asteptime,
                              const Integer level, const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::PreStepTrans" );

    PDE_.lasttimecalc_= asteptime;

    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    PDE_.algsys_->InitRHS();

    if (PDE_.geoUpdate_) {
      PDE_.algsys_->InitRHS();
      PDE_.algsys_->InitSol();
      PDE_.assemble_->InitMatrices();

      PDE_.assemble_->SetReassemble();   
    }
  }


  void BaseSolveStep::SolveStepTrans( const Integer kstep, const Double asteptime, 
                                const Integer level, const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::SolveStepTrans" );

    PDE_.lasttimecalc_= asteptime;
    PDE_.recalc_ = FALSE;

    if (PDE_.laststepcalc_ == kstep && kstep != 1) {
      PDE_.recalc_ = TRUE;
    }
    else {
      PDE_.laststepcalc_= kstep;
    }

    if (PDE_.nonLin_) {
      StepTransNonLin(kstep, asteptime, level, reset);
    }
    else {
      StepTransLin(kstep, asteptime, level, reset);
    }
  }


  //! \todo delete job parameter or replace it by a 
  //! meaningfull attribute
  void BaseSolveStep::StepTransLin( const Integer kstep, const Double asteptime,
                              const Integer level, const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::StepTransLin" );

    Double * ptsol;
    Integer job;

    //account for RHS
    PDE_.assemble_->AssembleSrcRHS(level,PDE_.lasttimecalc_);

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(PDE_.sol_);

    if ( PDE_.pdeIsCoupled_ == FALSE || PDE_.iterCoupledCounter_ == 0 ) {        
      Vector<Double> & solvector= solhelp->GetAlgSysVector();
      PDE_.TS_alg_->Predictor(solvector);
    }

    if ( PDE_.laststepcalc_ == 1 ) {
      job = 3;

      // why is the first statement checking for 'pdeIsCoupled'?
      if ( PDE_.pdeIsCoupled_ == FALSE || PDE_.iterCoupledCounter_ == 0 
           || PDE_.geoUpdate_ == TRUE ) {
        job = 1;
        PDE_.assemble_->AssembleMatrices(level);
        PDE_.algsys_->ConstructEffectiveMatrix(PDE_.matrix_factor_);
      }  
    }
    else if (reset) {
      job = 1;

      PDE_.algsys_->InitMatrix(SYSTEM);
      PDE_.algsys_->InitMatrix(STIFFNESS);
      PDE_.algsys_->InitMatrix(MASS);
      if (PDE_.dampingType_) {
        PDE_.algsys_->InitMatrix(DAMPING);
      }
      PDE_.algsys_->ConstructEffectiveMatrix(PDE_.matrix_factor_);
    }
    else {
      job = 3;

      // The following section is only an experiment up to now
      if ( PDE_.geoUpdate_ == TRUE && PDE_.pdeIsCoupled_ == TRUE ) {
        if (PDE_.isIncrFormulation_) {
          Error( "Incremental formulation and geoUpdate are currently not "
                 "working together" );
        }
        job = 1;
        PDE_.assemble_->AssembleMatrices(level);
        PDE_.algsys_->ConstructEffectiveMatrix(PDE_.matrix_factor_);      
      }
    }

    if (PDE_.isIncrFormulation_) {
      Vector<Double> & solvector =
        dynamic_cast<NodeStoreSol<Double>*>(PDE_.sol_)->GetAlgSysVector();

      solvector *= -1;
      PDE_.algsys_->UpdateRHS(SYSTEM,solvector.GetPointer());
    }

    PDE_.TS_alg_->UpdateRHS();

    PDE_.SetBCs( level, PDE_.lasttimecalc_ );
    PDE_.algsys_->BuildInDirichlet();

    if ( job == 1 ) {
      PDE_.algsys_->SetupPrecond( job );
      PDE_.algsys_->SetupSolver( job );
    }

    PDE_.algsys_->Solve();
    ptsol = PDE_.algsys_->GetSolutionVal();

    if ( PDE_.isIncrFormulation_ ) {

      //what a heuristic!!!!!
      Double relaxVal = 0.5;

      if (PDE_.iterCoupledCounter_ == 0)
        relaxVal = 0.1;
      if (PDE_.iterCoupledCounter_ == 1)
        relaxVal = 0.25;

      PDE_.StoreAlgsysToVec(PDE_.solIncr_, ptsol);
      if (PDE_.iterCoupledCounter_ == 0)
        PDE_.actSol_ = PDE_.solIncr_*relaxVal;
      else 
        PDE_.actSol_ += PDE_.solIncr_*relaxVal;

      PDE_.sol_->SetAlgSysVector(PDE_.actSol_);
    }
    else
      PDE_.sol_->CopyFromAlgSysDataPointer(ptsol);

    Vector<Double> & solvector =\
      dynamic_cast<NodeStoreSol<Double>*>(PDE_.sol_)->GetAlgSysVector();

    if (!PDE_.pdeIsCoupled_)
      PDE_.TS_alg_->Corrector(solvector);
  }


  void BaseSolveStep::PostStepTrans( const Integer kstep, const Double asteptime,
                               const Integer level ) {

    ENTER_FCN( "BaseSolveStep::PostStepTrans" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(PDE_.sol_);
    
    if ( PDE_.pdeIsCoupled_ ) {

      //save solution
      Vector<Double> & solvector= solhelp->GetAlgSysVector();

      //perform corrector step
      PDE_.TS_alg_->Corrector(solvector); 
    }
  
    if (PDE_.pdeIsCoupled_) {
      PDE_.iterCoupledCounter_++;
    }
  }



  // ======================================================
  // Solve Step Harmonic  SECTION  
  // ======================================================

  void BaseSolveStep::PreStepHarmonic( const Integer freqStep,
                                 const Double frequency, const Integer level,
                                 const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::PreStepHarmonic" );

    PDE_.actFrequency_ = frequency;
    PDE_.actFreqStep_  = freqStep;
    PDE_.assemble_->SetFrequency(frequency);
    PDE_.algsys_->InitRHS();

    if (reset) {
      PDE_.assemble_->InitMatrices();
    }
  }


  void BaseSolveStep::SolveStepHarmonic( const Integer freqStep,
                                   const Double frequency, const Integer level,
                                   const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::SolveStepHarmonic" );

    if ( PDE_.nonLin_ ) {
      StepHarmonicNonLin( freqStep, frequency, level, reset );
    }
    else {
      StepHarmonicLin( freqStep, frequency, level, reset );
    }
  }


  void BaseSolveStep::StepHarmonicLin( const Integer freqStep,
                                 const Double frequency, const Integer level,
                                 const Boolean reset ) {

    ENTER_FCN( "BaseSolveStep::StepHarmonicLin" );

    Integer job;
    Double * ptsol;

    if ( reset ) {
      PDE_.assemble_->AssembleMatrices( level );
    }


    //this has to be done each time!
    PDE_.assemble_->AssembleSrcRHS(level, frequency);

    if ( reset ) {

      //account for bcs
      PDE_.SetBCs(level, frequency);

      // calc new preconditioner
      job = 1;
    }
    else {
      job = 3;
    }

    PDE_.algsys_->BuildInDirichlet();
    
    if (job == 1) {
      PDE_.algsys_->SetupPrecond( job );
      PDE_.algsys_->SetupSolver( job );
    }

    PDE_.algsys_->Solve();

    ptsol = PDE_.algsys_->GetSolutionVal();

    PDE_.sol_->CopyFromAlgSysDataPointer(ptsol);
  }


  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double BaseSolveStep::SetLinRHS(const Integer level)
  {
    ENTER_FCN( "BaseSolveStep::SetLinRHS" );

    Double RhsLinL2Norm;  


    // to incorporate loads
    assemble_->AssembleSrcRHS(level, lasttimecalc_); 

    // Stores rhs vector into extForces and returns that L2-norm
    StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

    RhsLinL2Norm = RhsLinVal_.NormL2();
 
    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      Warning("Zero external force vector!! ", __FILE__,__LINE__);
    }

    return RhsLinL2Norm;
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void BaseSolveStep::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "BaseSolveStep::StoreAlgsysToVec" );

    //const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    if ( vec.GetSize() == 0 ) {
      vec.Resize(numElems);
    }

    for (Integer i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }


  // calculates L2-norm of RHS regarding dirichlet entries due to penalty 
  // formulation by setting them 0
  Double BaseSolveStep::RhsL2Norm(Vector<Double>& actRHS)
  {
    ENTER_FCN( "BaseSolveStep::RhsL2Norm" );

    Integer node, dof, eqnNr, eqnDof;
  
    std::list<Integer> nodes;
  
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    // Eliminate dirichlet node from RHS (due to penalty formulation)
    for (Integer i=0; i< PDE_.bcs_hd_.GetSize(); i++)
      {
	dof = 1;
	if ( PDE_.dofspernode_ > 1 ) {
	  dof = PDE_.GetBCDof( PDE_.homDirichDof_[i] );
	}

	nodes=PDE_.ptBCs_->GetNodesLevel(PDE_.bcs_hd_[i]);
      
	for (std::list<Integer>::const_iterator p=nodes.begin(); p!=nodes.end(); p++)
	  {
            node=*p;
            PDE_.eqnData_->Node2EQN(node,dof,eqnNr,eqnDof);
            if (eqnNr != 0){
              actRHS[(eqnNr-1)*dofsPerEQN + eqnDof-1] = 0;
            }
	  }
      }
    return actRHS.NormL2();
  }



  Double BaseSolveStep::LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
			     Double& etaLineSearch, Integer level, Boolean trans)
  {
    ENTER_FCN( "BaseSolveStep::LineSearch" );

    Vector<Double> solOld(actSol);
    const Integer nrEtas = 4;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125};
    Double etaOpt;
    Double residualL2NormOpt = 1e15;
  
    for(Integer i=0; i<nrEtas; i++)
      {
	actSol = solIncrement * eta[i];
	actSol += solOld;

	sol_->SetAlgSysVector(actSol);
	//StoreVecToSolArray(actSol);

	// recalculate RHS with new values to get new residual (f^(k+1))========
	algsys_->InitRHS(RhsLinVal_.GetPointer());

	if(trans)
	  {
	    assemble_->AssembleNLRHS(level, lasttimecalc_);
	    TS_alg_->UpdateRHS(actSol);
	  }
	else
	  assemble_->AssembleNLRHS(level);



	// =====================================================================
	// calculation of error norms
	// =====================================================================
	Vector<Double> actRHS;
	StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );

	// calculation of residual error =======================================
	Double residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )

	if (residualL2Norm < residualL2NormOpt)
	  {
	    residualL2NormOpt = residualL2Norm;
	    etaOpt = eta[i];
	  }
      }

    etaLineSearch = etaOpt;
  
    actSol =  solIncrement * etaOpt;
    actSol += solOld;

    return residualL2NormOpt;
  }


  // returns that L2-norm of an algsys vector
  Double BaseSolveStep::AlgsysL2Norm(Double * pt)
  {
    ENTER_FCN( "BaseSolveStep::AlgsysL2Norm" );

    //const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();
    Double quadSum = 0;
  
    for (Integer i=0; i<numElems; i++)   
      quadSum += pt[i]*pt[i];

    return sqrt(quadSum);
  }
  

  void BaseSolveStep::WriteClaNlNorms(const Integer iterationCounter, 
				      const Double residualL2Norm, 
				      const Double extForcesL2Norm,
				      const Double residualErr, 
				      const Double solIncrL2Norm, 
				      const Double actSolL2Norm, 
				      const Double incrementalErr)
  {
    ENTER_FCN( "BaseSolveStep::WriteClaNlNorms" );
  
    *cla << std::endl << " ======================================================= " 
	 << std::endl
	 << " NONLINEAR ITERATION " << iterationCounter << std::endl
	 << " ======================================================= " << std::endl;
    *cla << " === Residual norm:          " << residualL2Norm << std::endl;
    *cla << "     Norm of ext. forces:    " << extForcesL2Norm << std::endl;
    *cla << "     Residual error          " << residualErr << std::endl;
  
    *cla << " === Incremental sol L2Norm: " << solIncrL2Norm << std::endl;
    *cla << "     Actual solution L2Norm: " << actSolL2Norm << std::endl;
    *cla << "     Incremental error       " << incrementalErr << std::endl;      
  }



  // ======================================================
  // get BasePDE-methods  
  // ======================================================


  void BaseSolveStep::SetBCs( const Integer level,  const Double time )
  {
   PDE_.SetBCs(level, time);
  }

  void BaseSolveStep::GetElemCoords(const StdVector< Integer > connect,
					   Matrix< Double > &coordMat,
					   const Integer level)
  {
    PDE_.GetElemCoords(connect, coordMat, level);
  }

  void BaseSolveStep::GetSolVecOfElement(Vector<Double>& sol, 
					 StdVector<Integer>& connect_PDE)
  {
    PDE_.GetSolVecOfElement(sol, connect_PDE);
  }

  void BaseSolveStep::GetDerivSolVecOfElement(Vector<Double>& sol, 
						     StdVector<Integer>& connect_PDE)
  {
    PDE_.GetDerivSolVecOfElement(sol, connect_PDE);
  }

  void BaseSolveStep::GetDeriv2SolVecOfElement(Vector<Double>& sol, 
						      StdVector<Integer>& connect_PDE)
  {
    PDE_.GetDeriv2SolVecOfElement(sol, connect_PDE);
  }


  //   Default Destructor
  // **********************
  BaseSolveStep::~BaseSolveStep() {

    ENTER_FCN( "BaseSolveStep::~BaseSolveStep" );
 
  }


} // end of namespace

