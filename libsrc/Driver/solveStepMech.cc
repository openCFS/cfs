#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMech.hh"


namespace CoupledField {

  SolveStepMech::SolveStepMech(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepMech::SolveStepMech" );
  }

  
  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================


  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================

  void SolveStepMech:: PreStepStatic(const Integer kstep, const Double asteptime,
			       const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepMech::PreStepStatic" );

    if (pdeIsCoupled_)
      // init RHS at this place, because forces of other PDEs are added to RHS afterwards
      algsys_->InitRHS();     
    else {
      algsys_->InitRHS();        
      // if PDE is coupled, the solution of the prior outer loops must be kept
      algsys_->InitSol();
    }
 
  }


  void SolveStepMech::StepStaticLin( const Integer kstep, const Double aTime,
                               const Integer level, const Boolean reset ) {

    ENTER_FCN( "SolveStepMech::StepStaticLin" );

    Integer job = 3; // only update BCs
    Double * ptsol;
    lasttimecalc_ = aTime; // for correct output in unv-file


    // If the geometry has changed or the system matrix
    // is calculated for the first time,
    // the matrices have to be reassembled and therfore
    // the preconditioner has to be recalculated

    if ( geoUpdate_ == TRUE ||  firstTimeStepStatic_ == TRUE) {
      assemble_->AssembleMatrices(level);
      assemble_->AssembleSprings(level, lasttimecalc_);
      job = 1; // calc new preconditioner
    }
  
    // The RHS-sources have to be reassembled each time
    assemble_->AssembleSrcRHS(level, aTime);

    SetBCs(level, aTime);

    // Incorporate Boundary conitions and
    // recalc the prconditioner eventually
    algsys_->BuildInDirichlet();
    algsys_->SetupPrecond( job );
    algsys_->SetupSolver( job );

    // Solve problem
    algsys_->Solve();

    // Get the solution and store it
    ptsol = algsys_->GetSolutionVal();
    sol_->CopyFromAlgSysDataPointer(ptsol);

    firstTimeStepStatic_ = FALSE;
  }


  void SolveStepMech::StepStaticNonLin(const Integer kstep, const Double aTime,
				 const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepMech::SolveStepStaticNonLin" );

    const Integer job = 1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
    Vector<Double>  actSol = solHelp.GetAlgSysVector();

     Vector<Double> solIncrement;
    solIncrement.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());

    PDE_.SetBCs(level, 0);

    // store linear part of RHS
    Double extForcesL2Norm = SetLinRHS(level); 

    assemble_->AssembleNLRHS(level);

   std::cout << "In SolveStepStaticNonLin" << std::endl;
    do
      {
	iterationCounter++;

	std::cout << std::endl << "Nonlinear Mechanics: Perform internal loop nr. " 
		  << iterationCounter << std::endl;      

#ifdef DEBUG
	*debug << std::endl
	       << "====================================================== "
	       << std::endl
	       << "Nonlinear Mechanics: Perform internal loop nr. "
	       << iterationCounter << std::endl;
#endif


	// setup and solve new system (rhs is already set) =====================
	assemble_->InitNonLinMatrices();
	assemble_->AssembleMatrices(level);
      
	algsys_->BuildInDirichlet();

	if (job == 1) {
	  algsys_->SetupPrecond(job);
	  algsys_->SetupSolver(job);
	}

	algsys_->Solve();

      

	// new solution is only an increment of the full solution =============
	StoreAlgsysToVec(solIncrement, algsys_->GetSolutionVal() );

	Double residualL2Norm = 0;
	Double etaLineSearch=0;

	if ( lineSearch_ != "no" )
	  actSol += solIncrement;
	else
	  // TRUE is for transient simulation
	  residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch, level);
      
	sol_->SetAlgSysVector(actSol);

	// recalculate RHS with new values to get new residual (f^(k+1))========
	algsys_->InitRHS(RhsLinVal_.GetPointer());

	assemble_->AssembleNLRHS(level);  // inner forces due to nonlin formulation


	// =====================================================================
	// calculation of error norms
	// =====================================================================
	if ( lineSearch_ != "no" )
	  {
	    Vector<Double> actRHS;
	    StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
          
	    // calculation of residual error =======================================
	    residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
	  }


	// calculation of residual error =======================================
	Double residualErr = residualL2Norm / extForcesL2Norm;


	// calculate incremental error ========================================
	Double solIncrL2Norm = solIncrement.NormL2();
	Double actSolL2Norm  = actSol.NormL2();
      
	Double incrementalErr;
      
	if (actSolL2Norm)
	  incrementalErr = solIncrL2Norm / actSolL2Norm;
	else
	  {
	    incrementalErr = solIncrL2Norm;
	    Warning("Zero solution vector!! ", __FILE__,__LINE__);      
	  }
      


	// =====================================================================
	// output of norms and data
	// =====================================================================


	WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
			solIncrL2Norm, actSolL2Norm, incrementalErr);

      
	Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, 
			      etaLineSearch);


	// boolean variable, holds condition if another iteration step is necessary
	performOneMoreStep = 
	  (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);      
      
      
	if (!(performOneMoreStep && iterationCounter < nonLinMaxIter_))
	  mycout << "incrementalErr " << incrementalErr << myendl 
		 << "incStopCrit_ " << incStopCrit_ << myendl
		 << "residualErr " << residualErr  << myendl
		 << "residualStopCrit_ " << residualStopCrit_ << myendl;

      }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  }


  void SolveStepMech :: PostStepStatic(const Integer kstep, const Double asteptime,
				 const Integer level)
  {
    ENTER_FCN( "SolveStepMech::PostStepStatic" );

    if (pdeIsCoupled_)
      (*iterCoupledCounter_)++;

  }



  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepMech::PreStepTrans( const Integer kstep, const Double asteptime,
                              const Integer level, const Boolean reset ) {

    ENTER_FCN( "SolveStepMech::PreStepTrans" );

    lasttimecalc_= asteptime;

    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();

    if (geoUpdate_) {
      algsys_->InitRHS();
      algsys_->InitSol();
      assemble_->InitMatrices();

      assemble_->SetReassemble();  
    }
  }


  void SolveStepMech::SolveStepTrans( const Integer kstep, const Double asteptime, 
                                const Integer level, const Boolean reset ) {

    ENTER_FCN( "SolveStepMech::SolveStepTrans" );

    lasttimecalc_= asteptime;
    recalc_ = FALSE;

    if (laststepcalc_ == kstep && kstep != 1) {
      recalc_ = TRUE;
    }
    else {
      laststepcalc_= kstep;
    }

    if (nonLin_) {
      StepTransNonLin(kstep, asteptime, level, reset);
    }
    else {
      StepTransLin(kstep, asteptime, level, reset);
    }
  }


  //! \todo delete job parameter or replace it by a 
  //! meaningfull attribute
  void SolveStepMech::StepTransLin( const Integer kstep, const Double asteptime,
                              const Integer level, const Boolean reset ) {

    ENTER_FCN( "SolveStepMech::StepTransLin" );

    Double * ptsol;
    Integer job;
    lasttimecalc_= asteptime;

    //account for RHS
    assemble_->AssembleSrcRHS(level,lasttimecalc_);

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    if ( pdeIsCoupled_ == FALSE || *iterCoupledCounter_ == 0 ) {        
      Vector<Double> & solvector= solhelp->GetAlgSysVector();
      TS_alg_->Predictor(solvector);
    }

    if ( laststepcalc_ == 1 ) {
      job = 3;

      // why is the first statement checking for 'pdeIsCoupled'?
      if ( pdeIsCoupled_ == FALSE || *iterCoupledCounter_ == 0 
           || geoUpdate_ == TRUE ) {
        job = 1;
        assemble_->AssembleMatrices(level);
        assemble_->AssembleSprings(level, lasttimecalc_);
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
      }  
    }
    else if (reset) {
      job = 1;

      algsys_->InitMatrix(SYSTEM);
      algsys_->InitMatrix(STIFFNESS);
      algsys_->InitMatrix(MASS);
      if (dampingType_) {
	algsys_->InitMatrix(DAMPING);
      }
      assemble_->AssembleSprings(level, lasttimecalc_);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
    }
    else {
      job = 3;

      // The following section is only an experiment up to now
      if ( geoUpdate_ == TRUE && pdeIsCoupled_ == TRUE ) {
        if (isIncrFormulation_) {
          Error( "Incremental formulation and geoUpdate are currently not "
                 "working together" );
        }
        job = 1;
        assemble_->AssembleMatrices(level);
        assemble_->AssembleSprings(level, lasttimecalc_);
        algsys_->ConstructEffectiveMatrix(matrix_factor_);      
      }
    }

    if (isIncrFormulation_) {
      Vector<Double> & solvector =
        dynamic_cast<NodeStoreSol<Double>*>(sol_)->GetAlgSysVector();

      solvector *= -1;
      algsys_->UpdateRHS(SYSTEM,solvector.GetPointer());
    }

    TS_alg_->UpdateRHS();

    PDE_.SetBCs( level, lasttimecalc_);
    algsys_->BuildInDirichlet();

    if ( job == 1 ) {
      algsys_->SetupPrecond( job );
      algsys_->SetupSolver( job );
    }

    algsys_->Solve();
    ptsol = algsys_->GetSolutionVal();

    if ( isIncrFormulation_ ) {

      //what a heuristic!!!!!
      Double relaxVal = 0.5;

      if (*iterCoupledCounter_ == 0)
        relaxVal = 0.1;
      if (*iterCoupledCounter_ == 1)
        relaxVal = 0.25;

      StoreAlgsysToVec(solIncr_, ptsol);
      if (*iterCoupledCounter_ == 0)
        actSol_ = solIncr_*relaxVal;
      else 
        actSol_ += solIncr_*relaxVal;

      sol_->SetAlgSysVector(actSol_);
    }
    else
      sol_->CopyFromAlgSysDataPointer(ptsol);

    Vector<Double> & solvector =\
      dynamic_cast<NodeStoreSol<Double>*>(sol_)->GetAlgSysVector();

    if (!pdeIsCoupled_)
      TS_alg_->Corrector(solvector);
  }




  void SolveStepMech::PostStepTrans( const Integer kstep, const Double asteptime,
                               const Integer level ) {

    ENTER_FCN( "SolveStepMech::PostStepTrans" );

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
    if ( pdeIsCoupled_ ) {

      //save solution
      Vector<Double> & solvector= solhelp->GetAlgSysVector();

      //perform corrector step
      TS_alg_->Corrector(solvector); 
    }
  
    if (pdeIsCoupled_) {
      iterCoupledCounter_++;
    }
  }







  void SolveStepMech::StepTransNonLin(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepMech::StepTransNonLin" );

    const Integer job = 1;
  
    static Integer timeStepCounter=1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;

    Vector<Double> actSol;
    Vector<Double> solIncrement;
  
    //is already done in PreStepTrans (basePDE.cc!!)
    //  algsys_->InitRHS();

    // Cast BaseStoreSol into StoreSol<Double>,
    // since this function is only called
    // in the transient case
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  

    actSol = solhelp->GetAlgSysVector();
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    //! store linear part of RHS
    Double extForcesL2Norm = SetLinRHS(level);

    timeStepCounter++;

    //update RHS due to time stepping
    TS_alg_->UpdateRHS(actSol);

    //! account for Dirichlet BCs
    SetBCs(level, lasttimecalc_);

  
    do
      {
	iterationCounter++;
	std::cout << std::endl << "Nonlinear Mechanics: Perform internal loop nr. " 
		  << iterationCounter << std::endl;

#ifdef DEBUG
	*debug << std::endl << "====================================================== " << std::endl
	       << "Nonlinear Mechanics: Perform internal loop nr. " << iterationCounter << std::endl;      
#endif
 


	// setup and solve new system (rhs is already set) =====================
	assemble_->InitNonLinMatrices();
	assemble_->AssembleMatrices(level);
	algsys_->ConstructEffectiveMatrix(matrix_factor_);

	TS_alg_->UpdateRHS(actSol);
	SetBCs(level, lasttimecalc_);

	algsys_->BuildInDirichlet();

	if (job == 1) {
	  algsys_->SetupPrecond(job);
	  algsys_->SetupSolver(job);
	}

	algsys_->Solve();

	// new solution is only an increment of the full solution =============
	StoreAlgsysToVec(solIncrement, algsys_->GetSolutionVal() );

	Double residualL2Norm;
	Double etaLineSearch = 0;
      
	if ( lineSearch_ != "no" )
	  actSol += solIncrement;
	else
	  // TRUE is for transient simulation
	  residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch, level, TRUE);

	sol_->SetAlgSysVector(actSol);

	// recalculate RHS with new values to get new residual (f^(k+1))========
	algsys_->InitRHS();
	assemble_->AssembleSrcRHS(level,lasttimecalc_); 
	assemble_->AssembleNLRHS(level, lasttimecalc_);  // inner forces due to nonlin formulation
	TS_alg_->UpdateRHS(actSol);




	// =====================================================================
	// calculation of error norms
	// =====================================================================

	if ( lineSearch_ != "no" )
	  {
	    Vector<Double> actRHS;
	    StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
          
	    // calculation of residual error =======================================
	    residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
	  }
      
	Double residualErr    = residualL2Norm / extForcesL2Norm;
      
	if (extForcesL2Norm < NORM_EPS)
	  residualErr = residualL2Norm;  // take the absolute error instead of the relative


	// calculate incremental error ========================================
	Double solIncrL2Norm = solIncrement.NormL2();
	Double actSolL2Norm = actSol.NormL2();
	Double incrementalErr;

	if (actSolL2Norm)
	  incrementalErr = solIncrL2Norm / actSolL2Norm;
	else
	  {
	    Warning("Zero solution vector!! ", __FILE__,__LINE__);
	    incrementalErr = -EPS; // don't block the iteration loop
	  }
      
	if (actSolL2Norm < NORM_EPS)
	  incrementalErr = actSolL2Norm; // take the absolute error instead of the relative
       



	// =====================================================================
	// output of norms and data
	// =====================================================================


	WriteClaNlNorms(iterationCounter, residualL2Norm, extForcesL2Norm, residualErr,  
			solIncrL2Norm, actSolL2Norm, incrementalErr);
      
	Info->WriteNonLinIter(pdename_, iterationCounter, residualErr, incrementalErr, 
			      etaLineSearch);
      

	// boolean variable, holds condition if another iteration step is necessary
	performOneMoreStep = 
	  (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  
    //perform corrector step  
    TS_alg_->Corrector(solhelp->GetAlgSysVector());
  }



  //   Default Destructor
  // **********************
  SolveStepMech::~SolveStepMech() {

    ENTER_FCN( "SolveStepMech::~SolveStepMech" );
 
  }

} // end of namespace

