#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMag.hh"


namespace CoupledField {

  SolveStepMag::SolveStepMag(BasePDE& apde) : BaseSolveStep(apde)
  {

    ENTER_FCN( "SolveStepMag::SolveStepMag" );
  }

  
  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================


  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================
  void SolveStepMag :: PreStepStatic(const Integer kstep, const Double asteptime,
				     const Integer level, const Boolean reset) {
    ENTER_FCN( "SolveStepMag::PreStepStatic" );
    if (pdeIsCoupled_) 
      algsys_->InitSol();

  }


  void SolveStepMag::PostStepStatic(const Integer level) {
    ENTER_FCN( "SolveStepMag::PostStepStatic" );
    if (pdeIsCoupled_) (*iterCoupledCounter_)++;
  }


  void SolveStepMag::StepStaticNonLin(const Integer kstep, const Double aTime,
                                const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepMag::SolveStepStaticNonLin" );

    const Integer job = 1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;
  
    Vector<Double> solInc(eqnData_->GetNumEQNs());
    Vector<Double> actSol(eqnData_->GetNumEQNs());

    sol_->GetAlgSysVector(actSol);

    SetBCs(level, 0);

    // setup right hand side ==============================================
    Double RhsLinL2Norm = SetLinRHS(level); 

    //add nonlinear part to RHS 
    assemble_->AssembleNLRHS(level);

    do {
      iterationCounter++;

      std::cout << std::endl << "Nonlinear Magnetics: Perform internal loop "
                << "nr. " << iterationCounter << std::endl;      

#ifdef DEBUG
      *debug << std::endl << "=============================================="
             << "========\n"
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
      StoreAlgsysToVec(solInc, algsys_->GetSolutionVal() );
      
      actSol += solInc;
      sol_->SetAlgSysVector(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS(RhsLinVal_.GetPointer());

      assemble_->AssembleNLRHS(level);  // nonlinear part of RHS


      // calculation of residual error (takes care for Dirichlet BCs========
      Vector<Double> actRHS;
      StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
          
      Double residualL2Norm;
      residualL2Norm = RhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )


      // calculation of residual error =======================================
      Double residualErr;
      if (RhsLinL2Norm > 1.0)
        residualErr = residualL2Norm / RhsLinL2Norm;
      else
        residualErr = residualL2Norm;

      // calculate incremental error ========================================
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm  = actSol.NormL2();
      
      Double incrementalErr;      
      if (actSolL2Norm > 1.0)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;

      // =====================================================================
      // output of norms and data
      // =====================================================================
      Double etaLineSearch = 1;
      if ( nonLinLogging_ == TRUE ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                              incrementalErr, etaLineSearch);
      }
      

      // boolean variable, holds condition if another
      // iteration step is necessary
      performOneMoreStep = 
        (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
      
      if (!(performOneMoreStep && iterationCounter < nonLinMaxIter_))
        mycout << "incrementalErr " << incrementalErr << myendl 
               << "incStopCrit_ " << incStopCrit_ << myendl
               << "residualErr " << residualErr  << myendl
               << "residualStopCrit_ " << residualStopCrit_ << myendl;

    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

  }


  void SolveStepMag::StepTransNonLin(const Integer kstep, const Double asteptime,
                               const Integer level, const Boolean reset) {

    ENTER_FCN( "SolveStepMag::StepTransNonLin" );

    lasttimecalc_ = asteptime;
    laststepcalc_ = kstep;

    const Integer job = 1;
  
    static Integer timeStepCounter=1;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;

    Vector<Double> actSol(numPDENodes_);
    Vector<Double> solInc(numPDENodes_);
  
    // Cast BaseStoreSol into StoreSol<Double>,
    // since this function is only called
    // in the transient case
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //set actual solution  
    actSol = solhelp->GetAlgSysVector();

    if ( TS_alg_== NULL)
      Error("TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin");

    //compute predictors
    TS_alg_->Predictor(actSol);

    //now set up RHS: all linear source terms
    Double RhsLinL2Norm = SetLinRHS(level); 

    // inner forces due to nonlin formulation
    assemble_->AssembleNLRHS(level, lasttimecalc_);  

    //Update RHS (mass matrix on right hand side)
    TS_alg_->UpdateRHS(solhelp->GetAlgSysVector());
  
    timeStepCounter++;
    do {
      iterationCounter++;
      std::cout << std::endl
                << "Nonlinear Magnetics: Perform internal loop nr. " 
                << iterationCounter << std::endl;

#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             << "Nonlinear Magnetics: Perform internal loop nr. "
             << iterationCounter << std::endl;      
#endif

      // setup and solve new system (rhs is already set) ====================
      assemble_->InitNonLinMatrices();
      assemble_->AssembleMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      SetBCs(level, lasttimecalc_);

      algsys_->BuildInDirichlet();

      if (job == 1) {
        algsys_->SetupPrecond(job);
        algsys_->SetupSolver(job);
      }

      algsys_->Solve();

      // new solution is only an increment of the full solution =============
      StoreAlgsysToVec(solInc, algsys_->GetSolutionVal() );

      Double residualL2Norm;
      Double etaLineSearch = 0;
      
      if ( lineSearch_ != "no" ) {
        actSol += solInc;
      }
      else {
        Error("Currently no LINE-SEARCH for magnetic PDE");
      }

      //store A_/n+1) in the solution-object sol_
      sol_->SetAlgSysVector(actSol);

      // recalculate RHS with new values to get new residual (f^(k+1))=======
      algsys_->InitRHS(RhsLinVal_.GetPointer());

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(actSol);

      // inner forces due to nonlin formulation
      assemble_->AssembleNLRHS(level, lasttimecalc_);
 

      // ====================================================================
      // calculation of error norms
      // ====================================================================

      if ( lineSearch_ != "no" ) {

        Vector<Double> actRHS;
        StoreAlgsysToVec(actRHS, algsys_->GetRHSVal() );       
          
        // ------------------------------------------------------------------
        // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
        // ------------------------------------------------------------------
        residualL2Norm = RhsL2Norm(actRHS);
      }
      
      Double residualErr;
      if ( RhsLinL2Norm > 1)
        residualErr    = residualL2Norm /  RhsLinL2Norm;
      else
        residualErr    = residualL2Norm;

      // --------------------------------------------------------------------
      // calculate incremental error
      // --------------------------------------------------------------------
      Double solIncrL2Norm = solInc.NormL2();
      Double actSolL2Norm = actSol.NormL2();
      Double incrementalErr;
      
      if (actSolL2Norm > 1)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;

      // --------------------------------------------------------------------
      // output of norms and data
      // --------------------------------------------------------------------
      if ( nonLinLogging_ == TRUE ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualErr,
                              incrementalErr, etaLineSearch);
      }

      // boolean variable, holds condition if another iteration step
      // is necessary
      performOneMoreStep = 
        (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

    //perform corrector step  
    TS_alg_->Corrector(actSol);
  }



  //   Default Destructor
  // **********************
  SolveStepMag::~SolveStepMag() {

    ENTER_FCN( "SolveStepMag::~SolveStepMag" );
 
  }

} // end of namespace

