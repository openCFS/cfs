#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMech.hh"
#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

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

  void SolveStepMech:: PreStepStatic( const Boolean reset )
  {
    ENTER_FCN( "SolveStepMech::PreStepStatic" );

    if (isIterCoupled_)
      // init RHS at this place, because forces of other PDEs are added to RHS afterwards
      algsys_->InitRHS();     
    else {
      algsys_->InitRHS();        
      // if PDE is coupled, the solution of the prior outer loops must be kept
      algsys_->InitSol();
    }
 
  }


  void SolveStepMech::StepStaticNonLin( const Boolean reset )
  {
    ENTER_FCN( "SolveStepMech::SolveStepStaticNonLin" );

    Boolean performOneMoreStep;
    UInt iterationCounter=0;
    NodeStoreSol<Double>  & solHelp = dynamic_cast<NodeStoreSol<Double>&>(*sol_);
  
    Vector<Double>  actSol = solHelp.GetAlgSysVector();

     Vector<Double> solIncrement;
    solIncrement.Resize(eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN());

    PDE_.SetBCs(0);

    // store linear part of RHS
    Double loadFactor = 1.0;
    Double extForcesL2Norm = SetLinRHS(loadFactor); 

    PDE_.AssembleNLRHS();

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
	PDE_.InitNonLinMatrices();
	PDE_.AssembleMatrices();
      
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
	algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver();

	algsys_->Solve();   

	// new solution is only an increment of the full solution =============
	Double *solPtr;
	algsys_->GetSolutionVal( solPtr );
	StoreAlgsysToVec(solIncrement, solPtr);

	Double residualL2Norm = 0;
	Double etaLineSearch=0;

	if ( lineSearch_ != "no" )
	  actSol += solIncrement;
	else
	  // TRUE is for transient simulation
	  residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch);
      
	sol_->SetAlgSysVector(actSol);

	// recalculate RHS with new values to get new residual (f^(k+1))========
	algsys_->InitRHS(RhsLinVal_.GetPointer());

	PDE_.AssembleNLRHS();  // inner forces due to nonlin formulation


	// =====================================================================
	// calculation of error norms
	// =====================================================================
	if ( lineSearch_ != "no" )
	  {
	    Vector<Double> actRHS;
	    algsys_->GetRHSVal(solPtr);
	    StoreAlgsysToVec(actRHS, solPtr );       
          
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


  void SolveStepMech :: PostStepStatic()
  {
    ENTER_FCN( "SolveStepMech::PostStepStatic" );

    if (isIterCoupled_)
      (*iterCoupledCounter_)++;

  }



  void SolveStepMech::StepTransNonLin( const Boolean reset )
  {
    ENTER_FCN( "SolveStepMech::StepTransNonLin" );

    static UInt timeStepCounter=1;
    Boolean performOneMoreStep;
    UInt iterationCounter=0;

    Vector<Double> actSol;
    Vector<Double> solIncrement;
    Double *solPtr;
  
    //is already done in PreStepTrans (basePDE.cc!!)
    //  algsys_->InitRHS();

    // Cast BaseStoreSol into StoreSol<Double>,
    // since this function is only called
    // in the transient case
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  

    actSol = solhelp->GetAlgSysVector();
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    //! store linear part of RHS
    Double loadFactor = 1.0;
    Double extForcesL2Norm = SetLinRHS(loadFactor);

    timeStepCounter++;

    //update RHS due to time stepping
    TS_alg_->UpdateRHS(actSol);

    //! account for Dirichlet BCs
    SetBCs( actTime_ );

  
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
	PDE_.InitNonLinMatrices();
	PDE_.AssembleMatrices();
	assemble_->AssembleMatrices();
	algsys_->ConstructEffectiveMatrix(matrix_factor_);

	TS_alg_->UpdateRHS(actSol);
	SetBCs( actTime_ );

	algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver();

	algsys_->Solve();

	// new solution is only an increment of the full solution =============
	algsys_->GetSolutionVal( solPtr );
	StoreAlgsysToVec( solIncrement, solPtr );

	Double residualL2Norm;
	Double etaLineSearch = 0;
      
	if ( lineSearch_ != "no" )
	  actSol += solIncrement;
	else
	  // TRUE is for transient simulation
	  residualL2Norm = LineSearch(solIncrement, actSol, etaLineSearch, TRUE);

	sol_->SetAlgSysVector(actSol);

	// recalculate RHS with new values to get new residual (f^(k+1))========
	algsys_->InitRHS();
	assemble_->AssembleSrcRHS(actTime_); 
	assemble_->AssembleNLRHS(actTime_);  // inner forces due to nonlin formulation
	TS_alg_->UpdateRHS(actSol);




	// =====================================================================
	// calculation of error norms
	// =====================================================================

	if ( lineSearch_ != "no" )
	  {

	    // TODO: Replace the get of the RHS by getting directly the norm of the RHS
	    Vector<Double> actRHS;
	    Double *solPtr;
	    algsys_->GetRHSVal( solPtr );  
	    StoreAlgsysToVec(actRHS, solPtr );       
          
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

