#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcoustic.hh"
#include "Forms/forms_header.hh"

namespace CoupledField {

  SolveStepAcoustic::SolveStepAcoustic(BasePDE& apde) : BaseSolveStep(apde)
  {

    ENTER_FCN( "SolveStepAcoustic::SolveStepAcoustic" );
  }

  
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepAcoustic::StepTransNonLin(const Integer kstep, const Double asteptime,
				    const Integer level, const Boolean reset) {

    ENTER_FCN( "SolveStepAcoustic::StepTransNonLin" );

    laststepcalc_ = kstep;

    Integer job;
    Boolean performOneMoreStep;
    Integer iterationCounter=0;
  
    Vector<Double> newSol(numPDENodes_);
    Vector<Double> oldSol(numPDENodes_);
  
    // just update dirichlet values
    job = 3;
  
    // if first time step, setup system matrix
    if (laststepcalc_ == 1) {
      assemble_->AssembleMatrices(level);
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
        
      // set BCs 
      SetBCs(level, lasttimecalc_);
        
      //set job to 1: build in dirichlet BCs and compute preconditioner
      job = 1;
    }  
  
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //compute predictors
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    //set old solution  
    newSol = solhelp->GetAlgSysVector();

    //Update RHS (mass matrix on right hand side)
    TS_alg_->UpdateRHS(solhelp->GetAlgSysVector());
  
    // stores linear part of RHS
    StoreAlgsysToVec(RhsLinVal_, algsys_->GetRHSVal() );

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
	std::cout << std::endl << std::endl << "Time step:     "  << kstep 
		  << "  , internal loop no. " << iterationCounter;
      else 
	std::cout << std::endl << "                  " 
		  << "  , internal loop no. " << iterationCounter;

#ifdef DEBUG
      *debug << std::endl
	     << "====================================================== "
	     << std::endl
	     <<   "Nonlinear Acoustics: Perform internal loop no. "
	     << iterationCounter << std::endl;      
#endif
        
      //set solution of previous iteration
      oldSol = newSol;
        
      if ( iterationCounter>1 ) 
	job = 3;
        
      // Set linear part of RHS
      algsys_->InitRHS(RhsLinVal_.GetPointer());
        
        
      //put nonlinear part to RHS
      AddNonlinearRHS();

      algsys_->BuildInDirichlet();
      if ( job == 1 ) {
	algsys_->SetupSolver(job);
	algsys_->SetupPrecond(job);
      }

      algsys_->Solve();
        
      // store new solution in newSol
      StoreAlgsysToVec(newSol, algsys_->GetSolutionVal() );

      //put new solution to sol_
      sol_->SetAlgSysVector(newSol);  

#ifdef DEBUG
      *debug << std::endl
	     << "New Solution:" << std::endl << newSol << std::endl;
#endif

      //perform corrector step  
      TS_alg_->Corrector(newSol);

      // compute L2-Norm of error between last incremental solution and
      //   actual incremental solution
      Double solIncrL2Norm=0;
      for (Integer i=0; i<newSol.GetSize(); i++)
	solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
        
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newSol.NormL2();

      Double incrementalErr;
      if (actSolL2Norm > 1)
	incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
	incrementalErr = solIncrL2Norm;
        
      // output of norms and data
      if ( nonLinLogging_ == TRUE ) {
	Info->WriteNonLinIter(pdename_, iterationCounter, incrementalErr,
			      incrementalErr);
      }
        
      // boolean variable, holds condition if another iteration step
      //   is necessary
      performOneMoreStep =    (incrementalErr > incStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
  

  }

  void SolveStepAcoustic::AddNonlinearRHS() {
  
    ENTER_FCN( "SolveStepAcoustic::AddNonlinearRHS" );
  
    Double coeffN1 = 1.0, coeffN2 = 1.0;
    Double density, compressibility, c0, BoverA;

    BaseForm * rhsInt = new  nLinKuznetsovRHSInt(coeffN1, isaxi_);
    //   BaseForm * N1 = new nLinAcoustic1(coeffN1, isaxi_);
    //   BaseForm * N2 = new nLinAcoustic2(coeffN2, isaxi_);

    Matrix<Double> ptCoord, elemmat;
    Vector<Double> sol, solderiv1, solderiv2, rhs, rhsElem;
    BaseFE         * ptElem;
    StdVector<Integer> connect, Eqns;  
    Integer numEl=0;

    for (Integer actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      StdVector<Elem*> elemssd;
      ptgrid_->GetElemSD(elemssd,subdoms_[actSD],actlevel_);
        
      //set correct factors for bilinear forms
      density         = materialData_[actSD].GetDensity();
      compressibility = materialData_[actSD].GetCompressibility();
      c0 = sqrt(compressibility/density);
      BoverA = materialData_[actSD].GetBoverA();

      coeffN1 = density * BoverA / pow(c0,4);
      rhsInt->SetFactor(coeffN1);
      //      N1->SetFactor(coeffN1);
        
      coeffN2 = density * 2.0 / (c0*c0);
      rhsInt->SetSecondFactor(coeffN2);
      //      N2->SetFactor(coeffN1);
        
      for (Integer j=0; j < elemssd.GetSize(); j++) {
	ptElem  = elemssd[j]->ptElem;
	connect = elemssd[j]->connect;
          
	GetElemCoords(connect, ptCoord, actlevel_);
          
	GetSolVecOfElement(sol, connect);
	GetDerivSolVecOfElement(solderiv1, connect);
	GetDeriv2SolVecOfElement(solderiv2, connect);


	//        sol[0] = 1.00000000000000e-4;
	//        sol[1] = 0;
	//        sol[2] = 0;
	//        sol[3] = 1.00000000000000e-4;
	//        solderiv1[0] = 20000.00;
	//        solderiv1[1] = 0;
	//        solderiv1[2] = 0;
	//        solderiv1[3] = 20000.00;
	//        solderiv2[0] = 4000000000000.00;
	//        solderiv2[1] = 0;
	//        solderiv2[2] = 0;
	//        solderiv2[3] = 4000000000000.00;

	rhsInt->SetElemPtr(ptElem);
	rhsInt->SetActElemSol(sol);
	rhsInt->SetActElemSolDeriv1(solderiv1);
	rhsInt->SetActElemSolDeriv2(solderiv2);
	rhsInt->CalcElemVector(ptCoord, rhs);
          
	//        N1->SetElemPtr(ptElem);
	//        N1->SetActElemSolDeriv1(solderiv1);
	//        N1->CalcElementMatrix(ptCoord, elemmat);
	// #ifdef DEBUG
	//        *debug << std::endl
	//                       << "N1-Matrix: " << std::endl << elemmat << std::endl;
	// #endif
	//        rhs = elemmat * solderiv2;
	//        N2->SetElemPtr(ptElem);
	//        N2->SetActElemSol(sol);
	//        N2->CalcElementMatrix(ptCoord, elemmat);
	// #ifdef DEBUG
	//        *debug << std::endl
	//                       << "N2-Matrix: " << std::endl << elemmat << std::endl;
	// #endif
	//        rhs += rhsElem;

	//get equation numbers 
	StdVector<Integer> connect_PDE;
	eqnData_->Node2EQN(connect, connect_PDE);

	//        std::cerr << "Sol:\n" << sol << std::endl;
	//        std::cerr << "SolD1:\n" << solderiv1 << std::endl;
	//        std::cerr << "SolD2:\n" << solderiv2 << std::endl;
	//        std::cerr << "RHS:\n" << rhs << std::endl;

	//assemble
	algsys_->SetElementRHS(&rhs[0], connect_PDE.GetPointer(), connect_PDE.GetSize());
      }  
    }

  }



  //   Default Destructor
  // **********************
  SolveStepAcoustic::~SolveStepAcoustic() {

    ENTER_FCN( "SolveStepAcoustic::~SolveStepAcoustic" );
 
  }

} // end of namespace

