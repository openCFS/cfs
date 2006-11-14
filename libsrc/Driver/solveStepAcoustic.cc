#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcoustic.hh"

#include "assemble.hh"
#include "Forms/linearForm.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

namespace CoupledField {

  SolveStepAcoustic::SolveStepAcoustic(StdPDE& apde) : StdSolveStep(apde) 
  {
    ENTER_FCN( "SolveStepAcoustic::SolveStepAcoustic" );
  }

  SolveStepAcoustic::~SolveStepAcoustic() {
    ENTER_FCN( "SolveStepAcoustic::~SolveStepAcoustic" );
  }

 
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepAcoustic::StepTransNonLin() {

    ENTER_FCN( "SolveStepAcoustic::StepTransNonLin" );

    Double *solPtr;
  
    UInt job;
    bool performOneMoreStep;
    UInt iterationCounter=0;
  
    Vector<Double> newSol(numPDENodes_);
    Vector<Double> oldSol(numPDENodes_);
  
    // just update dirichlet values
    job = 3;
  
    // if first time step, setup system matrix
    if (actStep_ == 1) {
      assemble_->AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
               
      //set job to 1: build in dirichlet BCs and compute preconditioner
      job = 1;
    }  
  
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //compute predictors
    TS_alg_->Predictor(solhelp->GetAlgSysVector());

    // set BCs, if effective mass matrix formulation, values of BCs depend on 
    //  predictors, so predictors have to be computed beforehand
    PDE_.SetBCs( actTime_ );

    // set old solution  
    newSol = solhelp->GetAlgSysVector();

    // Update RHS (mass matrix and damping matrix on right hand side)
    TS_alg_->UpdateRHS();

    // stores this as linear part of RHS
    algsys_->GetRHSVal( solPtr );
    StoreAlgsysToVec(RhsLinVal_, solPtr );

    do {
      iterationCounter++;
      // write out number of iteration loops to standard out
	  if ( (actStep_ < 100) || (actStep_%1000) == 0 ) {
		if (iterationCounter == 1)
		  std::cout << std::endl << "Time step:   "  << actStep_ 
					<< "  ,Iterations: " << iterationCounter;
		else 
		  std::cout << "  " << iterationCounter;
	  }

#ifdef DEBUG
      *debug << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear Acoustics: Perform internal loop no. "
             << iterationCounter << std::endl;      
#endif
        
      // set solution of previous iteration
      oldSol = newSol;
        
      if ( iterationCounter>1 ) 
        job = 3;
        
      // Set linear part of RHS
      algsys_->InitRHS(RhsLinVal_.GetPointer());
        
      // put nonlinear part to RHS
      AddNonLinRHS();

      algsys_->BuildInDirichlet();
      if ( job == 1 ) {
        algsys_->SetupSolver();
        algsys_->SetupPrecond();
      }

      algsys_->Solve();

      // store new solution in newSol
      algsys_->GetSolutionVal( solPtr );
      StoreAlgsysToVec(newSol, solPtr );
          
      // perform corrector step, if effective mass formulation is used,
      //   we need the Corrector step before we store newsol to sol_,
      //   because newsol is second time derivative at first!
      TS_alg_->Corrector(newSol);

#ifdef DEBUG
      *debug << std::endl
             << "New Solution:" << std::endl << newSol << std::endl;
#endif

      //put new solution to sol_
      sol_->SetAlgSysVector(newSol);  

      // compute L2-Norm of error between last incremental solution and
      //   actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newSol.GetSize(); i++)
        solIncrL2Norm += (newSol[i]-oldSol[i])*(newSol[i]-oldSol[i]);
        
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newSol.NormL2();

      Double incrementalErr;
      if (actSolL2Norm > 1)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;
        
      // output of norms and data
      if ( nonLinLogging_ == true ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, incrementalErr,
                              incrementalErr);
      }
        
      // boolean variable, holds condition if another iteration step
      //  is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

  }


  void SolveStepAcoustic::AddNonLinRHS() {
  
    ENTER_FCN( "SolveStepAcoustic::AddNonLinRHS" );

    Vector<Double>     sol, solderiv1, solderiv2, rhs;
    BaseFE             * ptElem;
    StdVector<UInt> connect;
    StdVector<Integer> connect_PDE;
  
    Double coeff1, coeff2;
    Double density, compressibility, c0, BoverA;

    LinearForm * rhsInt = NULL;

    if ( subdoms_.GetSize() != nonLinPDEName_.GetSize() )
      Error("nonLinPDEName_ does not match size of subdoms_",__FILE__,__LINE__);

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[actSD] );
      
      // get material data
      materialData_[actSD]->GetScalar(density,DENSITY,REAL);
      materialData_[actSD]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
      c0 = sqrt(compressibility/density);
      materialData_[actSD]->GetScalar(BoverA,BOVERA,REAL);

      if ( nonLinPDEName_[actSD] == WESTERVELT ) {
        rhsInt = new nLinWesterveltRHSInt(1.0, isaxi_);

        // set correct factors for bilinear forms
        coeff1 = (1+0.5*BoverA) / pow(c0,4);
        rhsInt->SetFactor(coeff1);      
      }
      if ( nonLinPDEName_[actSD] == KUZNETSOV ) {
        rhsInt = new nLinKuznetsovRHSInt(1.0, isaxi_);

        // set correct factors for bilinear forms
        coeff1 = density * BoverA / pow(c0,4);
        rhsInt->SetFactor(coeff1);
        coeff2 = density * 2.0 / (c0*c0);
        rhsInt->SetSecondFactor(coeff2);
      }
        
      // Get iterator
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++) {
        
        ptElem  = it.GetElem()->ptElem;
        connect = it.GetElem()->connect;

        PDE_.GetSolVecOfElement(sol, it);
        PDE_.GetDerivSolVecOfElement(solderiv1, it);
        PDE_.GetDeriv2SolVecOfElement(solderiv2, it);
        
        rhsInt->SetActElemSol(sol);
        rhsInt->SetActElemSolDeriv1(solderiv1);
        rhsInt->SetActElemSolDeriv2(solderiv2);
        rhsInt->CalcElemVector(rhs, it);
        
        //get equation numbers 
        
        eqnMap_->GetNodeEqn( connect, connect_PDE );

        //assemble
        algsys_->SetElementRHS(&rhs[0], pdeId1_, connect_PDE.GetPointer(), 
                               connect_PDE.GetSize());
      }
    }
    delete rhsInt;
  }


} // end of namespace

