#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcoustic.hh"
#include "Forms/forms_header.hh"

#include "assemble.hh"

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

  void SolveStepAcoustic::StepTransNonLin( const Boolean reset ) {

    ENTER_FCN( "SolveStepAcoustic::StepTransNonLin" );

    Double *solPtr;
  
    UInt job;
    Boolean performOneMoreStep;
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
    SetBCs(actTime_);

    // set old solution  
    newSol = solhelp->GetAlgSysVector();

    // Update RHS (mass matrix and damping matrix on right hand side)
    TS_alg_->UpdateRHS();

    // stores this as linear part of RHS
    algsys_->GetRHSVal( solPtr );
    StoreAlgsysToVec(RhsLinVal_, solPtr );

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
      if (iterationCounter == 1)
        std::cout << std::endl << "Time step:   "  << actStep_ 
                  << "  ,Iterations: " << iterationCounter;
      else 
        std::cout     << "  " << iterationCounter;

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
      if ( nonLinLogging_ == TRUE ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, incrementalErr,
                              incrementalErr);
      }
        
      // boolean variable, holds condition if another iteration step
      //   is necessary
      performOneMoreStep =    (incrementalErr > incStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

  }


  void SolveStepAcoustic::AddNonLinRHS() {
  
    ENTER_FCN( "SolveStepAcoustic::AddNonLinRHS" );

    Matrix<Double>     ptCoord;
    Vector<Double>     sol, solderiv1, solderiv2, rhs;
    BaseFE             * ptElem;
    StdVector<UInt> connect;
    StdVector<Integer> connect_PDE;
  
    Double coeff1, coeff2;
    Double density, compressibility, c0, BoverA;

    BaseForm * rhsInt = NULL;

    if ( subdoms_.GetSize() != nonLinPDEName_.GetSize() )
      Error("nonLinPDEName_ does not match size of subdoms_",__FILE__,__LINE__);

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {

      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd,subdoms_[actSD]);
        
      // get material data
      density         = materialData_[actSD].GetDensity();
      compressibility = materialData_[actSD].GetCompressibility();
      c0 = sqrt(compressibility/density);
      BoverA = materialData_[actSD].GetBoverA();

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
        
      for (UInt j=0; j < elemssd.GetSize(); j++) {

        ptElem  = elemssd[j]->ptElem;
        connect = elemssd[j]->connect;
        GetElemCoords(connect, ptCoord);

        GetSolVecOfElement(sol, connect);
        GetDerivSolVecOfElement(solderiv1, connect);
        GetDeriv2SolVecOfElement(solderiv2, connect);

        rhsInt->SetElemPtr(ptElem);

        rhsInt->SetActElemSol(sol);
        rhsInt->SetActElemSolDeriv1(solderiv1);
        rhsInt->SetActElemSolDeriv2(solderiv2);
        rhsInt->CalcElemVector(ptCoord, rhs);
          
        //get equation numbers 
        eqnData_->Node2EQN(connect, connect_PDE);

        //assemble
        algsys_->SetElementRHS(&rhs[0], pdeId1_, connect_PDE.GetPointer(), 
                               connect_PDE.GetSize());
      }
    }
    delete rhsInt;
  }


} // end of namespace

