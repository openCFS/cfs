// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "OLAS/algsys/basesystem.hh"
#include "solveStepMicroPiezo.hh"
#include "Utils/preisach.hh"
#include "Utils/jiles.hh"
#include "assemble.hh"
#include "Forms/gradfieldop.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "Domain/domain.hh"
#include "Utils/result.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "Utils/piezoMicroModel.hh"
#include "Utils/piezoMicroModelBK.hh"

namespace CoupledField {

  SolveStepMicroPiezo::SolveStepMicroPiezo(StdPDE& apde) : StdSolveStep(apde)
  {

    doInit_ = true;
  }


  SolveStepMicroPiezo::~SolveStepMicroPiezo() {
  }
 

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMicroPiezo:: PreStepTrans()
  {


    // due to coupling-pdes, the RHS has to be initialized BEFORE 
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();


    // save solution of time step n as previous solution
    //save solution of previous time step 
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

    //store solution for (n)
    PDE_.SavePrevSolution( solHelp.GetPointer(), solHelp.GetSize() );

  }



  void SolveStepMicroPiezo::SolveStepTrans(InfoNode* analysis_id) {

    StepTransNonLin(analysis_id);

  }


  void SolveStepMicroPiezo::StepTransNonLin(InfoNode* analysis_id) {


    //std::cout << "\n In :StepTransNonLinEpsDiff  \n " << std::endl;
 
    bool performOneMoreStep;
    UInt iterationCounter=0;
  
    Vector<Double> newSol( numEqns_ ); 
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> solPrev( numEqns_ );
 
    // get second order derivative of previous time step;
    Vector<Double> solDeriv2Prev = PDE_.getS2();

    Vector<Double> coeff;

    oldSol.Init(0);
    newSol.Init(0);
  
    //save solution of previous time step  
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    solPrev = solHelp;
    //    std::cout << "Previous Solution:\n" << solPrev << std::endl;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Exception( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin");
    }
    else {
      //compute predictors
      TS_alg_->Predictor(solPrev);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( );    //clear RHS

    do {
      iterationCounter++;
      // for every time step write out number of iteration loops to standard out
//       if (iterationCounter == 1)
//         std::cout << std::endl << "Time step:   "  << actStep_
//                   << "  ,Iterations: " << iterationCounter << std::endl;
//       else 
//         std::cout << "Iter:  " << iterationCounter << std::endl;
    

#ifdef DEBUG
      std::cout << std::endl
             << "====================================================== "
             << std::endl
             <<   "Nonlinear MicroPiezo: Perform internal loop no. "
             << iterationCounter << std::endl;      
#endif
        

      // set solution of previous iteration
      if (iterationCounter == 1 ) {
        oldSol = solPrev;
      }
      else {
        oldSol = newSol;
      }
    
      //clear RHS
      algsys_->InitRHS();

      // setup RHS to incorporate loads and linear-Forms
      assemble_->AssembleLinRHS(); 

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS();

      //perform new assembly
      assemble_->AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      // K*(u_n,Ve_n)
      //algsys_->UpdateRHS(STIFFNESS, solPrev.GetPointer());
      algsys_->UpdateRHS(STIFFNESS, solPrev);

      // M*(d2u_n,Ve_n)
      //algsys_->UpdateRHS(MASS, solDeriv2Prev.GetPointer());
      algsys_->UpdateRHS(MASS, solDeriv2Prev);

      // build in the Dirichlet vales in system mmatrix and rhs
      algsys_->BuildInDirichlet();

      //get RHS
      Vector<Double> RHS;
      algsys_->GetRHSVal(RHS);
      Double residualNorm = PDE_.GetRhsL2Norm( RHS );
      residualNorm = 0.0;

      algsys_->SetupSolver(analysis_id);
      algsys_->SetupPrecond();
    
      algsys_->Solve(analysis_id);
      algsys_->GetSolutionVal(newSol); 

      //store solution for (n+1)
      PDE_.SaveSolution( newSol.GetPointer(), newSol.GetSize() );

      // compute L2-Norm of error between last incremental solution and
      // actual incremental solution
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
      nonLinLogging_ = true;
      if ( nonLinLogging_ == true ) {
        Info->WriteNonLinIter(pdename_, iterationCounter, residualNorm,
                              incrementalErr);
      }

      //      std::cout << "Norm=" << incrementalErr << std::endl << std::endl;
    
      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

//     if ( iterationCounter >= nonLinMaxIter_ ) {
//       (*error) << "Number of nonlinear iterations exceeds limit "
//                << "nonLinearMaxIter_ = "
//                << nonLinMaxIter_;
//       Error( __FILE__, __LINE__ );
//     }

    //set the current values to the previous for the next time step!
    SetPreviousVals();

    //perform corrector step  
    TS_alg_->Corrector(newSol);
  }


  void SolveStepMicroPiezo::SetPreviousVals() {
  
    // get solution object of PDE elec
    StdPDE * pdeDirectCoupledPiezo = domain->GetStdPDE( "electrostatic-mechanic" );

    StdVector<BasePairCoupling*> directCouplings =   
      *(pdeDirectCoupledPiezo->GetCouplingsObject());


    std::map<RegionIdType, BaseMaterial*> piezoMat = 
      directCouplings[0]->getPDEMaterialData();


//    PiezoMicroModelBK* microPiezo = piezoMat[0]->GetMicroPiezoModel(); unused, Fabian

    BaseMaterial * actSDMat = NULL;

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = piezoMat.begin(); it != piezoMat.end(); it++ ) {
      actSDMat = it->second;

      PiezoMicroModelBK* microPiezo = actSDMat->GetMicroPiezoModel();
      if (  microPiezo!= NULL ) {
        microPiezo->SetPreviousVolFrac();
        microPiezo-> SetPreviousIrreversibleValues();
      }
    }
  
  }


} // end of namespace

