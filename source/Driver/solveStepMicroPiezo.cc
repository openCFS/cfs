// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "OLAS/algsys/algebraicSys.hh"
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
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Utils/piezoMicroModel.hh"
#include "Utils/piezoMicroModelBK.hh"
#include "DataInOut/ResultCache.hh"

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
    ResultCache::SetStepValue( actTime_ );

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



  void SolveStepMicroPiezo::SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams) {

    StepTransNonLin(analysis_id);

  }


  void SolveStepMicroPiezo::StepTransNonLin(PtrParamNode analysis_id) {


    //    std::cout << "\n In :StepTransNonLinEpsDiff  \n " << std::endl;
 
    bool performOneMoreStep;
    UInt iterationCounter=0;
    Double etaLineSearch, residualNorm;
    Double incrementalErrPrev = 1e20;
    Double incrementalErr;

    Vector<Double> newSol( numEqns_ ); 
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> solPrev( numEqns_ );
    //    Vector<Double> incrSol( numEqns_ );
    //    Vector<Double> oldSolMech, newSolMech;

    // get second order derivative of previous time step;
    Vector<Double> solDeriv2Prev = PDE_.getDeriv(SECOND_DERIV);

    Vector<Double> coeff;

    oldSol.Init(0);
    newSol.Init(0);
    //    incrSol.Init(0);
  
    //save solution of previous time step  
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
    solPrev = solHelp;

    //get just mechanical solution
    DirectCoupledPDE& coupledPDE = dynamic_cast<DirectCoupledPDE&>(PDE_);
    StdVector<SinglePDE*> allSingelPDEs =  coupledPDE.GetSinglePDEs();
    BaseNodeStoreSol * solPDEMech =allSingelPDEs[1]->getPDESolution();
    NodeStoreSol<Double> *solhelp1 =
      dynamic_cast<NodeStoreSol<Double>*>(solPDEMech);

    Vector<Double> oldMech;
    solhelp1->GetGlobalSolVector(MECH_DISPLACEMENT, oldMech);

    //    StdPDE& mechPDE =  dynamic_cast<StdPDE&>(*allSingelPDEs[0]); 
    //    Vector<Double>& oldSolMech =  dynamic_cast<Vector<Double>&>(haha); 

    //    std::cout << "Oldmech:\n" << oldMech << "\n\n" << std::endl;

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

#ifndef NDEBUG
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

      //Vector<Double> RHS;
      //algsys_->GetRHSVal(RHS);
      //std::cout << "RHS:\n " << RHS << std::endl;

      //perform new assembly
      assemble_->AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);


      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS();

      // K*(u_n,Ve_n)
      algsys_->UpdateRHS(STIFFNESS, solPrev);
  
      // M*(d2u_n,Ve_n)
      algsys_->UpdateRHS(MASS, solDeriv2Prev);

      // build in the Dirichlet vales in system matrix and rhs
      algsys_->BuildInDirichlet();

      algsys_->SetupSolver(analysis_id);
      algsys_->SetupPrecond();
    
      algsys_->Solve(analysis_id);
      algsys_->GetSolutionVal(newSol); 
      //std::cout << "new Sol:\n" << newSol << std::endl;

//       //perform a line search
//       LineSearchPM( newSol, oldSol, solPrev, solDeriv2Prev, etaLineSearch, 
//                     residualNorm );

      //store solution for (n+1)
      PDE_.SaveSolution( newSol.GetPointer(), newSol.GetSize() );

      BaseNodeStoreSol * solPDEMechNew =allSingelPDEs[1]->getPDESolution();
      NodeStoreSol<Double> *solhelpNew =
        dynamic_cast<NodeStoreSol<Double>*>(solPDEMechNew);

      Vector<Double> newMech;
      solhelpNew->GetGlobalSolVector(MECH_DISPLACEMENT, newMech);

      // compute L2-Norm of error between last incremental solution and
      // actual incremental solution
      Double solIncrL2Norm=0;
      for (UInt i=0; i<newMech.GetSize(); i++)
        solIncrL2Norm += (newMech[i]-oldMech[i])*(newMech[i]-oldMech[i]);
    
      solIncrL2Norm = sqrt(solIncrL2Norm);
      Double actSolL2Norm = newMech.NormL2();
    
      if (actSolL2Norm > 1)
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      else
        incrementalErr = solIncrL2Norm;
 
      etaLineSearch = 1.0;
//       if (  incrementalErr > incrementalErrPrev ) {
//         //perform a line search

//         LineSearchPM( newSol, oldSol, solPrev, solDeriv2Prev, etaLineSearch, 
//                       residualNorm );

//         //store solution for (n+1)
//         PDE_.SaveSolution( newSol.GetPointer(), newSol.GetSize() );

//         solIncrL2Norm=0;
//         for (UInt i=0; i<newMech.GetSize(); i++)
//           solIncrL2Norm += (newMech[i]-oldMech[i])*(newMech[i]-oldMech[i]);
    
//         solIncrL2Norm = sqrt(solIncrL2Norm);
//         actSolL2Norm = newMech.NormL2();
    
//         if (actSolL2Norm > 1)
//           incrementalErr = solIncrL2Norm / actSolL2Norm;
//         else
//           incrementalErr = solIncrL2Norm;

//       }


      incrementalErrPrev = incrementalErr;
      residualNorm = incrementalErr;

      // output of norms and data
      nonLinLogging_ = true;
      if ( nonLinLogging_ == true )
        WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualNorm, incrementalErr, etaLineSearch);

      //      std::cout << "Norm=" << incrementalErr << std::endl << std::endl;
    
      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
//       if ( iterationCounter < nonLinMaxIter_ -1 ) 
//         performOneMoreStep = true;

      oldMech = newMech;

    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  

//     if ( iterationCounter >= nonLinMaxIter_ ) {
//       EXCEPTION("Number of nonlinear iterations exceeds limit "
//                 << "nonLinearMaxIter_ = "
//                 << nonLinMaxIter_);
//     }

    if ( incrementalErr > incStopCrit_ ) 
      std::cout << "Not converged, norm is: " <<incrementalErr << std::endl; 

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
        microPiezo->SetPreviousIrreversibleValues();
      }
    }
  
  }


  void SolveStepMicroPiezo::LineSearchPM( Vector<Double>& newSol, 
                                          Vector<Double>& oldSol,
                                          Vector<Double>& solPrev,
                                          Vector<Double>& solDeriv2Prev,
                                          Double& etaLineSearch, 
                                          Double& residualL2NormOpt )
  {

    Vector<Double> solIncrement = newSol;
    solIncrement -= oldSol;

    const UInt nrEtas = 4;
    const Double eta[nrEtas] = {1.0, 0.5, 0.3, 0.1 };
    //    const Double eta[nrEtas] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0};
		// initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    residualL2NormOpt = 1e15;

    Vector<Double> actSol;
    for( UInt i=0; i<nrEtas; i++) {
      actSol = oldSol + solIncrement * eta[i];
      
      //store new solution
      PDE_.SaveSolution( actSol.GetPointer(), actSol.GetSize() );

      // Recalculate residual, f-Cu-Mu-K*u
      algsys_->InitRHS();

      assemble_->AssembleLinRHS();

      // assemble!
      assemble_->AssembleMatrices();

      // time step scheme
      TS_alg_->UpdateRHS(actSol);

      // K*(u_n,Ve_n)
      algsys_->UpdateRHS(STIFFNESS, solPrev);
      
      // M*(d2u_n,Ve_n)
      algsys_->UpdateRHS(MASS, solDeriv2Prev);


      // substract K^* u^k from RHS
      //      TS_alg_->SubstractStiffnessFromRHS(actSol);
      Vector<Double> actSolMinus = actSol;
      actSolMinus *= -1.0;
      algsys_->UpdateRHS(STIFFNESS, actSolMinus);

      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  (f-Ku )

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;

    actSol  = oldSol + solIncrement * etaOpt;

  }


} // end of namespace

