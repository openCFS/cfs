// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepAcoustic.hh"

#include "assemble.hh"
#include "Forms/linearForm.hh"
#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  SolveStepAcoustic::SolveStepAcoustic(StdPDE& apde, bool justInterpolate)
    : StdSolveStep(apde),
      justInterpolate_(justInterpolate)
  {
  }

  SolveStepAcoustic::~SolveStepAcoustic() {
  }

 
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepAcoustic::StepTransNonLin(PtrParamNode analysis_base) 
  {
    UInt job;
    bool performOneMoreStep;
    UInt iterationCounter=0;
  
    Vector<Double> newSol(numEqns_);
    Vector<Double> oldSol(numEqns_);
  
    // just update dirichlet values
    job = 3;
  
    // if first time step, setup system matrix
    if (actStep_ == 1) {
      assemble_->AssembleMatrices();
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
               
      //set job to 1: build in dirichlet BCs and compute preconditioner
      job = 1;
    }  
  
    //get actual solution
    Vector<Double>  solhelp =
      dynamic_cast<Vector<Double>&>(*(PDE_.GetSolutionVector()));

    //  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    //compute predictors
    TS_alg_->Predictor(solhelp);

    // set BCs, if effective mass matrix formulation, values of BCs depend on 
    //  predictors, so predictors have to be computed beforehand
    PDE_.SetBCs();

    // set old solution  
    newSol = solhelp; //>GetAlgSysVector();

    // Update RHS (mass matrix and damping matrix on right hand side)
    TS_alg_->UpdateRHS();

    //account for RHS-forms
    assemble_->AssembleLinRHS();

    // stores this as linear part of RHS
    algsys_->GetRHSVal( RhsLinVal_ );
    

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
      
      PtrParamNode analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);
      
      // set solution of previous iteration
      oldSol = newSol;
        
      if ( iterationCounter>1 ) 
        job = 3;
        
      // Set linear part of RHS
      algsys_->InitRHS( RhsLinVal_ );
        
      // put nonlinear part to RHS
      AddNonLinRHS();

      algsys_->BuildInDirichlet();
      if ( job == 1 ) {
        algsys_->SetupSolver(analysis_id);
        algsys_->SetupPrecond();
      }

      algsys_->Solve(analysis_id);

      // store new solution in newSol
      algsys_->GetSolutionVal( newSol );
          
      // perform corrector step, if effective mass formulation is used,
      //   we need the Corrector step before we store newsol to sol_,
      //   because newsol is second time derivative at first!
      TS_alg_->Corrector(newSol);

      //put new solution to sol_
      PDE_.SaveSolution( newSol.GetPointer(), newSol.GetSize() );

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
      if ( nonLinLogging_ == true )
        WriteNonLinIterToInfoXML(pdename_, iterationCounter, incrementalErr, incrementalErr);
        
      // boolean variable, holds condition if another iteration step
      //  is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_);
      
    } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

  }


  void SolveStepAcoustic::AddNonLinRHS() {
  

    Vector<Double>     sol, solderiv1, solderiv2, rhs;
    BaseFE             * ptElem;
    StdVector<UInt> connect;
    StdVector<Integer> connect_PDE;
  
    Double coeff1, coeff2;
    Double density, compressibility, c0, BoverA;

    for (UInt actSD=0; actSD<subdoms_.GetSize(); actSD++) {

      std::auto_ptr<LinearForm> rhsInt;
      std::auto_ptr<nLinWesterveltRHSInt> rhsWest;
      std::auto_ptr<nLinKuznetsovRHSInt> rhsKuz;
      
      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[actSD] );
      RegionIdType actRegionId = subdoms_[actSD];

      // get material data
      materialData_[actRegionId]->GetScalar(density,DENSITY,Global::REAL);
      materialData_[actRegionId]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
      c0 = sqrt(compressibility/density);
      materialData_[actRegionId]->GetScalar(BoverA,BOVERA,Global::REAL);

      if ( regionNonLinType_[actRegionId] == WESTERVELT ) {
        rhsWest = std::auto_ptr<nLinWesterveltRHSInt>(new nLinWesterveltRHSInt( isaxi_));

        // set correct factors for bilinear forms
        coeff1 = (1+0.5*BoverA) / pow(c0,4);
        rhsWest->SetFactor( coeff1 );
        rhsInt = rhsWest;
      }
      if ( regionNonLinType_[actRegionId] == KUZNETSOV ) {
        rhsKuz = std::auto_ptr<nLinKuznetsovRHSInt>(new nLinKuznetsovRHSInt( isaxi_ ));

        // set correct factors for bilinear forms
        coeff1 = density * BoverA / pow(c0,4);
        rhsKuz->SetFactor( coeff1 );
        coeff2 = density * 2.0 / (c0*c0);
        rhsKuz->SetSecondFactor( coeff2 );
        rhsInt = rhsKuz;

      }

      // In case there is something nonlinear go on
      if ( (regionNonLinType_[actRegionId] == WESTERVELT)
           || (regionNonLinType_[actRegionId] == KUZNETSOV) ) {

        // Get iterator
        EntityIterator it = actSDList.GetIterator();
        for ( it.Begin(); !it.IsEnd(); it++) {
          
          ptElem  = it.GetElem()->ptElem;
          connect = it.GetElem()->connect;
          
          PDE_.GetSolVecOfElement(sol, it, results_[0]);
          PDE_.GetDerivSolVecOfElement(solderiv1, it, results_[0]);
          PDE_.GetDeriv2SolVecOfElement(solderiv2, it, results_[0]);
          
          rhsInt->SetActElemSol(sol);
          rhsInt->SetActElemSolDeriv1(solderiv1);
          rhsInt->SetActElemSolDeriv2(solderiv2);
          rhsInt->CalcElemVector(rhs, it);
          
          //get equation numbers 
          
          eqnMap_->GetNodeEqn( connect, connect_PDE );
          
          //assemble
          algsys_->SetElementRHS( rhs, pdeId1_, connect_PDE );
        }
      }
    }
  }

  void SolveStepAcoustic::StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) {
    if(justInterpolate_) {
      //account for RHS
      assemble_->AssembleLinRHS(adjointParams);
      PDE_.ComputeRHS( actTime_ );
      
      // store rhs vector back to PDE
      Vector<Double> rhs;
      algsys_->GetRHSVal(rhs);
      PDE_.SaveRHS( rhs.GetPointer(), rhs.GetSize());
    }
    else {
      StdSolveStep::StepTransLin(analysis_id, adjointParams);
    }
  }


} // end of namespace

