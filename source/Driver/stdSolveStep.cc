// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <math.h>
#include <iostream>
#include <string>
#include <sstream>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ResultCache.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/domain.hh"
#include "Driver/baseSolveStep.hh"
#include "Driver/basedriver.hh"
#include "Driver/singleDriver.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "OLAS/algsys/basesystem.hh"
#include "Optimization/Optimization.hh"
#include "PDE/StdPDE.hh"
#include "PDE/timestepping.hh"
#include "Utils/biotSavart.hh" 
#include "assemble.hh"
#include "boost/algorithm/string.hpp"
#include "stdSolveStep.hh"

namespace CoupledField {

class Hysteresis;

  StdSolveStep::StdSolveStep(StdPDE & apde)
    :BaseSolveStep(),
     PDE_(apde)
  {



    pdename_      = PDE_.GetName();
    numPDENodes_  = PDE_.getPDE_numPDENodes();
    numPDEElems_  = PDE_.getPDE_numElems();
    isaxi_        = PDE_.GetIsaxi();
    subdoms_      = PDE_.getPDE_subdoms();
    materialData_ = PDE_.getPDEMaterialData();
    ptgrid_       = PDE_.getPDE_grid();
    algsys_       = PDE_.getPDE_algsys();
    sol_          = PDE_.getPDESolution();
    assemble_     = PDE_.getPDE_assemble();
    TS_alg_       = PDE_.getTimeStepping();

    if( TS_alg_ != NULL ) {
      matrix_factor_ = TS_alg_->GetEffSysMatFactors();
    }

    eqnMap_       = PDE_.GetEqnMap();
    results_      = PDE_.GetResultInfos();
    numEqns_      = PDE_.GetSolutionVector()->GetSize();

    // nonlinear parameters
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;
    nonLin_           = PDE_.IsNonLin();
    isHyst_           = PDE_.IsHysteresis();
    totalFormulation_ = PDE_.IsTotaFormulation();
    regionNonLinTypes_ = PDE_.GetNonLinRegionTypes();

    // for direct coupled PDEs
    pdeId1_   = NO_PDE_ID;
    pdeId2_   = NO_PDE_ID;

    startStep_ = 1;
    etaMinCount_ = 0;

    // In the end, read nonlinear data from xml-file
    if( nonLin_ ) {
      ReadNonLinData();
    }

    mHandle_ = domain->GetMathParser()->GetNewHandle();
    mParser_ = domain->GetMathParser();
    mParser_->SetExpr(mHandle_,"step");
  }


  //! Destructor
  StdSolveStep::~StdSolveStep() {
    mParser_->ReleaseHandle(mHandle_);
  }


  // ======================================================
  // STATIC SOLVING SECTION
  // ======================================================

  void StdSolveStep::PreStepStatic( ) {

    // init RHS at this place, because e.g. forces of other PDEs are added
    // to RHS afterwards
    algsys_->InitRHS();

  }


  void StdSolveStep::PostStepStatic() {

    // increment coupling counter
    if ( PDE_.IsIterCoupled() ) {
      UInt& counter = PDE_.GetIterCoupledCounter();
      counter++;
    }
    
    // check for Biot Savart
    if ( PDE_.IsBiotSavart() ) {
      Vector<Double> & sol = 
          dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      Vector<Double>& magVecBiotSavart = 
          PDE_.GetBiotSavart()->CalcFieldAllEqns(false); 
      sol += magVecBiotSavart;
    }

  }


  void StdSolveStep::SolveStepStatic(PtrParamNode analysis_id, AdjointParameters* adjointParams) {

    if (nonLin_) {
      if ( totalFormulation_ )
        StepStaticNonLinTotal(analysis_id);
      else
        StepStaticNonLin(analysis_id);
    }
    else {
      StepStaticLin(analysis_id, adjointParams);
    }
  }


  void StdSolveStep::StepStaticLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) {

    assemble_->AssembleMatrices();

    // The RHS-sources and boundary conditions
    // have to be reassembled each time
    assemble_->AssembleLinRHS(adjointParams);
    PDE_.SetBCs();

    // store rhs vector back to PDE
    Vector<Double> tmpRhs;;
    algsys_->GetRHSVal(tmpRhs);
    PDE_.SaveRHS( tmpRhs.GetPointer(), tmpRhs.GetSize() );

    // Only if the matrices have changed (e.g. due to updated lagrangian
    // formulation) the system matrix has to be rebuild
    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }

    // Incorporate Boundary conitions and
    // recalc the preconditioner eventually
    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();

      algsys_->SetupSolver(analysis_id);
    }

    // Solve problem
    algsys_->Solve(analysis_id);

    // Get the solution and store it
    Vector<Double> tmpSol;
    algsys_->GetSolutionVal(tmpSol);
    PDE_.SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize());
  }


  void StdSolveStep::StepStaticNonLin(PtrParamNode analysis_id)
  {

    bool performOneMoreStep;
    Double incrementalErr, residualErr;

    Vector<Double> solInc( numEqns_ );

    //get actual solution
    Vector<Double>  actSol =
      dynamic_cast<Vector<Double>&>(*(PDE_.GetSolutionVector()));


    // set the boundary conditions
    PDE_.SetBCs();

//     //perform the load-steps
//     const UInt nrLoads = 5;
//     const Double loadFactors[nrLoads] = {0.25, 0.5, 0.75, 0.9, 1.0};

//     // currently just for testing!!
//     // loop over load factor
//     for ( UInt iload=0; iload<=nrLoads; iload++ ) {
//       info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactors[iload]);

//       // setup right hand side
//       Double RhsLinL2Norm = SetLinRHS(loadFactors[iload]);

    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // assemble nonlinear pars to RHS
      assemble_->AssembleNonLinRHS();

      // set iteration counter
      UInt iterationCounter=0;

      do
        {
          iterationCounter++;
          // RHS is already set up!!
      
          PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "nonLin", iterationCounter);
          
          // setup and solve new system (rhs is already set) =====================
          //assemble_.InitNonLinMatrices();
          assemble_->AssembleMatrices();

          algsys_->ConstructEffectiveMatrix(matrix_factor_);
          algsys_->BuildInDirichlet();
          algsys_->SetupPrecond();
          algsys_->SetupSolver(child_id);
          algsys_->Solve(child_id);

          // new solution is only an increment of the full solution =============
          algsys_->GetSolutionVal( solInc );

          Double residualL2Norm = 0.0;
          Double etaLineSearch  = 1.0;
          if ( lineSearch_ == "none" ) {
            actSol += solInc;
          }
          else {
            // true is for transient simulation
            residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
          }

          // store the new solution
          PDE_.SaveSolution( actSol.GetPointer(), actSol.GetSize() );

          if ( lineSearch_ == "none" ) {
            // recalculate RHS with new values to get new residual (f^(k+1))========
            algsys_->InitRHS(RhsLinVal_);
            assemble_->AssembleNonLinRHS();  

            // compute the norm of the residual
            Vector<Double> actRHS;
            algsys_->GetRHSVal(actRHS);

            // calculation of residual error =======================================
            residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )
          }

          // calculation of residual error =======================================
          if ( RhsLinL2Norm > 1.0 )
            residualErr = residualL2Norm / RhsLinL2Norm;
          else
            residualErr = residualL2Norm;

          // calculate incremental error ========================================
      	  Double solIncrL2Norm = solInc.NormL2();
          Double actSolL2Norm  = actSol.NormL2();

          if ( actSolL2Norm )
            incrementalErr = solIncrL2Norm / actSolL2Norm;
          else {
            incrementalErr = solIncrL2Norm;
            WARN("Zero solution vector!! ");
          }

          // output of norms and data
          if ( nonLinLogging_ == true )
            WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);

          // boolean variable, holds condition if another iteration step is necessary
          performOneMoreStep =
            (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

        } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } // load step loop

    if ( performOneMoreStep ) {
      Info->PrintF( pdename_, "NL scheme not converged\n" );
      std::stringstream out;
      out << "Residual error: " <<  residualErr << " ; Incremental error: " << incrementalErr 
          << std::endl << std::endl;
      std::string strAll = out.str();
      Info->PrintF( pdename_, strAll.c_str() );
    }
    
  }


  void StdSolveStep::StepStaticNonLinTotal(PtrParamNode analysis_id)
  {

    bool performOneMoreStep;
    std::cout << "In  stdSolveStep-NL-Total" << std::endl;

    Vector<Double> prevSol( numEqns_ );

    //get actual solution
    Vector<Double>  actSol =
      dynamic_cast<Vector<Double>&>(*(PDE_.GetSolutionVector()));

    // set the boundary conditions
    PDE_.SetBCs();

    //perform the load-steps
    Double loadFactor = 0.0;

    // currently just for testing!!
    // loop over load factor
    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      Info->PrintF(pdename_, "\n");
      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // set iteration counter
      UInt iterationCounter=0;

      do
        {
          iterationCounter++;

          //init RHS with linear part!
          algsys_->InitRHS(RhsLinVal_);

          // copy current solution, with which we start, to previous solution
          prevSol = actSol;

          // RHS is already set up!!
          PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "nonLin", iterationCounter);
          
          // setup and solve new system (rhs is already set) =====================
          //assemble_.InitNonLinMatrices();
          assemble_->AssembleMatrices();

          algsys_->ConstructEffectiveMatrix(matrix_factor_);
          algsys_->BuildInDirichlet();
          algsys_->SetupPrecond();
          algsys_->SetupSolver(child_id);
          algsys_->Solve(child_id);

          // new solution 
          algsys_->GetSolutionVal( actSol );
          //          std::cout << "actSol:\n" << actSol << std::endl;

          // store the new solution
          PDE_.SaveSolution( actSol.GetPointer(), actSol.GetSize() );

          // compute L2-Norm of error between last incremental solution and
          // actual incremental solution
          Double solIncrL2Norm=0;
          for (UInt i=0; i<actSol.GetSize(); i++)
            solIncrL2Norm += (actSol[i]-prevSol[i])*(actSol[i]-prevSol[i]);
    
          solIncrL2Norm = sqrt(solIncrL2Norm);
          Double actSolL2Norm = actSol.NormL2();
    
          Double incrementalErr, residualErr;
          if (actSolL2Norm > 1)
            incrementalErr = solIncrL2Norm / actSolL2Norm;
          else
            incrementalErr = solIncrL2Norm;
 
          //just dummies for log-file
          Double etaLineSearch = RhsLinL2Norm;
          etaLineSearch = 1.0;
          residualErr   = incrementalErr;

          // output of norms and data
          if ( nonLinLogging_ == true )
            WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);

          // boolean variable, holds condition if another iteration step is necessary
          performOneMoreStep =
            (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);

          //set previous solution of iteration to actual solution
          //prevSol = actSol;

        } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } // load step loop

  }


  // ======================================================
  // Solve Step Transient SECTION
  // ======================================================

  void StdSolveStep::PreStepTrans() {
    ResultCache::SetStepValue( actTime_ );

    // due to coupling-pdes, the RHS has to be initialized BEFORE
    // the coupling forces are assembled to the RHS
    algsys_->InitRHS();

    UInt step = (UInt) mParser_->Eval(mHandle_);

    PDE_.ReadDisplacementAndUpdateGrid( step );
  }


  void StdSolveStep::SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams) {


    if ( isHyst_ ) {
      StepTransNonLinHysteresis(analysis_id);
    }
    // do a nonlinear time step
    else if (nonLin_){
      if ( totalFormulation_ )
        StepTransNonLinTotal(analysis_id);
      else
        StepTransNonLin(analysis_id);
    }
    // do a linear time step
    else {
      StepTransLin(analysis_id, adjointParams);
    }
  }


  void StdSolveStep::StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) 
  {
    //account for RHS
    assemble_->AssembleLinRHS(adjointParams);
    PDE_.ComputeRHS( actTime_ );

    // store rhs vector back to PDE
    Vector<Double> tmpRhs;
    algsys_->GetRHSVal(tmpRhs);
    PDE_.SaveRHS( tmpRhs.GetPointer(), tmpRhs.GetSize() );

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
    bool isIterCoupled    = PDE_.IsIterCoupled();

    if(domain->GetOptimization() && domain->GetOptimization()->IsFirstTransientStepStatic()){
      if(actStep_ == 2 && adjointParams == NULL){ // reset everything but the solution after the first step
        PDE_.getTimeStepping()->ReInit();
      }
      if(actStep_ == 1 && adjointParams != NULL){ // reset everything before first step backwards
        ReInit(); // this sets everything = 0, so predictor does just add 0
      }
    }

    // perform predictor step: if we have an iterative coupled
    // PDE-system, we should perform the predictor state just
    // in the first iteration
    if ( isIterCoupled == false || iterCoupledCounter == 0 ) {
      Vector<Double> & solHelp =
        dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      TS_alg_->Predictor(solHelp);
    }

    assemble_->AssembleMatrices();
    bool effectiveMatrixUpdated = false;

    // calculate first step static, and the adjoint will also be static and restore the parameters after that
    if(domain->GetOptimization() && domain->GetOptimization()->IsFirstTransientStepStatic() && actStep_ == 1){
      double tMass = matrix_factor_[MASS];        
      double tDamp = matrix_factor_[DAMPING];
      matrix_factor_[MASS] = 0.0;
      matrix_factor_[DAMPING] = 0.0;
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      effectiveMatrixIsWrong_ = true;
      matrix_factor_[MASS] = tMass;
      matrix_factor_[DAMPING] = tDamp;
      effectiveMatrixUpdated = true;
    }else if(assemble_->IsMatrixUpdated() || effectiveMatrixIsWrong_ ){
      // we have to construct the effective matrix in 2 cases
      // either, the matrix was updated
      // or we have the second step after the first step had been static (the time runs forward, not in adjoint case)
      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      effectiveMatrixIsWrong_ = false;
      effectiveMatrixUpdated = true;
    }

    // update Right Hand Side
    TS_alg_->UpdateRHS();

    // set boundary conditions
    PDE_.SetBCs();
    algsys_->BuildInDirichlet();

    if( effectiveMatrixUpdated ){
      algsys_->SetupPrecond( );
      algsys_->SetupSolver(analysis_id);
    }

    algsys_->Solve(analysis_id);

    Vector<Double> tmpSol;
    algsys_->GetSolutionVal( tmpSol );
    PDE_.SaveSolution(tmpSol.GetPointer(),tmpSol.GetSize());

    Vector<Double> & solHelp =
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());


    TS_alg_->Corrector(solHelp);

    if ( isIterCoupled ) {
      iterCoupledCounter++;
    }

  }



  void StdSolveStep::StepTransNonLin(PtrParamNode base_analysis_id) {


    bool performOneMoreStep;
    Vector<Double> solInc( numEqns_ );

    Double residualErr, incrementalErr;
    //get actual solution
    Vector<Double>  actSol =
      dynamic_cast<Vector<Double>&>(*(PDE_.GetSolutionVector()));

    // perform predictor step
    if ( TS_alg_== NULL ) {
      EXCEPTION ( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin" );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs();

    // currently just for testing!!
    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS();  

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS(actSol);

      // set iteration counter
      UInt iterationCounter=0;


      do {
        iterationCounter++;

        // RHS is already set up!!
        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(base_analysis_id, "load", iload, "nonLin", iterationCounter);

        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver(child_id);
        algsys_->Solve(child_id);

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solInc );

				// initialize residualL2Norm or receive compiler warning
        Double residualL2Norm = 0.0;
        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none" ) {
          actSol += solInc;
        }
        else {
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch, true);
        }

        //store A_(n+1) in the solution-object sol_
        PDE_.SaveSolution( actSol.GetPointer(), actSol.GetSize() );

        if ( lineSearch_ == "none" ) {
          // calculation of error norms
          // recalculate RHS with new values to get new residual (f^(k+1))=======
          algsys_->InitRHS(RhsLinVal_);

          //Update RHS (mass matrix on right hand side)
          TS_alg_->UpdateRHS(actSol);

          // inner forces due to nonlin formulation
          assemble_->AssembleNonLinRHS();

          Vector<Double> actRHS;
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
        }

        if ( RhsLinL2Norm > 1.0 )
          residualErr    = residualL2Norm /  RhsLinL2Norm;
        else
          residualErr    = residualL2Norm;

        // calculate incremental error
        Double solIncrL2Norm = solInc.NormL2();
        Double actSolL2Norm = actSol.NormL2();
 
        if ( actSolL2Norm > 1.0)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        // --------------------------------------------------------------------
        // output of norms and data
        // --------------------------------------------------------------------
        if ( nonLinLogging_ == true )
          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);

        // boolean variable, holds condition if another iteration step
        // is necessary
        performOneMoreStep =
          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } // load step loop

    if ( performOneMoreStep ) {
      UInt actStep = (UInt) mParser_->Eval(mHandle_);
      std::string strAll = "NL scheme for time step ";
      std::string nrStep;
      std::stringstream out;
      out << actStep;
      strAll +=  out.str() + " not converged!\n";
      Info->PrintF( pdename_, strAll.c_str() );
      out.str(std::string());
      out << "Residual error: " <<  residualErr << " ; Incremental error: " << incrementalErr 
          << std::endl << std::endl;
      strAll = out.str();
      Info->PrintF( pdename_, strAll.c_str() );
    }

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();
    if ( PDE_.IsIterCoupled() )
    {
      ++iterCoupledCounter;
    }

    //perform corrector step
    TS_alg_->Corrector(actSol);
  }


  void StdSolveStep::StepTransNonLinTotal(PtrParamNode base_analysis_id) {

    std::cout << "\n In :StepTransNonLinTotal  \n " << std::endl;
 
    bool performOneMoreStep;
    UInt iterationCounter=0;
    Double etaLineSearch;
    Double incrementalErr, residualErr;

    Vector<Double> newSol( numEqns_ ); 
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> solPrev( numEqns_ );
    Vector<Double> solMinus( numEqns_ );

    oldSol.Init(0);
    newSol.Init(0);
  
    //get solution of previous time step  
    Vector<Double> & solHelp = 
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

    solPrev = solHelp;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      Exception( "TS_alg has NULL-Pointer, in StdSolveStep::StepTransNonLin");
    }
    else {
      //compute predictors
      TS_alg_->Predictor(solPrev);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs( );    

    Double loadFactor = 0.0;
    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      Info->PrintF(pdename_, "\n");
      Info->PrintF(pdename_, "LoadFactor: %g \n", loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

 
      do {
        iterationCounter++;
        
        // set solution of previous iteration
        if (iterationCounter == 1 ) {
          oldSol = solPrev;
        }
        else {
          oldSol = newSol;
        }
        
        //init RHS with linear part!
        algsys_->InitRHS(RhsLinVal_);
 
        // inner forces due to nonlin formulation
        // currently just fixed point, so AssembleNonLinRHS() is not necessary
        // assemble_->AssembleNonLinRHS();  
       
        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(base_analysis_id, "load", iload, "nonLin", iterationCounter);
        
        //perform new assembly
        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        
        
        //Update RHS (mass matrix on right hand side)
        TS_alg_->UpdateRHS();
        
        // build in the Dirichlet vales in system matrix and rhs
        algsys_->BuildInDirichlet();
        
        algsys_->SetupSolver(child_id);
        algsys_->SetupPrecond();
        
        algsys_->Solve(child_id);
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
        
        if (actSolL2Norm > 1)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        //just dummy things
        etaLineSearch = RhsLinL2Norm;
        etaLineSearch = 1.0;
        residualErr   = incrementalErr;
        
//         //computation of residuum ====== Start
//         assemble_->AssembleMatrices();
//         algsys_->InitRHS(RhsLinVal_);

//         //RHS - timeparam*M ( solpred - actsol )
//         TS_alg_->UpdateRHS(newSol);
//         solMinus = -newSol;
//         algsys_->UpdateRHS(STIFFNESS,solMinus);

//         Vector<Double> actRHS;
//         algsys_->GetRHSVal( actRHS );
//         residualErr = actRHS.NormL2();
//         std::cout << "RhsLinL2Norm: " << RhsLinL2Norm << std::endl;

//         if (  RhsLinL2Norm > 1.0 )
//           residualErr /= RhsLinL2Norm;
        
//         //computation of residuum ====== END

        // output of norms and data
        nonLinLogging_ = true;
        if ( nonLinLogging_ == true )
          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);
        
        performOneMoreStep = (incrementalErr > incStopCrit_) ; //|| (residualErr > residualStopCrit_);      
        
      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
    } // load step loop

    if ( incrementalErr > incStopCrit_ ) 
      std::cout << "Not converged, norm is: " <<incrementalErr << std::endl; 

    //perform corrector step  
    TS_alg_->Corrector(newSol);
  }

  
  void StdSolveStep::StepTransNonLinHysteresis(PtrParamNode analysis_id) {
    bool performOneMoreStep;

    Vector<Double> solInc( numEqns_ );
    Vector<Double> actSol( numEqns_ );
    Vector<Double> oldSol( numEqns_ );
    Vector<Double> prevSol( numEqns_ );
    actSol.Init();

    // get actual solution
    algsys_->GetSolutionVal( actSol );

    //solution from previous time step
    prevSol = actSol;

    // perform predictor step
    if ( TS_alg_== NULL ) {
      EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin" );
    }
    else {
      //compute predictors
      TS_alg_->Predictor(actSol);
    }

    //! account for Dirichlet BCs
    PDE_.SetBCs();

    // currently just for testing!!
    // loop over load factor
    Double loadFactor = 0.0;

    for ( UInt iload=0; iload<1; iload++ ) {
      loadFactor += 1.0;
      info->Get("PDE")->Get(pdename_)->Get("load_factor")->SetValue(loadFactor);

      // setup right hand side
      Double RhsLinL2Norm = SetLinRHS(loadFactor);

      // inner forces due to nonlin formulation
      assemble_->AssembleNonLinRHS();  

      //Update RHS (mass matrix on right hand side)
      TS_alg_->UpdateRHS();

      // set iteration counter
      UInt iterationCounter=0;

      do {

        iterationCounter++;
        oldSol = actSol;

        PtrParamNode child_id = BaseDriver::CreateAnalysisIdChild(analysis_id, "load", iload, "nonLin", iterationCounter);

        //        RHS is already set up!!
        if ( iterationCounter > 0 ) {
          // setup linear part of right hand side
          algsys_->InitRHS(RhsLinVal_);

          // inner forces due to nonlin formulation
          assemble_->AssembleNonLinRHS();  

          //Update RHS (mass matrix on right hand side)
          TS_alg_->UpdateRHS();
        }

        assemble_->AssembleMatrices();
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();

        algsys_->SetupPrecond();
        algsys_->SetupSolver(child_id);
        algsys_->Solve(child_id);

        // new solution is only an increment of the full solution =============
        algsys_->GetSolutionVal( solInc );

        Double residualL2Norm;
        Double etaLineSearch = 1.0;

        if ( lineSearch_ == "none" ) {
          actSol = solInc;
        }
        else {
          EXCEPTION("Currently lineSreach not supported" );
          residualL2Norm = LineSearch(solInc, actSol, etaLineSearch);
        }

        //store actual solution to the solution-object sol_
        PDE_.SaveSolution(actSol.GetPointer(), actSol.GetSize() );

        if ( lineSearch_ == "none" ) {
          // calculation of error norms
          // recalculate RHS with new values to get new residual (f^(k+1))=======
          algsys_->InitRHS(RhsLinVal_);

          //substract from RHS: intFactro*MASS*acc - intFactor*DAMP*vel
          TS_alg_->UpdateRHS(actSol);

          // substracte stiff-matrix from RHS
          TS_alg_->SubstractStiffnessFromRHS(actSol);

          // inner forces due to nonlin formulation
          //assemble_->AssembleNonLinRHS( actTime_ );

          Vector<Double> actRHS;
          algsys_->GetRHSVal( actRHS );

          // calculation of residual error: L2Norm of ( f_i^(k+1) - f_a )
          residualL2Norm = PDE_.GetRhsL2Norm(actRHS);
        }

        Double residualErr;
        if ( RhsLinL2Norm > 1.0 )
          residualErr    = residualL2Norm /  RhsLinL2Norm;
        else
          residualErr    = residualL2Norm;


        residualErr = 0.0;

        // calculate incremental error

        // compute L2-Norm of error between last incremental solution and
        // actual incremental solution
        Double solIncrL2Norm=0;
        for (UInt i=0; i<actSol.GetSize(); i++)
          solIncrL2Norm += (actSol[i]-oldSol[i])*(actSol[i]-oldSol[i]);

        solIncrL2Norm = sqrt(solIncrL2Norm);
        Double actSolL2Norm = actSol.NormL2();

        Double incrementalErr;
        if (actSolL2Norm > 1)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          incrementalErr = solIncrL2Norm;

        // --------------------------------------------------------------------
        // output of norms and data
        // --------------------------------------------------------------------
        if ( nonLinLogging_ == true )
          WriteNonLinIterToInfoXML(pdename_, iterationCounter, residualErr, incrementalErr, etaLineSearch);

        // boolean variable, holds condition if another iteration step
        // is necessary
        performOneMoreStep =
          (incrementalErr > incStopCrit_)||(residualErr > residualStopCrit_);

      } while(performOneMoreStep && iterationCounter < nonLinMaxIter_);

    } // load step loop

    //perform corrector step
    TS_alg_->Corrector(actSol);
  }


  void StdSolveStep::PostStepTrans( ) {


    Vector<Double> & solHelp =
      dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());

    // Following method is essential for fractional damping model
    TS_alg_->AdvanceTimestep(solHelp);
    
    // check for Biot Savart
    if ( PDE_.IsBiotSavart() ) {
      Vector<Double> & sol = 
          dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector());
      Vector<Double>& magVecBiotSavart = 
          PDE_.GetBiotSavart()->CalcFieldAllEqns(false);
      sol += magVecBiotSavart;
    }
  }



  // ======================================================
  // Solve Step Harmonic  SECTION
  // ======================================================w

  void StdSolveStep::PreStepHarmonic() {
    ResultCache::SetStepValue( actFreq_ );
    algsys_->InitRHS();
  }


  void StdSolveStep::SolveStepHarmonic(PtrParamNode analysis_id, AdjointParameters* adjointParams) {
    if ( nonLin_ ) {
      StepHarmonicNonLin(analysis_id);
    }
    else {
      StepHarmonicLin(analysis_id, adjointParams);
    }
  }


  void StdSolveStep::StepHarmonicLin(PtrParamNode analysis_id, AdjointParameters* adjointParams) {


    //this has to be done each frequency!
    assemble_->AssembleLinRHS(adjointParams);

    if ( PDE_.IsComputeRHS4HarmSet() ) {
      // Evaluating RHS with nodal srcs for harmonic flownoise problems
      PDE_.ComputeRHS(actFreq_);
    }

    assemble_->AssembleMatrices( );

    PDE_.SetBCs();

    // store rhs vector back to PDE
    Vector<Complex> tmpSol;
    algsys_->GetRHSVal( tmpSol );
    PDE_.SaveRHS( tmpSol.GetPointer(), tmpSol.GetSize() );

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->ConstructEffectiveMatrix( matrix_factor_ );
    }

    algsys_->BuildInDirichlet();

    if( assemble_->IsMatrixUpdated() ) {
      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
    }

    algsys_->Solve(analysis_id);

    algsys_->GetSolutionVal(tmpSol);
    PDE_.SaveSolution(tmpSol.GetPointer(), tmpSol.GetSize());
  }


  // ======================================================
  // METHODS FOR EIGENVALUE COMPUTATION
  // ======================================================

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Double> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, false );

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds);

    return frequencies.GetSize();
  }

  UInt StdSolveStep::CalcEigenFrequencies( Vector<Complex> & frequencies,
                                           Vector<Double> & errBounds,
                                           UInt numFreq, Double shift ) {

    // Init algsys data structures
    algsys_->InitRHS();
    algsys_->InitSol();
    algsys_->InitMatrix();

    assemble_->AssembleMatrices();

    // Setup solver
    algsys_->SetupEigenSolver( numFreq, shift, true );

    // Calculate eigenfrequencies
    algsys_->CalcEigenFrequencies( frequencies, errBounds );

    return frequencies.GetSize();
  }

  void StdSolveStep::CalcEigenMode( UInt numMode ) {


    algsys_->CalcEigenMode( numMode );

    // Get the solution and store it
    Vector<Double> tmpSol;
    algsys_->GetSolutionVal(tmpSol);
    PDE_.SaveSolution( tmpSol.GetPointer(), tmpSol.GetSize());
  }




  // ======================================================
  // METHODS FOR NONLINEAR ANALYSIS
  // ======================================================

  // sets excitation coil and returns L2Norm of them
  Double StdSolveStep::SetLinRHS( Double loadFactor)
  {

    Double RhsLinL2Norm;


    // to incorporate loads
    assemble_->AssembleLinRHS(); 

    // Stores rhs vector into extForces and returns that L2-norm
    algsys_->GetRHSVal( RhsLinVal_ );
    RhsLinVal_ *= loadFactor;

    RhsLinL2Norm = RhsLinVal_.NormL2();

    // If extForcesL2Norm is 0, no residual norm can be calculated
    if (!RhsLinL2Norm) {
      WARN("Zero external force vector!! ");
    }

    return RhsLinL2Norm;
  }

  UInt StdSolveStep::SetDeltaLinRHS()
  {

    // to incorporate loads
    assemble_->AssembleLinRHS(); 
    
    Vector<Double> newRhsLinVal; //!< external forces (for nonlin simulations)
    algsys_->GetRHSVal( newRhsLinVal );
    DeltaRhsLinVal_=newRhsLinVal-tmpOldRhsLinVal_;
    
    RhsLinVal_=tmpOldRhsLinVal_;

    Double DeltaNorm=DeltaRhsLinVal_.NormL2();
    Double oldNorm  =tmpOldRhsLinVal_.NormL2();
    Double aux;
    UInt nrLoadSteps;

    UInt minNrLoadSteps = 1;
    UInt maxNrLoadSteps = 1;

    if(oldNorm<1e-13)
      nrLoadSteps=maxNrLoadSteps;
    else{
      aux= Double(maxNrLoadSteps)*(DeltaNorm/oldNorm);
      nrLoadSteps=UInt(aux);
    }
    if(nrLoadSteps<minNrLoadSteps)
      nrLoadSteps=minNrLoadSteps;
    else if (nrLoadSteps>maxNrLoadSteps)
      nrLoadSteps=maxNrLoadSteps;

    oldRhsLinVal_=tmpOldRhsLinVal_;
    tmpOldRhsLinVal_=newRhsLinVal;

    return nrLoadSteps;
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void StdSolveStep::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {


    //const UInt numElems = numPDENodes_ * dofspernode_;
    vec.Resize( numEqns_ );
    vec.Init();

    for (UInt i=0; i<numEqns_; i++) {
      vec[i] = pt[i];
    }
  }



  Double StdSolveStep::LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol,
                                  Double& etaLineSearch, bool trans)
  {

    Vector<Double> solOld(actSol);
    const UInt nrEtas = 5;
    const Double eta[nrEtas] = {1, 0.5, 0.25, 0.125, 0.01};
		// initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    for( UInt i=0; i<nrEtas; i++) {
      actSol = solIncrement * eta[i];
      actSol += solOld;

      //store new solution
      PDE_.SaveSolution(actSol.GetPointer(),actSol.GetSize());

      // recalculate RHS with new values to get new residual (f^(k+1))========
      algsys_->InitRHS(RhsLinVal_ );

      if( trans ) {
        assemble_->AssembleNonLinRHS();
        TS_alg_->UpdateRHS(actSol);
      }
      else {
        assemble_->AssembleNonLinRHS();
      }


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  ( f_i^(k+1) - f_a )

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

//     if ( etaOpt == eta[nrEtas-1] ) 
//       etaMinCount_++;
//     else 
//       etaMinCount_ = 0;
//     if ( etaMinCount_ == 10 ) {
//       //now try to boost the current solution out of a local minima
//       etaOpt = eta[0];
//       etaMinCount_ = 0;
//     }

    etaLineSearch = etaOpt;

    actSol  = solIncrement * etaOpt;
    actSol += solOld;

    return residualL2NormOpt;
  }


  Double StdSolveStep::LineSearchMaterial(Vector<Double>& solIncrement, Vector<Double>& actSol,
                                          Double& etaLineSearch, Double& RhsLinL2Norm, bool trans)
  {

    Vector<Double> solOld(actSol);
    const UInt nrEtas = 3;
    const Double eta[nrEtas] = {0.9, 0.5, 0.3};
    //    const Double eta[nrEtas] = {0.1, 0.2, 0.4, 0.5, 0.7, 0.9, 1.0};
		// initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;

    Vector<Double> tmpSol;
    algsys_->GetSolutionVal( tmpSol );

    for( UInt i=0; i<nrEtas; i++) {
      //       if (i>0)
      //         actSol=solOld;
      actSol = solOld + solIncrement * eta[i];
      //      actSol -= solOld;


      //store new solution
      tmpSol = actSol;
      PDE_.SaveSolution( tmpSol.GetPointer(), tmpSol.GetSize() );

      // Recalculate residual, f-Cu-Mu-K*u
      algsys_->InitRHS();

      assemble_->AssembleLinRHS();

      // assemble!
      assemble_->AssembleMatrices();

      // account for Dirichlet BCs
      PDE_.SetBCs();

      algsys_->ConstructEffectiveMatrix(matrix_factor_);

      algsys_->BuildInDirichlet();

      TS_alg_->UpdateRHS(actSol);
      // substract K^* u^k from RHS

      algsys_->RemoveIDBCInfoFromMatrix();

      TS_alg_->SubstractStiffnessFromRHS(actSol);


      // =====================================================================
      // calculation of error norms
      // =====================================================================
      Vector<Double> actRHS;
      algsys_->GetRHSVal( actRHS );

      // calculation of residual error =======================================
      Double residualL2Norm = PDE_.GetRhsL2Norm(actRHS); // L2Norm of  (f-Ku )

      // Double residualErr; // TODO: Unused variable residualErr

      // if ( RhsLinL2Norm > 1.0 )
      //   residualErr    = residualL2Norm /  RhsLinL2Norm;
      // else
      //  residualErr    = residualL2Norm;

      if (residualL2Norm < residualL2NormOpt) {
        residualL2NormOpt = residualL2Norm;
        etaOpt = eta[i];
      }
    }

    etaLineSearch = etaOpt;

    actSol  = solOld + solIncrement * etaOpt;

    return residualL2NormOpt;
  }


  // returns that L2-norm of an algsys vector
  Double StdSolveStep::AlgsysL2Norm(Double * pt)
  {

    //const UInt numElems = numPDENodes_ * dofspernode_;
    Double quadSum = 0;

    for (UInt i=0; i<numEqns_; i++)
      quadSum += pt[i]*pt[i];

    return sqrt(quadSum);
  }

  // read nonlinear parameters from xml file
  void StdSolveStep::ReadNonLinData() {

    // Get ParamNode of pde
    PtrParamNode nonLinNode;
    std::string pdeName = PDE_.GetName();

    // SPECIAL: If PDE-name contains the "-" sign,
    // we assume that we have a coupledPDE
    if(boost::find_first(pdeName, "-")) {
      // go to direct-coupling List
      // iterate over all nodes and pick the one
      // which has a "nonLinear" element
      // get current multiSequenceStep
      UInt actMsStep =
        domain->GetSingleDriver()->GetActSequenceStep();
      StdVector<PtrParamNode> pairCouplings;
      PtrParamNode cplList =
        param->GetByVal("sequenceStep",std::string("index"), actMsStep)
        ->Get("couplingList", ParamNode::PASS);
      if( cplList ) {
        StdVector<PtrParamNode> pairCouplings =
          cplList->Get("direct")->GetChildren();
        // look for each direct coupling, if it has a "nonLin"-node
        for( UInt iCpl = 0; iCpl < pairCouplings.GetSize(); iCpl++ ) {
          nonLinNode = pairCouplings[iCpl]->Get("nonLinear", ParamNode::PASS );
          if( nonLinNode ) {
            break;
          }
        }
      }
      if( !nonLinNode ) {
        // in this case we iterate over all single PDEs and try to
        // find the related entry
        ParamNodeList pdes =
          param->GetByVal("sequenceStep",std::string("index"), actMsStep)
          ->Get("pdeList")->GetChildren();
        for (UInt iPde = 0; iPde < pdes.GetSize(); iPde++ ) {
          if( boost::find_first(pdeName, pdes[iPde]->GetName() ) ) {
            nonLinNode = pdes[iPde]->Get("nonLinear", ParamNode::PASS );
            if( nonLinNode ) {
              break;
            }
          }
        }
      }
    } else {
      nonLinNode = PDE_.GetParamNode()->Get("nonLinear", ParamNode::PASS );
    }

    // Check, if any nonlinear node was found
    if( !nonLinNode ) {
      WARN("Taking default parameters for nonlinear data" );
    }

    // Read data, if "nonLinear" element was found
    if( nonLinNode ) {

      // solution method
      nonLinNode->GetValue( "method", nonLinMethod_, ParamNode::PASS );

      // perform logging?
      nonLinNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );

      // type of line search
      nonLinNode->Get( "lineSearch")->GetValue( "type", lineSearch_ );

      // incremental stopping criterion
      nonLinNode->GetValue( "incStopCrit", incStopCrit_ );

      // residual stopping criterion
      nonLinNode->GetValue( "resStopCrit", residualStopCrit_ );

      // maximal number of NL-iterations
      nonLinNode->GetValue( "maxNumIters", nonLinMaxIter_ );
    }

  }

  void StdSolveStep::ReInit(){
    dynamic_cast<Vector<Double>&>(*PDE_.GetSolutionVector()).Init();
    PDE_.getTimeStepping()->ReInit();
  }
  
  void StdSolveStep::WriteNonLinIterToInfoXML(const std::string& pdeName, const UInt iterationCounter,
          const Double residualErr, const Double incrementalErr, double etaLineSearch)
  {
    PtrParamNode iter = info->Get("PDE")->Get(pdeName)->Get("nonlinear_iteration", ParamNode::APPEND); 
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("residualErr")->SetValue(residualErr);
    iter->Get("incrementalErr")->SetValue(incrementalErr);
    if(etaLineSearch)
      iter->Get("eta_linesearch")->SetValue(etaLineSearch);
    // update convergen process
    info->ToFile();
  }

} // end of namespace

