// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "solveStepFluidMech.hh"

#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/basesystem.hh"

namespace CoupledField {

  // declare logging stream
  DECLARE_LOG(solvestepfluidmech)
  DEFINE_LOG(solvestepfluidmech, "solvestepfluidmech")

    SolveStepFluidMech::SolveStepFluidMech(StdPDE& apde) : StdSolveStep(apde) 
  {

    isInstationary_ = PDE_.IsInstationary();
    //isTrapezoidal_  = PDE_.IsTrapezoidal();
    //isNewton_       = PDE_.IsNewton();

  }

  SolveStepFluidMech::~SolveStepFluidMech() {

  }

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepFluidMech::StepTransNonLin(InfoNode* analysis_base)
  {

    assemble_->AssembleLinRHS();

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();

    //std::cout << "iterCoupledCounter:" << iterCoupledCounter << std::endl;

    bool isIterCoupled    = PDE_.IsIterCoupled();
    Double incrementalErr, incrementalVeloErr, incrementalPresErr;

    std::string pdeNameLong(pdename_);
    
    pdeNameLong += "-PDE: ";

    bool performOneMoreStep;
    // TODO: Check if this is still needed
    // Integer eqnNr;
    // UInt  eqnDof;
  
    Vector<Double>  actSol;
    Vector<Double>  actVelo, actPres, tmpVelo, tmpPres;

    sol_->GetAlgSysVector(actSol);
    sol_->GetGlobalSolVector(FLUIDMECH_VELOCITY,actVelo);
    sol_->GetGlobalSolVector(FLUIDMECH_PRESSURE,actPres);

    LOG_DBG2(solvestepfluidmech) << "actVelo\n" << actVelo << std::endl;

    LOG_DBG2(solvestepfluidmech) << "actPres\n" << actPres << std::endl;

    Vector<Double> newSol;
    Vector<Double> solIncrement, solVelocityInc, solPressureInc;
    newSol.Resize( numEqns_ );
    solIncrement.Resize( numEqns_ );

    // perform predictor step
    if ( isInstationary_ ) {
      if ( TS_alg_== NULL ) {
        EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin");
      }
      else {
        if ( isIterCoupled == false || iterCoupledCounter == 0 ) {        
          TS_alg_->Predictor(actSol);
        }
      }
    }
    
    PDE_.SetBCs();

    // inner forces due to nonlin formulation
    //assemble_->AssembleNonLinRHS( actTime_ );  

    UInt iterationCounter=0;
    Double iterTime=0;

    PDE_.InitStabParams();
    //Update RHS (mass matrix on right hand side)
    if ( isInstationary_ ) {
      assemble_->AssembleMatrices();
      TS_alg_->UpdateRHS();
    }

    do
      {
        iterationCounter++;
        iterTime+=0.01;

        LOG_DBG(solvestepfluidmech) << "loop=" << iterationCounter << " newton=" << isNewton_;

        // setup and solve new system (rhs is already set) =====================
        SETPROFILE("Before AssembleMatrices");
        if ( !isInstationary_ || iterationCounter != 1){
          assemble_->AssembleMatrices();
          PDE_.PrintStabParams();
        }
        
        InfoNode* analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);
        
        SETPROFILE("After AssembleMatrices");
        algsys_->ConstructEffectiveMatrix(matrix_factor_);
        algsys_->BuildInDirichlet();

        SETPROFILE("Before Solve");
        algsys_->SetupPrecond();
        algsys_->SetupSolver(analysis_id);
        algsys_->Solve(analysis_id);   
        SETPROFILE("After Solve");

        // new solution is NOT only an increment of the full solution =============
        algsys_->GetSolutionVal( newSol );
        LOG_DBG(solvestepfluidmech) << "newSol\n" << newSol << std::endl;
        
        sol_->SetAlgSysVector(newSol);

        sol_->GetGlobalSolVector(FLUIDMECH_VELOCITY,tmpVelo);
        sol_->GetGlobalSolVector(FLUIDMECH_PRESSURE,tmpPres);

        LOG_DBG2(solvestepfluidmech) << "newVelo\n" << tmpVelo << std::endl;
        LOG_DBG2(solvestepfluidmech) << "newPres\n" << tmpPres << std::endl;

        solVelocityInc=tmpVelo-actVelo;
        solPressureInc=tmpPres-actPres;
        solIncrement=newSol-actSol;

        actVelo=tmpVelo;
        actPres=tmpPres;
        actSol=newSol;

        //PDE_.WriteResultsInFile(iterationCounter, actTime_+iterTime);

        // calculate incremental error ========================================
        Double solVelocityIncL2Norm = solVelocityInc.NormL2();
        Double solPressureIncL2Norm = solPressureInc.NormL2();
        Double solIncrL2Norm = solIncrement.NormL2();

        Double actVeloL2Norm  = tmpVelo.NormL2();
        Double actPresL2Norm  = tmpPres.NormL2();
        Double actSolL2Norm  = actSol.NormL2();

        // TODO: Check if this is still needed
        // Double etaLineSearch = 0;
        // Double residualErr = 0;

        if (actVeloL2Norm > 1.0){
          incrementalVeloErr = solVelocityIncL2Norm / actVeloL2Norm;
        }
        else {
          incrementalVeloErr = solVelocityIncL2Norm;
          //Warning("Zero velocity solution vector!! ", __FILE__,__LINE__);      
        }
          
        if (actPresL2Norm > 1.0){
          incrementalPresErr = solPressureIncL2Norm / actPresL2Norm;
        }
        else {
          incrementalPresErr = solPressureIncL2Norm;
          //Warning("Zero pressure solution vector!! ", __FILE__,__LINE__);      
        }

        if (actSolL2Norm > 1.0)
          incrementalErr = solIncrL2Norm / actSolL2Norm;
        else
          {
            incrementalErr = solIncrL2Norm;
            //Warning("Zero solution vector!! ", __FILE__,__LINE__);      
          }
        
        if ( nonLinLogging_ == true ) {

          *(Info->GetInfoStreamPointer()) << std::endl << pdeNameLong << "NONLINEAR ITERATION "
                                          << iterationCounter 
                                          << " ==========================================\n"
                                          << pdeNameLong << "=== Norm of Solution     " << actSolL2Norm
                                          << std::endl
                                          << pdeNameLong << "=== Norm of Velo Sol     " << actVeloL2Norm
                                          << std::endl
                                          << pdeNameLong << "=== Norm of Pres Sol     " << actPresL2Norm
                                          << std::endl
                                          << pdeNameLong << "=== Incr. Sol error      " << incrementalErr 
                                          << std::endl
                                          << pdeNameLong << "=== Incr. Vel error      " << incrementalVeloErr 
                                          << std::endl
                                          << pdeNameLong << "=== Incr. Pres error     " << incrementalPresErr 
                                          << std::endl;
        }
        
        // boolean variable, holds condition if another iteration step is necessary
        performOneMoreStep = (incrementalErr > incStopCrit_ || incrementalVeloErr > incStopCrit_ || incrementalPresErr > incStopCrit_);
        //performOneMoreStep = (incrementalErr > incStopCrit_ || incrementalVeloErr > incStopCrit_);      

        //           mycout << "Solu = " << actSolL2Norm  
        //                  << "; Velo = " << actSolL2Norm  
        //                  << "; Pres = " << actSolL2Norm  
        //                  << "\tinc Err = " << incrementalErr
        //                  << "; inc Velo Err = " << incrementalVeloErr
        //                  << "; inc Pres Err = " << incrementalPresErr
        //                  << myendl;

      }while(performOneMoreStep && iterationCounter < nonLinMaxIter_);  
      
    if (incrementalErr > 50*incStopCrit_ || incrementalVeloErr > 50*incStopCrit_ || incrementalPresErr > 50*incStopCrit_)
      Warning("FluidMech did not converged!!! Simulation SHOULD be aborted",__FILE__,__LINE__);
    else if (incrementalErr > incStopCrit_ || incrementalVeloErr > incStopCrit_ || incrementalPresErr > incStopCrit_){
      Warning("FluidMech did not converged!!!",__FILE__,__LINE__);
    }
           
    if ( isInstationary_ ) {
      TS_alg_->Corrector(actSol);
    }

    LOG_DBG(solvestepfluidmech) << "actSol\n" << actSol << std::endl;

    if ( isIterCoupled ) {
      iterCoupledCounter++;
    }
    PDE_.AcouSourceCalc();
  }
} // end of namespace

