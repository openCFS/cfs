#include <fstream>
#include <iostream>
#include <string>

#include "solveStepFluidMech.hh"

#include "assemble.hh"

#include "Utils/nodestoresol.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "OLAS/algsys/basesystem.hh"

namespace CoupledField
{
  // declare logging stream
  DECLARE_LOG(solvestepfluidmech)
  DEFINE_LOG(solvestepfluidmech, "solvestepfluidmech")

  SolveStepFluidMech::SolveStepFluidMech(StdPDE& apde) : StdSolveStep(apde)
  {
    isInstationary_ = PDE_.IsInstationary();
  }

  SolveStepFluidMech::~SolveStepFluidMech()
  {
  }

  // ======================================================
  // Solve Step Static SECTION
  // ======================================================

  void SolveStepFluidMech::StepTransNonLin(PtrParamNode analysis_base)
  {
    assemble_->AssembleLinRHS();

    UInt& iterCoupledCounter = PDE_.GetIterCoupledCounter();

    bool isIterCoupled    = PDE_.IsIterCoupled();
    Double incrementalErr = 0.0;
    Double incrementalVeloErr = 0.0;
    Double incrementalPresErr = 0.0;

    std::string pdeNameLong(pdename_);

    pdeNameLong += "-PDE: ";

    bool performOneMoreStep = true;

    sol_->GetAlgSysVector(actSol_);
    sol_->GetGlobalSolVector(FLUIDMECH_VELOCITY,actVelo_);
    sol_->GetGlobalSolVector(FLUIDMECH_PRESSURE,actPres_);

    LOG_DBG2(solvestepfluidmech) << "actVelo_\n" << actVelo_ << std::endl;

    LOG_DBG2(solvestepfluidmech) << "actPres_\n" << actPres_ << std::endl;

    newSol_.Resize( numEqns_ );
    solIncrement_.Resize( numEqns_ );

    // perform predictor step
    if ( isInstationary_ )
    {
      if ( TS_alg_== NULL )
      {
        EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin");
      } else {
        if ( isIterCoupled == false || iterCoupledCounter == 0 )
        {
          TS_alg_->Predictor(actSol_);
        }
      }
    }

    UInt iterationCounter=0;

    PDE_.InitStabParams();
    //Update RHS (mass matrix on right hand side)
    if ( isInstationary_ )
    {
      assemble_->AssembleMatrices();
      TS_alg_->UpdateRHS();
    }

    PDE_.SetBCs();

    while (performOneMoreStep && iterationCounter < nonLinMaxIter_)
    {
      iterationCounter++;

      // setup and solve new system (rhs is already set) =====================
      if ( !isInstationary_ || iterationCounter != 1)
      {
        algsys_->InitRHS();
        assemble_->AssembleMatrices();
        TS_alg_->UpdateRHS();
      }

      PtrParamNode analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->BuildInDirichlet();

      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
      algsys_->Solve(analysis_id);

      // new solution is NOT only an increment of the full solution =============
      algsys_->GetSolutionVal( newSol_ );
      LOG_DBG(solvestepfluidmech) << "newSol_\n" << newSol_ << std::endl;

      sol_->SetAlgSysVector(newSol_);

      sol_->GetGlobalSolVector(FLUIDMECH_VELOCITY,tmpVelo_);
      sol_->GetGlobalSolVector(FLUIDMECH_PRESSURE,tmpPres_);

      LOG_DBG2(solvestepfluidmech) << "newVelo\n" << tmpVelo_ << std::endl;
      LOG_DBG2(solvestepfluidmech) << "newPres\n" << tmpPres_ << std::endl;

      solVelocityInc_ = tmpVelo_ - actVelo_;
      solPressureInc_ = tmpPres_ - actPres_;
      solIncrement_   = newSol_  - actSol_;

      actVelo_ = tmpVelo_;
      actPres_ = tmpPres_;
      actSol_  = newSol_;

      // calculate incremental error ========================================
      Double solVelocityIncL2Norm = NormL2(solVelocityInc_);
      Double solPressureIncL2Norm = NormL2(solPressureInc_);
      Double solIncrL2Norm = NormL2(solIncrement_);

      Double actVeloL2Norm = NormL2(tmpVelo_);
      Double actPresL2Norm = NormL2(tmpPres_);
      Double actSolL2Norm  = NormL2(actSol_);

      if (actVeloL2Norm > 1.0)
      {
        incrementalVeloErr = solVelocityIncL2Norm / actVeloL2Norm;
      } else {
        incrementalVeloErr = solVelocityIncL2Norm;
      }

      if (actPresL2Norm > 1.0)
      {
        incrementalPresErr = solPressureIncL2Norm / actPresL2Norm;
      } else {
        incrementalPresErr = solPressureIncL2Norm;
      }

      if (actSolL2Norm > 1.0)
      {
        incrementalErr = solIncrL2Norm / actSolL2Norm;
      } else {
        incrementalErr = solIncrL2Norm;
      }

      // boolean variable, holds condition if another iteration step is necessary
      performOneMoreStep = (incrementalErr > incStopCrit_ ||
                            incrementalVeloErr > incStopCrit_ ||
                            incrementalPresErr > incStopCrit_);

#if 1 // INFO: For observation, helps to rely on the solution
      std::cerr << iterationCounter << " Solu = " << actSolL2Norm
        << "; Velo = " << actSolL2Norm
        << "; Pres = " << actSolL2Norm
        << "\tinc Err = " << incrementalErr
        << "; inc Velo Err = " << incrementalVeloErr
        << "; inc Pres Err = " << incrementalPresErr
        << std::endl;
#endif
    }

    if (incrementalErr > 50*incStopCrit_ ||
        incrementalVeloErr > 50*incStopCrit_ ||
        incrementalPresErr > 50*incStopCrit_)
    {
      WARN("FluidMech did not converge!!! Simulation SHOULD be aborted");
    } else {
      if (incrementalErr > incStopCrit_ ||
          incrementalVeloErr > incStopCrit_ ||
          incrementalPresErr > incStopCrit_)
      {
        WARN("FluidMech did not converge!!!");
      }
    }

    if ( isInstationary_ )
    {
      TS_alg_->Corrector(actSol_);
    }

    LOG_DBG(solvestepfluidmech) << "actSol_\n" << actSol_ << std::endl;

    if ( isIterCoupled )
    {
      iterCoupledCounter++;
    }
  }
} // end of namespace

