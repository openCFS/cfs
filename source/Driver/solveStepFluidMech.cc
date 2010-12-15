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
    if ( isInstationary_ )
    {
      if ( TS_alg_== NULL )
      {
        EXCEPTION( "TS_alg has NULL-Pointer, in SolveStepMag::StepTransNonLin");
      } else {
        if ( isIterCoupled == false || iterCoupledCounter == 0 )
        {
          TS_alg_->Predictor(actSol);
        }
      }
    }

    PDE_.SetBCs();

    UInt iterationCounter=0;

    PDE_.InitStabParams();
    //Update RHS (mass matrix on right hand side)
    if ( isInstationary_ )
    {
      assemble_->AssembleMatrices();
      TS_alg_->UpdateRHS();
      PDE_.SetBCs();
    }

    while (performOneMoreStep && iterationCounter < nonLinMaxIter_)
    {
      iterationCounter++;

      // setup and solve new system (rhs is already set) =====================
      if ( !isInstationary_ || iterationCounter != 1)
      {
        algsys_->InitRHS();
        assemble_->AssembleMatrices();
        TS_alg_->UpdateRHS();
        PDE_.SetBCs();
      }

      PtrParamNode analysis_id = BaseDriver::CreateAnalysisIdChild(analysis_base, "nonLin", iterationCounter);

      algsys_->ConstructEffectiveMatrix(matrix_factor_);
      algsys_->BuildInDirichlet();

      algsys_->SetupPrecond();
      algsys_->SetupSolver(analysis_id);
      algsys_->Solve(analysis_id);

      // new solution is NOT only an increment of the full solution =============
      algsys_->GetSolutionVal( newSol );
      LOG_DBG(solvestepfluidmech) << "newSol\n" << newSol << std::endl;

      sol_->SetAlgSysVector(newSol);

      sol_->GetGlobalSolVector(FLUIDMECH_VELOCITY,tmpVelo);
      sol_->GetGlobalSolVector(FLUIDMECH_PRESSURE,tmpPres);

      LOG_DBG2(solvestepfluidmech) << "newVelo\n" << tmpVelo << std::endl;
      LOG_DBG2(solvestepfluidmech) << "newPres\n" << tmpPres << std::endl;

      solVelocityInc = tmpVelo - actVelo;
      solPressureInc = tmpPres - actPres;
      solIncrement   = newSol  - actSol;

      actVelo = tmpVelo;
      actPres = tmpPres;
      actSol  = newSol;

      // calculate incremental error ========================================
      Double solVelocityIncL2Norm = NormL2(solVelocityInc);
      Double solPressureIncL2Norm = NormL2(solPressureInc);
      Double solIncrL2Norm = NormL2(solIncrement);

      Double actVeloL2Norm = NormL2(tmpVelo);
      Double actPresL2Norm = NormL2(tmpPres);
      Double actSolL2Norm  = NormL2(actSol);

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
      TS_alg_->Corrector(actSol);
    }

    LOG_DBG(solvestepfluidmech) << "actSol\n" << actSol << std::endl;

    if ( isIterCoupled )
    {
      iterCoupledCounter++;
    }
    PDE_.AcouSourceCalc();
  }
} // end of namespace

