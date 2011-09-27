#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/algsys/algebraicSys.hh"
#include "bdf2.hh"

namespace CoupledField
{
  DECLARE_LOG(bdf2)
  DEFINE_LOG(bdf2, "bdf2")

  Bdf2 :: Bdf2( AlgebraicSys * algebraicsystem)
    :TimeStepping( algebraicsystem )
  {
    // Commented out the warning, since defaults are not bad at all and the
    // average student user gets not disturbed by any warnings
    //check if integration parameters are defined in conf-file
    //Info->WARN( "Bdf2: Using defaults for gamma!" );
    firstTime_=true;
  }

  Bdf2 :: ~Bdf2()
  {
  }

  void Bdf2::Init( Double dt, UInt rhsSize )
  {
    dt_ = dt;
    rhsSize_ = rhsSize;
    SBM_Vector dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();
    coeffMass_ = dummyVec;

    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0;

    //not used matrices
    matrix_factors_[MASS] = 0.0;
    matrix_factors_[CONVECTION] = 0.0;
    matrix_factors_[DAMPING] = 0.0;
    matrix_factors_[AUXILIARY] = 0.0;

    if ( !is_SolTimeStep_set(TIMESTEP_1) )
    {
      sol_timeStepVec_[TIMESTEP_1] = dummyVec;
    }
    if ( !is_SolTimeStep_set(TIMESTEP_2) )
    {
      sol_timeStepVec_[TIMESTEP_2] = dummyVec;
    }
    /* only needed if direct acoustic is calculated */
    //if ( !is_Deriv_set(FIRST_DERIV) )
    //{
    //  solDeriv_vec_[FIRST_DERIV] = dummyVec;
    //}
  }


  void Bdf2::Predictor(SBM_Vector& solold)
  {
    /** szoerner: put into AdvanceTime since bdf is not a predictor corrector
     * method
     * sol_tn_2_ = sol_tn_1_;
     * sol_tn_1_ = solold;
     */
  }


  void Bdf2::UpdateRHS()
  {
    // mass part
    // the last term occurs due to the newton method
    SBM_Vector currentSol;

    algsys_->GetSolutionVal(currentSol);
    algsys_->UpdateRHS(AUXILIARY, currentSol);

    // the minus occurs due to the newton method
    const Double coeff1 = - sol_timeStepCoeff_[TIMESTEP_1];
    const Double coeff2 = - sol_timeStepCoeff_[TIMESTEP_2];
    
    coeffMass_.Add( coeff1, sol_timeStepVec_[TIMESTEP_1],
                    coeff2, sol_timeStepVec_[TIMESTEP_2] );
    
//    coeffMass_  = sol_timeStepVec_[TIMESTEP_1] * coeff1;
//    coeffMass_ += sol_timeStepVec_[TIMESTEP_2] * coeff2;
    algsys_->UpdateRHS(MASS, coeffMass_);
  }


  void Bdf2::Corrector(SBM_Vector& solnew)
  {
    //in solderiv1_ and solderiv2_ for bdf2 not the derivatives are stored.
    // in order to realize restart functionality the last solutions [t_(n-1) t_(n-2)} are stored to retrace a simulation
    //solderiv1_=solnew;
    //solderiv2_=sol_tn_1_;
  }

  void Bdf2::AdvanceTimestep(SBM_Vector& solnew)
  {
    // TODO: were should the derivative be calculated??? Necessary only for
    // calculating acoustics
    //solderiv1_ =  (solnew - sol_tn_1_) * (1.0/dt_);

    sol_timeStepVec_[TIMESTEP_2] = sol_timeStepVec_[TIMESTEP_1];
    sol_timeStepVec_[TIMESTEP_1] = solnew;
  }

  void Bdf2 :: CalcParameters(Double dt)
  {
    //for u_{n}
    sol_timeStepCoeff_[TIMESTEP_0] = 3.0 / (2.0 * dt_); // 1.0
    //for u_{n-1}
    sol_timeStepCoeff_[TIMESTEP_1] = - 2.0/dt_;
    //a0_ = 4.0/3.0;

    //for u_{n-2}
    sol_timeStepCoeff_[TIMESTEP_2] = 1.0/(2.0*dt_);
    //a1_ = 1.0/3.0;
  }


} // end of namespace
