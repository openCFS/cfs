#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/algsys/basesystem.hh"
#include "bdf2.hh"

namespace CoupledField
{
  DECLARE_LOG(bdf2)
  DEFINE_LOG(bdf2, "bdf2")

  Bdf2 :: Bdf2( BaseSystem * algebraicsystem)
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
    Vector<Double> dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0; //dt_*2.0/3.0
    sol_timeStepCoeff_[COEFFRHS] = 3.0 / (2.0*dt_); // 1.0
    matrix_factors_[MASS] = sol_timeStepCoeff_[COEFFRHS]; // 1.0

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 0.0;       

    if ( !is_SolTimeStep_set(TIMESTEP_1) )
    {
      sol_timeStepVec_[TIMESTEP_1] = dummyVec;
    }
    if ( !is_SolTimeStep_set(TIMESTEP_2) )
    {
      sol_timeStepVec_[TIMESTEP_2] = dummyVec;
    }
    if ( !is_Deriv_set(FIRST_DERIV) )
    {
      solDeriv_vec_[FIRST_DERIV] = dummyVec;
    }
  }


  void Bdf2::Predictor(Vector<Double>& solold)
  {
  }


  void Bdf2::UpdateRHS()
  {
    Vector<Double> coeffMass;

    // mass part
    coeffMass = (sol_timeStepVec_[TIMESTEP_1] * sol_timeStepCoeff_[TIMESTEP_1]) \
                - (sol_timeStepVec_[TIMESTEP_2] * sol_timeStepCoeff_[TIMESTEP_2]);
    algsys_->UpdateRHS(MASS,coeffMass);
    LOG_DBG(bdf2) << "\n coeffMass: \n" << coeffMass << std::endl;
  }


  void Bdf2::Corrector(Vector<Double>& solnew)
  {
    //in solderiv1_ and solderiv2_ for bdf2 not the derivatives are stored.
    // in order to realize restart functionality the last solutions [t_(n-1) t_(n-2)} are stored to retrace a simulation
    //solderiv1_=solnew;
    //solderiv2_=sol_tn_1_;

    //std::cout << "\n Bdf2 Cor solnew: \n" << solnew << std::endl;
    //std::cout << "\n Bdf2 Cor solderiv2_: \n" << solderiv2_ << std::endl;
    //LOG_DBG(bdf2) << "\n Cor solnew: \n" << solnew << std::endl;
    //LOG_DBG(bdf2) << "\n Cor solderiv2_: \n" << solderiv2_ << std::endl;
  }

  void Bdf2::AdvanceTimestep(Vector<Double>& solnew)
  {
    // TODO: were should the derivative be calculated??? Necessary only for
    // calculating acoustics
    //solderiv1_ =  (solnew - sol_tn_1_) * (1.0/dt_);

    sol_timeStepVec_[TIMESTEP_2] = sol_timeStepVec_[TIMESTEP_1];
    sol_timeStepVec_[TIMESTEP_1] = solnew;
  }

  void Bdf2 :: CalcParameters(Double dt)
  {
    //for u_{n-1}
    sol_timeStepCoeff_[TIMESTEP_1] = 2.0/dt_;
    //a0_ = 4.0/3.0;

    //for u_{n-2}
    sol_timeStepCoeff_[TIMESTEP_2] = 1.0/(2.0*dt_);
    //a1_ = 1.0/3.0;
  }


} // end of namespace
