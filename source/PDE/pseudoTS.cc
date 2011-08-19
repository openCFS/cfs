#include <fstream>
#include <iostream>
#include <string>

//#include "DataInOut/WriteInfo.hh"
#include "pseudoTS.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{
  // declare logging stream
//  DECLARE_LOG(pseudoTS)
//  DEFINE_LOG(pseudoTS, "pseudoTS")

  PseudoTS::PseudoTS(BaseSystem * algebraicsystem )
    :TimeStepping(algebraicsystem )
  {
  }

  PseudoTS::~PseudoTS()
  {
  }

  void PseudoTS::Init( Double dt, UInt rhsSize ) {

    rhsSize_ = rhsSize;
    Vector<Double> dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();

    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
        
    matrix_factors_[STIFFNESS]  = 1.0;
    matrix_factors_[MASS]       = 0.0;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING]    = 0.0;

    //get the memory
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

  void PseudoTS::Predictor(Vector<Double>& solold)
  {
  }

  void PseudoTS::Corrector(Vector<Double>& solnew)
  {
    const Double inv_dt_2 = 1.0 / ( dt_ * 2.0);
    solDeriv_vec_[FIRST_DERIV] = (solnew * 3.0 - sol_timeStepVec_[TIMESTEP_1] * 4.0 + \
        sol_timeStepVec_[TIMESTEP_2] ) * inv_dt_2;
  }

  void PseudoTS::AdvanceTimestep(Vector<Double>& solnew)
  {
    sol_timeStepVec_[TIMESTEP_2] = sol_timeStepVec_[TIMESTEP_1];
    sol_timeStepVec_[TIMESTEP_1] = solnew;
  }

} // end of namespace
