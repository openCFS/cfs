#include <fstream>
#include <iostream>
#include <string>

#include "trapezoidal.hh"

namespace CoupledField
{

Trapezoidal :: Trapezoidal(std::string apdename, BaseSystem * algebraicsystem, Integer dofspernode, 
		   Integer numnode)
:TimeStepping(apdename, algebraicsystem)
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::Trapezoidal" << std::endl;
#endif

  gamma_ = 1;

  //check if integration parameters are defined in conf-file
  conf->ifget("gamma_P",gamma_,pdename_);

  //get the memory
  solderiv1_.reshape(dofspernode, numnode);  
  solderiv1_.init();

  solpred_.reshape(dofspernode, numnode); 
  solpred_.init();
}

Trapezoidal :: ~Trapezoidal()
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::~Trapezoidal" << std::endl;
#endif

}

void Trapezoidal::Init(Double * matrix_factors, Double dt)
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::Init" << std::endl;
#endif
  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix
  matrix_factors[3] = a1_;       // factor for mass matrix

  //not used matrices
  matrix_factors[1] = 0.0; 
  matrix_factors[2] = 0.0;       
}


void Trapezoidal::Predictor(Array<Double>& solold)
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::Predictor" << std::endl;
#endif

  solpred_ = solold + solderiv1_*a0_;
}


void Trapezoidal::UpdateRHS()
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::UpdateRHS" << std::endl;
#endif

  Vector<Double> coeffMass;

  // mass part
  coeffMass = solpred_*a1_;
  algsys_->UpdateRHS(MASS,coeffMass.get());
}


void Trapezoidal::Corrector(Array<Double>& solnew)
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::Corrector" << std::endl;
#endif

  solderiv1_ = (solnew - solpred_)*a1_;

}

void Trapezoidal :: CalcParameters(Double dt)
{
#ifdef TRACE
  (*trace) << "entering Trapezoidal::CalcParameters" << std::endl;
#endif

  //for predictors
  a0_ = (1-gamma_)*dt_;

  //for correctors, matrices
  a1_ = 1.0/(gamma_*dt_);
}


} // end of namespace
