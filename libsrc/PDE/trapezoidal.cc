#include <fstream>
#include <iostream>
#include <string>

#include "trapezoidal.hh"

namespace CoupledField
{

  Trapezoidal :: Trapezoidal(std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN)
:TimeStepping(apdename, algebraicsystem, ptEQN)
{
  ENTER_FCN( "Trapezoidal::Trapezoidal" );

  gamma_ = 1;

  // Commented out the warning, since defaults are not bad at all and the 
  // average student user gets not disturbed by any warnings
  //check if integration parameters are defined in conf-file
  //Info->Warning( "Trapezoidal: Using defaults for gamma!" );

  Integer numEQNs = ptEQN_->GetNumEQNs();
  Integer dofs = ptEQN_->GetNumDofsPerEQN();

  //get the memory
  solderiv1_.Resize(numEQNs * dofs);  
  solderiv1_.Init();

  solpred_.Resize(numEQNs * dofs); 
  solpred_.Init();
}

Trapezoidal :: ~Trapezoidal()
{
  ENTER_FCN( "Trapezoidal::~Trapezoidal" );

}

void Trapezoidal::Init(Double * matrix_factors, Double dt)
{
  ENTER_FCN( "Trapezoidal::Init" );

  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix
  matrix_factors[3] = a1_;       // factor for mass matrix

  //not used matrices
  matrix_factors[1] = 0.0; 
  matrix_factors[2] = 0.0;       
}


void Trapezoidal::Predictor(Vector<Double>& solold)
{
  ENTER_FCN( "Trapezoidal::Predictor" );

  solpred_ = solold + solderiv1_*a0_;
}


void Trapezoidal::UpdateRHS()
{
  ENTER_FCN( "Trapezoidal::UpdateRHS" );

  Vector<Double> coeffMass;

  // mass part
  coeffMass = solpred_*a1_;
  algsys_->UpdateRHS(MASS,coeffMass.GetPointer());
}


void Trapezoidal::UpdateRHS(Vector<Double>& actSol)
{
  ENTER_FCN( "Trapezoidal::UpdateRHS" );

  Vector<Double> coeffMass;

  // mass part
  coeffMass = (solpred_ - actSol) *a1_;
  algsys_->UpdateRHS(MASS,coeffMass.GetPointer());
}


void Trapezoidal::Corrector(Vector<Double>& solnew)
{
  ENTER_FCN( "Trapezoidal::Corrector" );

  solderiv1_ = (solnew - solpred_)*a1_;
}

void Trapezoidal :: CalcParameters(Double dt)
{
  ENTER_FCN( "Trapezoidal::CalcParameters" );

  //for predictors
  a0_ = (1-gamma_)*dt_;

  //for correctors, matrices
  a1_ = 1.0/(gamma_*dt_);
}


} // end of namespace
