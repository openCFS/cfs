#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "newmarkdamp.hh"

namespace CoupledField
{

NewmarkDamp :: NewmarkDamp(std::string apdename, BaseSystem * algebraicsystem, Integer dofspernode, 
		   Integer numnode, Integer adamping, Integer afrac_memory, Double damp)
:TimeStepping(apdename, algebraicsystem)
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::NewmarkDamp" << std::endl;
#endif

  damping_ = adamping;
  frac_memory_ = afrac_memory;
  y_ = damp;

  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  conf->ifget("alpha_NM",alpha_,pdename_); 
  conf->ifget("beta_NM",beta_,pdename_); 
  conf->ifget("gamma_NM",gamma_,pdename_);

  //get the memory
  solderiv1_.reshape(dofspernode, numnode);  
  solderiv1_.init();
  solderiv2_.reshape(dofspernode, numnode);  
  solderiv2_.init();
  

  solpred_.reshape(dofspernode, numnode); 
  solpred_.init();
  solderiv1pred_.reshape(dofspernode, numnode); 
  solderiv1pred_.init();
  if (frac_memory_==0)
    Error("Attenuation model needs frac_memory larger than 0",__FILE__,__LINE__);

  //get the memory
  solfrac_ = new Array<Double>[frac_memory_];
  for (Integer i=0; i<frac_memory_; i++)
    {
      solfrac_[i].reshape(dofspernode, numnode);  
      solfrac_[i].init();
    }
}

NewmarkDamp :: ~NewmarkDamp()
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::~NewmarkDamp" << std::endl;
#endif

}

void NewmarkDamp::Init(Double * matrix_factors, Double dt)
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::Init" << std::endl;
#endif
  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 
  if (damping_)
    matrix_factors[1] = 1/pow(dt_,(y_+1.0)); // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0*a2_;   // factor for mass matrix

}


void NewmarkDamp::Predictor(Array<Double>& solold)
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::Predictor" << std::endl;
#endif

  solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
}


void NewmarkDamp::UpdateRHS()
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::UpdateRHS" << std::endl;
#endif

  Vector<Double> coeffMass, coeffDamp;

  // mass part
  coeffMass = solpred_*a2_;
  algsys_->UpdateRHS(MASS,coeffMass.get());

  // damping part
  Double coeff;
  coeff = 1;
  for (Integer i=0; i<frac_memory_; i++)
    {
      coeff *= (i - y_-1)/(i+1);
      coeffDamp = solfrac_[i] * coeff/pow(dt_,(y_+1));
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());
    }
   
}


void NewmarkDamp::Corrector(Array<Double>& solnew)
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::Corrector" << std::endl;
#endif

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

  //solfrac[0]: p_n; solfrac[1]: p_(n-1); ..
  for (Integer i=frac_memory_-1; i>0; i--)
    solfrac_[i] = solfrac_[i-1];
      
  solfrac_[0] = solnew;
}

void NewmarkDamp :: CalcParameters(Double dt)
{
#ifdef TRACE
  (*trace) << "entering NewmarkDamp::CalcParameters" << std::endl;
#endif

  //for predictors
  a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
  a1_ = (1-gamma_)*dt_;

  //for correctors, matrices
  a2_ = 1.0/(beta_*dt_*dt_);
  a3_ = gamma_*dt_;

  //for RHS, Matrices
  a4_ = gamma_ / (beta_*dt_);
}


} // end of namespace
