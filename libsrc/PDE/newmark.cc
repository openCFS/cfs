#include <fstream>
#include <iostream>
#include <string>

#include "newmark.hh"

namespace CoupledField
{

Newmark :: Newmark(std::string apdename, BaseSystem * algebraicsystem, Integer dofspernode, 
		   Integer numnode, DampingType adamping)
:TimeStepping(apdename, algebraicsystem)
{
#ifdef TRACE
  (*trace) << "entering Newmark::Newmark" << std::endl;
#endif

  damping_ = adamping;

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

}

Newmark :: ~Newmark()
{
#ifdef TRACE
  (*trace) << "entering Newmark::~Newmark" << std::endl;
#endif

}

void Newmark::Init(Double * matrix_factors, Double dt)
{
#ifdef TRACE
  (*trace) << "entering Newmark::Init" << std::endl;
#endif
  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 
  if (damping_)
    matrix_factors[1] = 1.0*a4_; // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0*a2_;   // factor for mass matrix

}


void Newmark::Predictor(Array<Double>& solold)
{
#ifdef TRACE
  (*trace) << "entering Newmark::Predictor" << std::endl;
#endif

  solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
}


void Newmark::UpdateRHS()
{
#ifdef TRACE
  (*trace) << "entering Newmark::UpdateRHS" << std::endl;
#endif

  Vector<Double> coeffMass;

  // mass part
  coeffMass = solpred_*a2_;
  algsys_->UpdateRHS(MASS,coeffMass.get());

  // damping part
  if (damping_) 
    {
      Vector<Double> coeffDamp;

      coeffDamp = -solderiv1pred_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());

      coeffDamp = solpred_*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());
   }
}





void Newmark::UpdateRHS(std::vector<Double>& actSol)
{
#ifdef TRACE
  (*trace) << "entering Newmark::UpdateRHS" << std::endl;
#endif

  // mass part
  Vector<Double> coeffMass;
  coeffMass = (solpred_ - actSol) * a2_;
  algsys_->UpdateRHS(MASS, coeffMass.get());


  // damping part
  if (damping_) 
    {
      Vector<Double> coeffDamp;

      coeffDamp = -solderiv1pred_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());

      coeffDamp = (solpred_-actSol)*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.get());
   }
}





void Newmark::Corrector(Array<Double>& solnew)
{
#ifdef TRACE
  (*trace) << "entering Newmark::Corrector" << std::endl;
#endif

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;
}



void Newmark :: CalcParameters(Double dt)
{
#ifdef TRACE
  (*trace) << "entering Newmark::CalcParameters" << std::endl;
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
