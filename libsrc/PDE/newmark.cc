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


// ====================================================
// Effective Mass Matrix Formulation
// ====================================================



NewmarkEffMass :: NewmarkEffMass(std::string apdename, 
				 BaseSystem * algebraicsystem, 
				 Integer dofspernode, 
				 Integer numnode, 
				 DampingType adamping)
  :TimeStepping(apdename, algebraicsystem)
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::NewmarkEffMass" << std::endl;
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
  sol_.reshape(dofspernode, numnode);  
  sol_.init();
  solpred_.reshape(dofspernode, numnode); 
  solpred_.init();

  solderiv1_.reshape(dofspernode, numnode);  
  solderiv1_.init();
  solderiv1pred_.reshape(dofspernode, numnode); 
  solderiv1pred_.init();  

  solderiv2_.reshape(dofspernode, numnode);  
  solderiv2_.init();
}



NewmarkEffMass :: ~NewmarkEffMass()
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::~NewmarkEffMass" << std::endl;
#endif

}

void NewmarkEffMass::Init(Double * matrix_factors, Double dt)
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::Init" << std::endl;
#endif
  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0*a2_;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 
  if (damping_)
    matrix_factors[1] = 1.0*a3_; // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0;       // factor for mass matrix

}


void NewmarkEffMass::Predictor(Array<Double>& solold)
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::Predictor" << std::endl;
#endif

  solpred_ = sol_ + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
}



void NewmarkEffMass::Corrector(Array<Double>& aNew)
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::Corrector" << std::endl;
#endif

  sol_ = solpred_ + aNew * a2_;
  solderiv1_ = solderiv1pred_ + aNew*a3_;
  solderiv2_ = aNew;
}




void NewmarkEffMass::UpdateRHS()
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::UpdateRHS" << std::endl;
#endif

  Vector<Double> coeffStiff;
  coeffStiff = -solpred_;

  algsys_->UpdateRHS(STIFFNESS, coeffStiff.get());

  // damping part
  if (damping_)
    {
      Vector<Double> coeffDamp;
      coeffDamp = -solderiv1pred_;

      algsys_->UpdateRHS(DAMPING, coeffDamp.get());
    }
}



void const NewmarkEffMass :: StoreSol(Array<Double> & solArr) const
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::CalcParameters" << std::endl;
#endif

  for(Integer pos=0; pos<solArr.size(); pos++) 
    for(Integer actDim=0; actDim<solArr.dim(); actDim++) 
      solArr[actDim][pos] = sol_[0][pos * solArr.dim() + actDim];
}





void NewmarkEffMass :: CalcParameters(Double dt)
{
#ifdef TRACE
  (*trace) << "entering NewmarkEffMass::CalcParameters" << std::endl;
#endif

  //for predictors
  a0_ = 0.5*dt_*dt_*(1-2*beta_);
  a1_ = (1-gamma_)*dt_;
  
  //for correctors, matrices
  a2_ = beta_*dt_*dt_;
  a3_ = gamma_*dt_;
  
}


} // end of namespace
