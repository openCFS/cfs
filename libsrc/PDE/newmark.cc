#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/WriteInfo.hh"
#include "newmark.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

Newmark :: Newmark(std::string apdename, BaseSystem * algebraicsystem, NodeEQN * ptEQN, 
		   Boolean needDampingMatrix)
:TimeStepping(apdename, algebraicsystem,ptEQN)
{
  ENTER_FCN( "Newmark::Newmark" );

  damping_ = needDampingMatrix;

  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  std::string analysis;
  params->Get( "type", analysis, "analysis" );
  if(analysis != "paramIdent")
    Info->Warning( "Newmark: Using defaults for alpha, beta and gamma!" );

  Integer numEQNs = ptEQN_->GetNumEQNs();
  Integer dofs = ptEQN_->GetNumDofsPerEQN();
  
  //get the memory
  solderiv1_.Resize(numEQNs * dofs);
  solderiv1_.Init();
  solderiv2_.Resize(numEQNs * dofs);
  solderiv2_.Init();

  solpred_.Resize(numEQNs * dofs);
  solpred_.Init();
  solderiv1pred_.Resize(numEQNs * dofs);
  solderiv1pred_.Init();

}

Newmark :: ~Newmark()
{
  ENTER_FCN( "Newmark::~Newmark" );

}

void Newmark::Init(Double * matrix_factors, Double dt)
{
  ENTER_FCN( "Newmark::Init" );

  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 
  if (damping_)
    matrix_factors[1] = 1.0*a4_; // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0*a2_;   // factor for mass matrix

}


void Newmark::Predictor(Vector<Double>& solold)
{
  ENTER_FCN( "Newmark::Predictor" );

  solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
#ifdef DEBUG
  (*debug) << "Predictor:\n";
  (*debug) << "solpred: \n" << solpred_ << std::endl;
  (*debug) << " solderiv1pred: \n" <<  solderiv1pred_ << std::endl;
#endif
}


void Newmark::UpdateRHS()
{
  ENTER_FCN( "Newmark::UpdateRHS" );

  Vector<Double> coeffMass;

  // mass part
  coeffMass = solpred_*a2_;
  algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

//   // additional term for thermoviscous damping model
//   if ( dampType_ == THERMOVISCOUS ) {

// 	coeffDamp = -solderiv1pred_ + solpred_*a4_; // damping part
// 	algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
//   }

  // damping part
  if (damping_) 
    {
      Vector<Double> coeffDamp;

      coeffDamp = -solderiv1pred_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());

      coeffDamp = solpred_*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
   }
}





void Newmark::UpdateRHS(Vector<Double>& actSol)
{
  ENTER_FCN( "Newmark::UpdateRHS" );

  // mass part
  Vector<Double> coeffMass;
  coeffMass = (solpred_ - actSol) * a2_;
  algsys_->UpdateRHS(MASS, coeffMass.GetPointer());


  // damping part
  if (damping_) 
    {
      Vector<Double> coeffDamp;

      coeffDamp = -solderiv1pred_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());

      coeffDamp = (solpred_-actSol)*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
   }
}





void Newmark::Corrector(Vector<Double>& solnew)
{
  ENTER_FCN( "Newmark::Corrector" );

  solderiv2_ = (solnew - solpred_) * a2_;
  solderiv1_ = solderiv1pred_ + solderiv2_*a3_;
#ifdef DEBUG
  (*debug) << "Corrector:\n";
  (*debug) << "sol: \n" << solnew << std::endl;
  (*debug) << "solderiv1: \n" <<  solderiv1_ << std::endl;
  (*debug) << "solderiv2: \n" <<  solderiv2_ << std::endl;
#endif
}



void Newmark :: CalcParameters(Double dt)
{
  ENTER_FCN( "Newmark::CalcParameters" );

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
				 NodeEQN * ptEQN, 
				 Boolean adamping)
  :TimeStepping(apdename, algebraicsystem,ptEQN)
{ 
  ENTER_FCN( "NewmarkEffMass::NewmarkEffMass" );

  damping_ = adamping;

  alpha_ = 0.0;
  beta_  = 0.25;
  gamma_ = 0.5;

  //check if integration parameters are defined in conf-file
  std::string analysis;
  params->Get( "type", analysis, "analysis" );
  if(analysis != "paramIdent")
    Info->Warning( "Newmark: Using defaults for alpha, beta and gamma!" );

  Integer numEQNs = ptEQN_->GetNumEQNs();
  Integer dofs = ptEQN_->GetNumDofsPerEQN();


  //get the memory
  sol_.Resize(numEQNs*dofs);
  sol_.Init();
  solpred_.Resize(numEQNs*dofs);
  solpred_.Init();

  solderiv1_.Resize(numEQNs*dofs);
  solderiv1_.Init();
  solderiv2_.Resize(numEQNs*dofs);
  solderiv2_.Init();

  
  solderiv1pred_.Resize(numEQNs*dofs);
  solderiv1pred_.Init();
}



NewmarkEffMass :: ~NewmarkEffMass()
{
  ENTER_FCN( "NewmarkEffMass::~NewmarkEffMass" );

}

void NewmarkEffMass::Init(Double * matrix_factors, Double dt)
{
  ENTER_FCN( "NewmarkEffMass::Init" );

  dt_ = dt;
  CalcParameters(dt_);

  matrix_factors[0] = 1.0*a2_;       // factor for stiffness matrix

  matrix_factors[1] = 0.0; 
  if (damping_)
    matrix_factors[1] = 1.0*a3_; // factor for damping matrix

  matrix_factors[2] = 0.0;       // factor for convection matrix
  matrix_factors[3] = 1.0;       // factor for mass matrix

}


void NewmarkEffMass::Predictor(Vector<Double>& solold)
{
  ENTER_FCN( "NewmarkEffMass::Predictor" );

  solpred_ = sol_ + solderiv1_*dt_ + solderiv2_*a0_;
  solderiv1pred_ = solderiv1_ + solderiv2_*a1_;

#ifdef DEBUG
  (*debug) << "Predictor:\n";
  (*debug) << "solpred: \n" << solpred_ << std::endl;
  (*debug) << " solderiv1pred: \n" <<  solderiv1pred_ << std::endl;
#endif
}



void NewmarkEffMass::Corrector(Vector<Double>& aNew)
{
  ENTER_FCN( "NewmarkEffMass::Corrector" );

  sol_ = solpred_ + aNew * a2_;
  solderiv1_ = solderiv1pred_ + aNew*a3_;
  solderiv2_ = aNew;
#ifdef DEBUG
  (*debug) << "Corrector:\n";
  (*debug) << "sol: \n" << sol_ << std::endl;
  (*debug) << "solderiv1: \n" <<  solderiv1_ << std::endl;
  (*debug) << "solderiv2: \n" <<  solderiv2_ << std::endl;
#endif
}




void NewmarkEffMass::UpdateRHS()
{
  ENTER_FCN( "NewmarkEffMass::UpdateRHS" );

  Vector<Double> coeffStiff;
  coeffStiff = -solpred_;

  algsys_->UpdateRHS(STIFFNESS, coeffStiff.GetPointer());

  // damping part
  if (damping_)
    {
      Vector<Double> coeffDamp;
      coeffDamp = -solderiv1pred_;

      algsys_->UpdateRHS(DAMPING, coeffDamp.GetPointer());
    }
}



void const NewmarkEffMass :: StoreSolution(Vector<Double> & solArr) const
{
  ENTER_FCN( "NewmarkEffMass::CalcParameters" );

  solArr = sol_; 
 }





void NewmarkEffMass :: CalcParameters(Double dt)
{
  ENTER_FCN( "NewmarkEffMass::CalcParameters" );

  //for predictors
  a0_ = 0.5*dt_*dt_*(1-2*beta_);
  a1_ = (1-gamma_)*dt_;
  
  //for correctors, matrices
  a2_ = beta_*dt_*dt_;
  a3_ = gamma_*dt_;
  
}


} // end of namespace
