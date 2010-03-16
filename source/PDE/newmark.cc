// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "OLAS/algsys/basesystem.hh"

#include "PDE/StdPDE.hh"

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "newmark.hh"

namespace CoupledField
{

  Newmark::Newmark(BaseSystem * algebraicsystem, PtrParamNode systemNode )
    :TimeStepping(algebraicsystem )
  {

    // Default values for beta and gamma are:
    //   beta  = 0.25
    //   gamma = 0.5
    alpha_ = 0.0; //  alpha_ = -1/3;
    beta_ = (1-alpha_)*(1-alpha_) / 4.0;
    gamma_ = (1 - 2*alpha_) / 2.0;
    nu_ = 0.0;

    if ( systemNode->Has("timeSteppingParameters") ) {
      PtrParamNode myParam = systemNode->Get("timeSteppingParameters");
      myParam->GetValue("beta", beta_, ParamNode::PASS);
      myParam->GetValue("gamma", gamma_, ParamNode::PASS);
      myParam->GetValue("nu", nu_, ParamNode::PASS);
    }

    gamma_ = gamma_ + nu_;

    Info->PrintF("","In Newmark TimeStepping Scheme use:\n  beta=%f\n  gamma=%f\n",
                 beta_, gamma_);

  }

  Newmark::~Newmark()
  {

  }

  void Newmark::Init( Double dt, UInt rhsSize ) {

    rhsSize_ = rhsSize;

    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
    globalDeltaT_=dt;

    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*a2_;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*a4_;

    if( !isSolTN1Set_){
      sol_tn_1_.Resize(rhsSize_);
      sol_tn_1_.Init();
    }

    //get the memory
    if( !isDeriv1Set_ ) {
      solderiv1_.Resize(rhsSize_);
      solderiv1_.Init();
    }
    if( !isDeriv2Set_ ) {
      solderiv2_.Resize(rhsSize_);
      solderiv2_.Init();
    }

    solpred_.Resize(rhsSize_);
    solpred_.Init();
    solderiv1pred_.Resize(rhsSize_);
    solderiv1pred_.Init();

  }

  void Newmark::setSubSteps(  UInt subSteps ) {

    dt_ = globalDeltaT_/Double(subSteps);
    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*a2_;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*a4_;

    solpredSAVE_ = solpred_;
    solderiv1predSAVE_ = solderiv1pred_;
  }

  void Newmark::resetDeltaT( ) {

    dt_ = globalDeltaT_;
    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*a2_;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*a4_;

    solpred_ = solpredSAVE_;
    solderiv1pred_ = solderiv1predSAVE_;
  }

  void Newmark::Predictor(Vector<Double>& solold)
  {

    solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
 
  }

  void Newmark::UpdateRHS()
  {

    Vector<Double> coeffMass;

    // mass part
    coeffMass = solpred_*a2_;
    algsys_->UpdateRHS(MASS,coeffMass);

    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      Vector<Double> coeffDamp;
        
      coeffDamp = -solderiv1pred_ + solpred_*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp);
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      Vector<Double> coeffStiff;
      coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff);
    }

  }

  void Newmark::SubstractStiffnessFromRHS(Vector<Double>& actSol) {


    Vector<Double> actSolTemp;
    actSolTemp=-actSol;

    algsys_->UpdateRHS(STIFFNESS,actSolTemp);
  }

  void Newmark::UpdateRHS(Vector<Double>& actSol)
  {

    // mass part
    Vector<Double> coeffMass;
    coeffMass = (solpred_ - actSol) * a2_;
    algsys_->UpdateRHS(MASS, coeffMass);


    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      Vector<Double> coeffDamp;
        
      coeffDamp = -solderiv1pred_ + (solpred_-actSol)*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp);
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      Vector<Double> coeffStiff;
      coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff);
    }

  }

  void Newmark::Corrector(Vector<Double>& solnew)
  {

    solderiv2_ = (solnew - solpred_) * a2_;
    solderiv1_ = solderiv1pred_ + solderiv2_*a3_;

    sol_tn_1_=solnew;
  }



  void Newmark::CalcParameters(Double dt)
  {

    //for predictors
    a0_ = 0.5*(1-2.0*beta_)*dt_*dt_;
    a1_ = (1-gamma_)*dt_;

    //for correctors, matrices
    a2_ = 1.0/(beta_*dt_*dt_);
    a3_ = gamma_*dt_;

    //for RHS, Matrices
    a4_ = (1+alpha_) * gamma_ / (beta_*dt_);
  }


  // ====================================================
  // Effective Mass Matrix Formulation
  // ====================================================

  NewmarkEffMass::NewmarkEffMass(BaseSystem * algebraicsystem, PtrParamNode systemNode,
                                 bool intExplicit)
    :TimeStepping(algebraicsystem)
  { 

    intExplicit_ = intExplicit;

    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    nu_ = 0.0;

    if ( systemNode->Has("timeSteppingParameters") ) {
      PtrParamNode myParam = systemNode->Get("timeSteppingParameters");
      myParam->GetValue("beta", beta_, ParamNode::PASS);
      myParam->GetValue("gamma", gamma_, ParamNode::PASS);
      myParam->GetValue("nu", nu_, ParamNode::PASS);
    }

    gamma_ = gamma_ + nu_;

    Info->PrintF("","In Newmark TimeStepping Scheme use:\n  beta=%f\n  gamma=%f\n",
                 beta_, gamma_);
  }

  NewmarkEffMass::~NewmarkEffMass()
  {

  }

  void NewmarkEffMass::Init( Double dt, UInt rhsSize ) 
  {

    std::cout<<"Newmark uses Effective Mass Formulation " <<std::endl;

    rhsSize_ = rhsSize;

    
    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
    CalcParameters(dt_);

    if ( intExplicit_ == true ) {
      matrix_factors_[MASS] = 1.0; 
    }
    else {
      matrix_factors_[STIFFNESS] = 1.0*a2_;
      matrix_factors_[MASS] = 1.0;
      matrix_factors_[CONVECTION] = 0.0;
      
      matrix_factors_[DAMPING] = 1.0*a3_; 
    }

    //get the memory
    sol_.Resize(rhsSize_);
    sol_.Init();
    solpred_.Resize(rhsSize_);
    solpred_.Init();

    if( !isDeriv1Set_ ) {
      solderiv1_.Resize(rhsSize_);
      solderiv1_.Init();
    }

    if( !isDeriv2Set_ ) {
      solderiv2_.Resize(rhsSize_);
      solderiv2_.Init();
    }

  
    solderiv1pred_.Resize(rhsSize_);
    solderiv1pred_.Init();

  }


  void NewmarkEffMass::Predictor(Vector<Double>& solold)
  {

    solpred_ = sol_ + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;


  }

  void NewmarkEffMass::Corrector(Vector<Double>& aNew)
  {

    // after solving the algebraic system of equation, we obtain as solution
    // the 2nd time derivative: aNew .. 2nd second time derivative
    sol_ = solpred_ + aNew * a2_;
    solderiv1_ = solderiv1pred_ + aNew*a3_;
    solderiv2_ = aNew;

    //now overwrite the solution with the physical quantity itself
    aNew = sol_;

  }

  void NewmarkEffMass::UpdateRHS()
  {

    Vector<Double> coeffStiff;
    coeffStiff = -solpred_;

    algsys_->UpdateRHS(STIFFNESS, coeffStiff);

    // damping part
    if ( FeMatrixPresent( DAMPING ) )
      {
        Vector<Double> coeffDamp;
        coeffDamp = -solderiv1pred_;

        algsys_->UpdateRHS(DAMPING, coeffDamp);
      }
  }

  void NewmarkEffMass::CalcParameters(Double dt)
  {

    //for predictors
    a0_ = 0.5*dt_*dt_*(1-2*beta_);
    a1_ = (1-gamma_)*dt_;
  
    //for correctors, matrices
    a2_ = beta_*dt_*dt_;
    a3_ = gamma_*dt_;
  
  }

  Double NewmarkEffMass::DirichletBC4EffMassMatrix(Double val, Integer eq)
  {

    Double accval;

    accval = (val- solpred_[eq-1]) / a2_;
    return accval;
  }



} // end of namespace
