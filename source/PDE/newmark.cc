// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "newmark.hh"

namespace CoupledField
{

  Newmark::Newmark(BaseSystem * algebraicsystem, const std::string& sysName )
    :TimeStepping(algebraicsystem )
  {
    ENTER_FCN( "Newmark::Newmark" );

    // Default values for beta and gamma are:
    //   beta  = 0.25
    //   gamma = 0.5
    alpha_ = 0.0; //  alpha_ = -1/3;
    beta_ = (1-alpha_)*(1-alpha_) / 4.0;
    gamma_ = (1 - 2*alpha_) / 2.0;
    nu_ = 0.0;

    ParamNode* myParam = NULL;
    myParam = param->Get("linearSystems", false);
    if( myParam != NULL ) {
      myParam = myParam->Get("system", "name", sysName, false);
      if( myParam != NULL ) {
        if ( param->Has("timeSteppingParameters") ) {
          myParam =  myParam->Get("timeSteppingParameters");
          myParam->Get("beta", beta_, false);
          myParam->Get("gamma", gamma_, false);
          myParam->Get("nu", nu_, false);
        }
      }
    }

    gamma_ = gamma_ + nu_;

    Info->PrintF("","In Newmark TimeStepping Scheme use:\n  beta=%f\n  gamma=%f\n",
                 beta_, gamma_);

  }

  Newmark::~Newmark()
  {
    ENTER_FCN( "Newmark::~Newmark" );

  }

  void Newmark::Init( Double dt, UInt rhsSize ) {
    ENTER_FCN( "Newmark::Init" );

    rhsSize_ = rhsSize;

    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*a2_;   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*a4_;

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

  void Newmark::Predictor(Vector<Double>& solold)
  {
    ENTER_FCN( "Newmark::Predictor" );

    solpred_ = solold + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;
 
  }

  void Newmark::UpdateRHS()
  {
    ENTER_FCN( "Newmark::UpdateRHS" );

    Vector<Double> coeffMass;

    // mass part
    coeffMass = solpred_*a2_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());

    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      Vector<Double> coeffDamp;
        
      coeffDamp = -solderiv1pred_ + solpred_*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      Vector<Double> coeffStiff;
      coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff.GetPointer());
    }

  }

  void Newmark::SubstractStiffnessFromRHS(Vector<Double>& actSol) {

    ENTER_FCN( "Newmark::SubstractStiffness" );

    Vector<Double> actSolTemp;
    actSolTemp=-actSol;

    algsys_->UpdateRHS(STIFFNESS,actSolTemp.GetPointer());
  }

  void Newmark::UpdateRHS(Vector<Double>& actSol)
  {
    ENTER_FCN( "Newmark::UpdateRHS" );

    // mass part
    Vector<Double> coeffMass;
    coeffMass = (solpred_ - actSol) * a2_;
    algsys_->UpdateRHS(MASS, coeffMass.GetPointer());


    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      Vector<Double> coeffDamp;
        
      coeffDamp = -solderiv1pred_ + (solpred_-actSol)*a4_;
      algsys_->UpdateRHS(DAMPING,coeffDamp.GetPointer());
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      Vector<Double> coeffStiff;
      coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff.GetPointer());
    }

  }

  void Newmark::Corrector(Vector<Double>& solnew)
  {
    ENTER_FCN( "Newmark::Corrector" );

    solderiv2_ = (solnew - solpred_) * a2_;
    solderiv1_ = solderiv1pred_ + solderiv2_*a3_;
  }



  void Newmark::CalcParameters(Double dt)
  {
    ENTER_FCN( "Newmark::CalcParameters" );

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

  NewmarkEffMass::NewmarkEffMass(BaseSystem * algebraicsystem, const std::string& sysName,
                                 bool intExplicit)
    :TimeStepping(algebraicsystem)
  { 
    ENTER_FCN( "NewmarkEffMass::NewmarkEffMass" );

    intExplicit_ = intExplicit;

    alpha_ = 0.0;
    beta_  = 0.25;
    gamma_ = 0.5;

    nu_ = 0.0;


    ParamNode* myParam = NULL;
    myParam = param->Get("linearSystems", false);
    if( myParam != NULL ) {
      myParam = myParam->Get("system", "name", sysName, false);
      if( myParam != NULL ) {
        if ( param->Has("timeSteppingParameters") ) {
          myParam =  myParam->Get("timeSteppingParameters");
          myParam->Get("beta", beta_, false);
          myParam->Get("gamma", gamma_, false);
          myParam->Get("nu", nu_, false);
        }
      }
    }

    gamma_ = gamma_ + nu_;

    Info->PrintF("","In Newmark TimeStepping Scheme use:\n  beta=%f\n  gamma=%f\n",
                 beta_, gamma_);
  }

  NewmarkEffMass::~NewmarkEffMass()
  {
    ENTER_FCN( "NewmarkEffMass::~NewmarkEffMass" );

  }

  void NewmarkEffMass::Init( Double dt, UInt rhsSize ) 
  {
    ENTER_FCN( "NewmarkEffMass::Init" );

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
    ENTER_FCN( "NewmarkEffMass::Predictor" );

    solpred_ = sol_ + solderiv1_*dt_ + solderiv2_*a0_;
    solderiv1pred_ = solderiv1_ + solderiv2_*a1_;


  }

  void NewmarkEffMass::Corrector(Vector<Double>& aNew)
  {
    ENTER_FCN( "NewmarkEffMass::Corrector" );

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
    ENTER_FCN( "NewmarkEffMass::UpdateRHS" );

    Vector<Double> coeffStiff;
    coeffStiff = -solpred_;

    algsys_->UpdateRHS(STIFFNESS, coeffStiff.GetPointer());

    // damping part
    if ( FeMatrixPresent( DAMPING ) )
      {
        Vector<Double> coeffDamp;
        coeffDamp = -solderiv1pred_;

        algsys_->UpdateRHS(DAMPING, coeffDamp.GetPointer());
      }
  }

  void NewmarkEffMass::CalcParameters(Double dt)
  {
    ENTER_FCN( "NewmarkEffMass::CalcParameters" );

    //for predictors
    a0_ = 0.5*dt_*dt_*(1-2*beta_);
    a1_ = (1-gamma_)*dt_;
  
    //for correctors, matrices
    a2_ = beta_*dt_*dt_;
    a3_ = gamma_*dt_;
  
  }

  Double NewmarkEffMass::DirichletBC4EffMassMatrix(Double val, Integer eq)
  {
    ENTER_FCN( "NewmarkEffMass::DirichletBC4EffMassMatrix" );

    Double accval;

    accval = (val- solpred_[eq-1]) / a2_;
    return accval;
  }



} // end of namespace
