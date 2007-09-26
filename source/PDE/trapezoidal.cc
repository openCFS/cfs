// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "trapezoidal.hh"

namespace CoupledField
{

  Trapezoidal :: Trapezoidal( BaseSystem * algebraicsystem)
    :TimeStepping( algebraicsystem )
  {

    gamma_ = 1;

    // Commented out the warning, since defaults are not bad at all and the 
    // average student user gets not disturbed by any warnings
    //check if integration parameters are defined in conf-file
    //Info->Warning( "Trapezoidal: Using defaults for gamma!" );


  }

  Trapezoidal :: ~Trapezoidal()
  {

  }

  void Trapezoidal::Init( Double dt, UInt rhsSize ) {
    
    dt_ = dt;
    rhsSize_ = rhsSize;
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0;
    matrix_factors_[MASS] = a1_;

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 0.0;       

    //get the memory
    if( !isDeriv1Set_ ) {
      solderiv1_.Resize(rhsSize_);  
      solderiv1_.Init();
    }

    solpred_.Resize(rhsSize_); 
    solpred_.Init();
  }


  void Trapezoidal::Predictor(Vector<Double>& solold)
  {

    solpred_ = solold + solderiv1_*a0_;
  }


  void Trapezoidal::UpdateRHS()
  {

    Vector<Double> coeffMass;

    // mass part
    coeffMass = solpred_*a1_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());
  }


  void Trapezoidal::UpdateRHS(Vector<Double>& actSol)
  {

    Vector<Double> coeffMass;

    // mass part
    coeffMass = (solpred_ - actSol) *a1_;
    algsys_->UpdateRHS(MASS,coeffMass.GetPointer());
  }


  void Trapezoidal::Corrector(Vector<Double>& solnew)
  {

    solderiv1_ = (solnew - solpred_)*a1_;
  }

  void Trapezoidal :: CalcParameters(Double dt)
  {

    //for predictors
    a0_ = (1-gamma_)*dt_;

    //for correctors, matrices
    a1_ = 1.0/(gamma_*dt_);
  }



  //====================================================================
  //---------------------------- Effective Mass ------------------
  //====================================================================

  TrapezoidalEffMass :: TrapezoidalEffMass( BaseSystem * algebraicsystem)
    :TimeStepping( algebraicsystem )
  {

    gamma_ = 0.51;

  }

  TrapezoidalEffMass :: ~TrapezoidalEffMass()
  {

  }

  void TrapezoidalEffMass::Init( Double dt, UInt rhsSize ) {
    
    dt_ = dt;
    rhsSize_ = rhsSize;
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = a1_;
    matrix_factors_[MASS] = 1.0;

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0;
    matrix_factors_[DAMPING] = 0.0;

    //get the memory
    sol_.Resize(rhsSize_);
    sol_.Init();

    solderiv1_.Resize(rhsSize_);
    solderiv1_.Init();

    solpred_.Resize(rhsSize_); 
    solpred_.Init();

  }


  void TrapezoidalEffMass::Predictor(Vector<Double>& solold)
  {
    solpred_ = solold + solderiv1_*a0_;
  }


  void TrapezoidalEffMass::UpdateRHS()
  {

    Vector<Double> coeffStiff;

    // mass part
    coeffStiff = -solpred_;
    algsys_->UpdateRHS(STIFFNESS,coeffStiff.GetPointer());
  }


  void TrapezoidalEffMass::UpdateRHS(Vector<Double>& actSol)
  {
    UpdateRHS();
  }


  void TrapezoidalEffMass::Corrector(Vector<Double>& vNew)
  {
    // after solving the algebraic system of equation, we obtain as solution
    // the 1st time derivative: vNew .. 1st time derivative
    sol_ = solpred_ + vNew * a1_;
    solderiv1_ = vNew;

    //now overwrite the solution with the physical quantity itself
    vNew = sol_;
  }

  void TrapezoidalEffMass :: CalcParameters(Double dt)
  {

    //for predictors
    a0_ = (1-gamma_)*dt;

    //for correctors, matrices
    a1_ = (gamma_*dt);
  }


  Double TrapezoidalEffMass::DirichletBC4EffMassMatrix(Double val, Integer eq)
  {

    Double velVal;

    velVal = (val- solpred_[eq-1]) / a1_;
    return velVal;
  }




} // end of namespace
