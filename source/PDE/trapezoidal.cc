// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "trapezoidal.hh"

#include "OLAS/algsys/basesystem.hh"

namespace CoupledField
{

  Trapezoidal::Trapezoidal( BaseSystem * algebraicsystem, 
                            PtrParamNode systemNode )
    :TimeStepping( algebraicsystem )
  {

    gamma_ = 1;

    // Commented out the warning, since defaults are not bad at all and the 
    // average student user gets not disturbed by any warnings
    //check if integration parameters are defined in conf-file
    //Info->WARN( "Trapezoidal: Using defaults for gamma!" );
    if ( systemNode->Has("timeSteppingParameters") ) {
      PtrParamNode myParam = systemNode->Get("timeSteppingParameters");
      myParam->GetValue("omitInitialSol", omitFirstPredictor_, ParamNode::PASS);
      myParam->GetValue("gamma", gamma_, ParamNode::PASS);
    }

  }

  Trapezoidal::~Trapezoidal()
  {

  }

  void Trapezoidal::Init( Double dt, UInt rhsSize ) {
    
    dt_ = dt;
    rhsSize_ = rhsSize;
    Vector<Double> dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0;
    matrix_factors_[MASS] = sol_timeStepCoeff_[CORRECTOR_1];

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 0.0;       

    //get the memory
    if ( !is_SolTimeStep_set(TIMESTEP_1) )
    {
      sol_timeStepVec_[TIMESTEP_1] = dummyVec;
    }
    if ( !is_Deriv_set(FIRST_DERIV) )
    {
      solDeriv_vec_[FIRST_DERIV] = dummyVec;
    }
    if ( !is_Deriv_set(SECOND_DERIV) )
    {
      solDeriv_vec_[SECOND_DERIV] = dummyVec;
    }

    solpred_.Resize(rhsSize_); 
    solpred_.Init();
  }


  void Trapezoidal::Predictor(Vector<Double>& solold)
  {
    if( !omitFirstPredictor_) {
      solpred_ = solold + solDeriv_vec_[FIRST_DERIV]*sol_timeStepCoeff_[PREDICTOR_1];
    } else {
      omitFirstPredictor_ = false;
    }
  }


  void Trapezoidal::UpdateRHS()
  {

    Vector<Double> coeffMass;

    // mass part
    coeffMass = solpred_*sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS,coeffMass);
  }


  void Trapezoidal::UpdateRHS(Vector<Double>& actSol)
  {

    Vector<Double> coeffMass;

    // mass part
    coeffMass = (solpred_ - actSol) *sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS,coeffMass);
  }


  void Trapezoidal::Corrector(Vector<Double>& solnew)
  {

    solDeriv_vec_[FIRST_DERIV] = (solnew - solpred_)*sol_timeStepCoeff_[CORRECTOR_1];
  }

  void Trapezoidal::CalcParameters(Double dt)
  {

    //for predictors
    sol_timeStepCoeff_[PREDICTOR_1] = (1-gamma_)*dt_;

    //for correctors, matrices
    sol_timeStepCoeff_[CORRECTOR_1] = 1.0/(gamma_*dt_);
  }



  //====================================================================
  //---------------------------- Effective Mass ------------------
  //====================================================================

  TrapezoidalEffMass::TrapezoidalEffMass( BaseSystem * algebraicsystem,
                                          PtrParamNode systemNode)
    :TimeStepping( algebraicsystem )
  {
    if ( systemNode->Has("timeSteppingParameters") ) {
      PtrParamNode myParam = systemNode->Get("timeSteppingParameters");
      myParam->GetValue("omitInitialSol", omitFirstPredictor_, ParamNode::PASS);
    }

    gamma_ = 0.51;

  }

  TrapezoidalEffMass::~TrapezoidalEffMass()
  {

  }

  void TrapezoidalEffMass::Init( Double dt, UInt rhsSize ) {
    
    dt_ = dt;
    rhsSize_ = rhsSize;
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = sol_timeStepCoeff_[CORRECTOR_1];
    matrix_factors_[MASS] = 1.0;

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0;
    matrix_factors_[DAMPING] = 0.0;

    //get the memory
    sol_.Resize(rhsSize_);
    sol_.Init();

    solDeriv_vec_[FIRST_DERIV].Resize(rhsSize_);
    solDeriv_vec_[FIRST_DERIV].Init();

    solpred_.Resize(rhsSize_); 
    solpred_.Init();

  }


  void TrapezoidalEffMass::Predictor(Vector<Double>& solold)
  {
    if( !omitFirstPredictor_) {
      solpred_ = solold + solDeriv_vec_[FIRST_DERIV]*sol_timeStepCoeff_[PREDICTOR_1];
    } else {
      omitFirstPredictor_ = false;
    }
  }


  void TrapezoidalEffMass::UpdateRHS()
  {

    Vector<Double> coeffStiff;

    // mass part
    coeffStiff = -solpred_;
    algsys_->UpdateRHS(STIFFNESS,coeffStiff);
  }


  void TrapezoidalEffMass::UpdateRHS(Vector<Double>& actSol)
  {
    UpdateRHS();
  }


  void TrapezoidalEffMass::Corrector(Vector<Double>& vNew)
  {
    // after solving the algebraic system of equation, we obtain as solution
    // the 1st time derivative: vNew .. 1st time derivative
    sol_ = solpred_ + vNew * sol_timeStepCoeff_[CORRECTOR_1];
    solDeriv_vec_[FIRST_DERIV] = vNew;

    //now overwrite the solution with the physical quantity itself
    vNew = sol_;
  }

  void TrapezoidalEffMass::CalcParameters(Double dt)
  {

    //for predictors
    sol_timeStepCoeff_[PREDICTOR_1] = (1-gamma_)*dt;

    //for correctors, matrices
    sol_timeStepCoeff_[CORRECTOR_1] = (gamma_*dt);
  }


  Double TrapezoidalEffMass::DirichletBC4EffMassMatrix(Double val, Integer eq)
  {

    Double velVal;

    velVal = (val- solpred_[eq-1]) / sol_timeStepCoeff_[CORRECTOR_1];
    return velVal;
  }




} // end of namespace
