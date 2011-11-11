#include <fstream>
#include <iostream>
#include <string>

#include "Trapezoidal.hh"

#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField
{

  Trapezoidal::Trapezoidal( AlgebraicSys * algebraicsystem, 
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
    }

  }

  Trapezoidal::~Trapezoidal()
  {

  }

  void Trapezoidal::Init( Double dt, UInt rhsSize ) {
    
    dt_ = dt;
    rhsSize_ = rhsSize;
    SBM_Vector dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();
    
    CalcParameters(dt_);

    matrix_factors_[STIFFNESS] = 1.0;
    matrix_factors_[MASS] = sol_timeStepCoeff_[CORRECTOR_1];

    //not used matrices
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 0.0;       

    //get the memory
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


  void Trapezoidal::Predictor(SBM_Vector& solold)
  {
    if( !omitFirstPredictor_) {
      solpred_.Add( 1.0, solold,
                    sol_timeStepCoeff_[PREDICTOR_1], solDeriv_vec_[FIRST_DERIV] );
//      solpred_ = solold + solDeriv_vec_[FIRST_DERIV]*sol_timeStepCoeff_[PREDICTOR_1];
    } else {
      omitFirstPredictor_ = false;
    }
  }


  void Trapezoidal::UpdateRHS()
  {

    SBM_Vector coeffMass;

    // mass part
    coeffMass = solpred_;
    coeffMass.ScalarMult( sol_timeStepCoeff_[CORRECTOR_1] );
    //coeffMass = solpred_*sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS,coeffMass);
  }


  void Trapezoidal::UpdateRHS(SBM_Vector& actSol)
  {

    SBM_Vector coeffMass;

    // mass part
    coeffMass.Add( sol_timeStepCoeff_[CORRECTOR_1], solpred_,
                   -sol_timeStepCoeff_[CORRECTOR_1], actSol );
    //coeffMass = (solpred_ - actSol) *sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS,coeffMass);
  }


  void Trapezoidal::Corrector(SBM_Vector& solnew)
  {

    solDeriv_vec_[FIRST_DERIV].Add( sol_timeStepCoeff_[CORRECTOR_1], solnew,
                                    - sol_timeStepCoeff_[CORRECTOR_1], solpred_ );
    //solDeriv_vec_[FIRST_DERIV] = (solnew - solpred_)*sol_timeStepCoeff_[CORRECTOR_1];
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

  TrapezoidalEffMass::TrapezoidalEffMass( AlgebraicSys * algebraicsystem,
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


  void TrapezoidalEffMass::Predictor(SBM_Vector& solold)
  {
    if( !omitFirstPredictor_) {
      solpred_.Add( 1.0, solold, 
                    sol_timeStepCoeff_[PREDICTOR_1], solDeriv_vec_[FIRST_DERIV] );
         
      //solpred_ = solold + solDeriv_vec_[FIRST_DERIV]*sol_timeStepCoeff_[PREDICTOR_1];
    } else {
      omitFirstPredictor_ = false;
    }
  }


  void TrapezoidalEffMass::UpdateRHS()
  {

    SBM_Vector coeffStiff;

    // mass part
    coeffStiff = solpred_;
    coeffStiff.ScalarMult( -1.0 );
//    coeffStiff = -solpred_;
    algsys_->UpdateRHS(STIFFNESS,coeffStiff);
  }


  void TrapezoidalEffMass::UpdateRHS(SBM_Vector& actSol)
  {
    UpdateRHS();
  }


  void TrapezoidalEffMass::Corrector(SBM_Vector& vNew)
  {
    // after solving the algebraic system of equation, we obtain as solution
    // the 1st time derivative: vNew .. 1st time derivative
    sol_.Add( 1.0, solpred_,
              sol_timeStepCoeff_[CORRECTOR_1], vNew ); 
    
//    sol_ = solpred_ + vNew * sol_timeStepCoeff_[CORRECTOR_1];
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
    EXCEPTION("This method won't work correctly. We should consider "
                 "modifying the Dirichlet BC handling in the algebraic system "
                 "to associate the IDBC-matrix with the STIFFNESS-matrix, so "
                 "that is gets automatically multiplicated within the "
                 "AlgebraicSys::ConstructEffectiveMatrix() method" );
//    Double velVal;
//
//    velVal = (val- solpred_[eq-1]) / sol_timeStepCoeff_[CORRECTOR_1];
//    return velVal;
    return 0;
  }




} // end of namespace
