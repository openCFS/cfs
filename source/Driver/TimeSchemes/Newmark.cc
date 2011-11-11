#include <fstream>
#include <iostream>
#include <string>

#include "Newmark.hh"


#include "OLAS/algsys/AlgebraicSys.hh"
//#include "PDE/StdPDE.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"


namespace CoupledField
{

  Newmark::Newmark(AlgebraicSys * algebraicsystem, PtrParamNode systemNode )
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
      myParam->GetValue("omitInitialSol", omitFirstPredictor_, ParamNode::PASS);
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
    SBM_Vector dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();

    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
    globalDeltaT_=dt;

    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*sol_timeStepCoeff_[CORRECTOR_1];   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*sol_timeStepCoeff_[COEFFRHS];

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
    solderiv1pred_.Resize(rhsSize_);
    solderiv1pred_.Init();

  }
  
  void Newmark::ReInit(){
    std::map<TIMEStepType, SBM_Vector >::iterator iterTime = sol_timeStepVec_.begin();
    std::map<DERIVType, SBM_Vector >::iterator iterDeriv = solDeriv_vec_.begin();
    for (;iterTime != sol_timeStepVec_.end(); ++iterTime)
    {
      iterTime->second.Init();
    }
    for (;iterDeriv != solDeriv_vec_.end(); ++iterDeriv)
    {
      iterDeriv->second.Init();
    }

    solpred_.Init();
    solderiv1pred_.Init();
  }

  void Newmark::setSubSteps(  UInt subSteps ) {

    dt_ = globalDeltaT_/Double(subSteps);
    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*sol_timeStepCoeff_[CORRECTOR_1];   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*sol_timeStepCoeff_[COEFFRHS];

    solpredSAVE_ = solpred_;
    solderiv1predSAVE_ = solderiv1pred_;
  }

  void Newmark::resetDeltaT( ) {

    dt_ = globalDeltaT_;
    CalcParameters(dt_);
        
    matrix_factors_[STIFFNESS] = (1.0 + alpha_);
    matrix_factors_[MASS] = 1.0*sol_timeStepCoeff_[CORRECTOR_1];   
    matrix_factors_[CONVECTION] = 0.0; 
    matrix_factors_[DAMPING] = 1.0*sol_timeStepCoeff_[COEFFRHS];

    solpred_ = solpredSAVE_;
    solderiv1pred_ = solderiv1predSAVE_;
  }

  void Newmark::Predictor(SBM_Vector& solold)
  {
    if( !omitFirstPredictor_) {

      solpred_ = solold;
      solpred_.Add( dt_, solDeriv_vec_[FIRST_DERIV] );
      solpred_.Add( sol_timeStepCoeff_[PREDICTOR_1], 
                    solDeriv_vec_[SECOND_DERIV] );
      
      solderiv1pred_ = solDeriv_vec_[FIRST_DERIV];
      solderiv1pred_.Add( sol_timeStepCoeff_[PREDICTOR_2],
                          solDeriv_vec_[SECOND_DERIV] );
      // old, Single-Vector based code
//      solpred_ = solold + solDeriv_vec_[FIRST_DERIV] * dt_ \
//               + solDeriv_vec_[SECOND_DERIV] * sol_timeStepCoeff_[PREDICTOR_1];
//      solderiv1pred_ = solDeriv_vec_[FIRST_DERIV] 
//                    + solDeriv_vec_[SECOND_DERIV] * sol_timeStepCoeff_[PREDICTOR_2];
    } else {
      omitFirstPredictor_ = false;
    }
  }

  void Newmark::UpdateRHS()
  {

    SBM_Vector coeffMass;

    // mass part
    coeffMass = solpred_;
    coeffMass.ScalarMult( sol_timeStepCoeff_[CORRECTOR_1] );
    
//    coeffMass = solpred_*sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS,coeffMass);

    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      SBM_Vector coeffDamp;
      coeffDamp.Add( -1.0, solderiv1pred_,
                     sol_timeStepCoeff_[COEFFRHS], solpred_ );
      
//      coeffDamp = -solderiv1pred_ + solpred_*sol_timeStepCoeff_[COEFFRHS];
      algsys_->UpdateRHS(DAMPING,coeffDamp);
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      SBM_Vector coeffStiff;
      coeffStiff = solpred_;
      coeffStiff.ScalarMult(alpha_);
//      coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff);
    }

  }

  void Newmark::SubstractStiffnessFromRHS(SBM_Vector& actSol) {


    SBM_Vector actSolTemp;
    actSolTemp = actSol;
    actSolTemp.ScalarMult(-1.0);
//    actSolTemp=-actSol;

    algsys_->UpdateRHS(STIFFNESS,actSolTemp);
  }

  void Newmark::UpdateRHS(SBM_Vector& actSol)
  {

    // mass part
    SBM_Vector coeffMass;
    coeffMass.Add( sol_timeStepCoeff_[CORRECTOR_1], solpred_,
                  -sol_timeStepCoeff_[CORRECTOR_1], actSol ); 
    
    //coeffMass = (solpred_ - actSol) * sol_timeStepCoeff_[CORRECTOR_1];
    algsys_->UpdateRHS(MASS, coeffMass);


    // damping part
    if ( FeMatrixPresent( DAMPING ) ) {
      SBM_Vector coeffDamp;
      coeffDamp.Add( sol_timeStepCoeff_[COEFFRHS], solpred_,
                     - sol_timeStepCoeff_[COEFFRHS], actSol );
      coeffDamp.Add( -1.0, solderiv1pred_);
      //coeffDamp = -solderiv1pred_ + (solpred_-actSol)*sol_timeStepCoeff_[COEFFRHS];
      algsys_->UpdateRHS(DAMPING,coeffDamp);
    }

    //just in case of alpha-Method
    if (abs(alpha_) > 0) {
      SBM_Vector coeffStiff;
      coeffStiff = solpred_;
      coeffStiff.ScalarMult( alpha_ );
      //coeffStiff = solpred_*alpha_;
      algsys_->UpdateRHS(DAMPING,coeffStiff);
    }

  }

  void Newmark::Corrector(SBM_Vector& solnew)
  {
    solDeriv_vec_[SECOND_DERIV].Add( sol_timeStepCoeff_[CORRECTOR_1], solnew,
                                     -sol_timeStepCoeff_[CORRECTOR_1], solpred_ );
    solDeriv_vec_[FIRST_DERIV].Add( 1.0, solderiv1pred_,
                                    sol_timeStepCoeff_[CORRECTOR_2], 
                                    solDeriv_vec_[SECOND_DERIV] );
//    solDeriv_vec_[SECOND_DERIV] = (solnew - solpred_) * sol_timeStepCoeff_[CORRECTOR_1];
//    solDeriv_vec_[FIRST_DERIV] = solderiv1pred_ +
//        solDeriv_vec_[SECOND_DERIV]*sol_timeStepCoeff_[CORRECTOR_2];
  }



  void Newmark::CalcParameters(Double dt)
  {

    //for predictors
    sol_timeStepCoeff_[PREDICTOR_1] = 0.5*(1-2.0*beta_)*dt_*dt_;
    sol_timeStepCoeff_[PREDICTOR_2] = (1-gamma_)*dt_;

    //for correctors, matrices
    sol_timeStepCoeff_[CORRECTOR_1] = 1.0/(beta_*dt_*dt_);
    sol_timeStepCoeff_[CORRECTOR_2] = gamma_*dt_;

    //for RHS, Matrices
    sol_timeStepCoeff_[COEFFRHS] = (1+alpha_) * gamma_ / (beta_*dt_);
  }

  void Newmark::AdvanceTimestep(SBM_Vector& solnew)
  {
    sol_timeStepVec_[TIMESTEP_1]=solnew;
  }


  // ====================================================
  // Effective Mass Matrix Formulation
  // ====================================================

  NewmarkEffMass::NewmarkEffMass(AlgebraicSys * algebraicsystem, PtrParamNode systemNode,
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
      myParam->GetValue("omitInitialSol", omitFirstPredictor_, ParamNode::PASS);
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
    SBM_Vector dummyVec;
    dummyVec.Resize(rhsSize_);
    dummyVec.Init();

    
    // Calculate parameters and store it in matrix_factors
    dt_ = dt;
    CalcParameters(dt_);

    if ( intExplicit_ == true ) {
      matrix_factors_[MASS] = 1.0; 
    }
    else {
      matrix_factors_[STIFFNESS] = 1.0*sol_timeStepCoeff_[CORRECTOR_1];
      matrix_factors_[MASS] = 1.0;
      matrix_factors_[CONVECTION] = 0.0;
      
      matrix_factors_[DAMPING] = 1.0*sol_timeStepCoeff_[CORRECTOR_2]; 
    }

    //get the memory
    sol_.Resize(rhsSize_);
    sol_.Init();
    solpred_.Resize(rhsSize_);
    solpred_.Init();

    if ( !is_Deriv_set(FIRST_DERIV) )
    {
      solDeriv_vec_[FIRST_DERIV] = dummyVec;
    }
    if ( !is_Deriv_set(SECOND_DERIV) )
    {
      solDeriv_vec_[SECOND_DERIV] = dummyVec;
    }

  
    solderiv1pred_.Resize(rhsSize_);
    solderiv1pred_.Init();

  }


  void NewmarkEffMass::Predictor(SBM_Vector& solold)
  {

    if( !omitFirstPredictor_) {
      solpred_ = sol_;
      solpred_.Add( dt_, solDeriv_vec_[FIRST_DERIV] );
      solpred_.Add( sol_timeStepCoeff_[PREDICTOR_1], 
                    solDeriv_vec_[SECOND_DERIV] );
      
      solderiv1pred_ = solDeriv_vec_[FIRST_DERIV];
      solderiv1pred_.Add( sol_timeStepCoeff_[PREDICTOR_2],
                          solDeriv_vec_[SECOND_DERIV] );
//      solpred_ = sol_ + solDeriv_vec_[FIRST_DERIV] * dt_ \
//        + solDeriv_vec_[SECOND_DERIV]*sol_timeStepCoeff_[PREDICTOR_1];
//      solderiv1pred_ = solDeriv_vec_[FIRST_DERIV] \
//        + solDeriv_vec_[SECOND_DERIV]*sol_timeStepCoeff_[PREDICTOR_2];
    } else {
      omitFirstPredictor_ = false;
    }


  }

  void NewmarkEffMass::Corrector(SBM_Vector& aNew)
  {

    // after solving the algebraic system of equation, we obtain as solution
    // the 2nd time derivative: aNew .. 2nd second time derivative
    sol_.Add( 1.0, solpred_, 
              sol_timeStepCoeff_[CORRECTOR_1], aNew );
    solDeriv_vec_[FIRST_DERIV].Add( 1.0, solderiv1pred_,
                                    sol_timeStepCoeff_[CORRECTOR_2], aNew );
//    sol_ = solpred_ + aNew * sol_timeStepCoeff_[CORRECTOR_1];
//    solDeriv_vec_[FIRST_DERIV] = solderiv1pred_ + aNew * sol_timeStepCoeff_[CORRECTOR_2];
    solDeriv_vec_[SECOND_DERIV] = aNew;

    //now overwrite the solution with the physical quantity itself
    aNew = sol_;

  }

  void NewmarkEffMass::UpdateRHS()
  {

    SBM_Vector coeffStiff;
    coeffStiff = solpred_;
    coeffStiff.ScalarMult( -1.0 );
//    coeffStiff = -solpred_;

    algsys_->UpdateRHS(STIFFNESS, coeffStiff);

    // damping part
    if ( FeMatrixPresent( DAMPING ) )
      {
        SBM_Vector coeffDamp;
        coeffDamp = solderiv1pred_;
        coeffDamp.ScalarMult(-1.0);
//        coeffDamp = -solderiv1pred_;

        algsys_->UpdateRHS(DAMPING, coeffDamp);
      }
  }

  void NewmarkEffMass::CalcParameters(Double dt)
  {

    //for predictors
    sol_timeStepCoeff_[PREDICTOR_1] = 0.5*dt_*dt_*(1-2*beta_);
    sol_timeStepCoeff_[PREDICTOR_2] = (1-gamma_)*dt_;
  
    //for correctors, matrices
    sol_timeStepCoeff_[CORRECTOR_1] = beta_*dt_*dt_;
    sol_timeStepCoeff_[CORRECTOR_2] = gamma_*dt_;
  
  }

  Double NewmarkEffMass::DirichletBC4EffMassMatrix(Double val, Integer eq)
  {
    EXCEPTION("This method won't work correctly. We should consider "
              "modifying the Dirichlet BC handling in the algebraic system "
              "to associate the IDBC-matrix with the STIFFNESS-matrix, so "
              "that is gets automatically multiplicated within the "
              "AlgebraicSys::ConstructEffectiveMatrix() method" );
//    Double accval;
//  
//    accval = (val- solpred_[eq-1]) / sol_timeStepCoeff_[CORRECTOR_1];
//    return accval;
    return 0.0;
  }



} // end of namespace
