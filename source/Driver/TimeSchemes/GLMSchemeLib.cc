// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     GLMSchemeLib.cc 
 *       \brief    Implementation file for the GLM Schemes mainly tableau definitions
 *
 *       \date     Jan 23, 2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#include "GLMSchemeLib.hh"
#include "Utils/mathParser/mathParser.hh"


namespace CoupledField{

//================================================================
//GENERAL SCHEME
//================================================================

GLMScheme::GLMScheme(){
  curTStepSize_ = 0;
}

GLMScheme::~GLMScheme(){

}

Double GLMScheme::TransformBC(const StdVector< SingleVector* > & glm,
                                 Double value,
                                 UInt valDerivOrder,
                                 Integer eqnNumber){
  //this algorithm is correct for trapezoidal, newmark and alpha method
  if(valDerivOrder == solDerivOrder_){
    return value;
  }else{
    Double retVal = 0.0;
    //temporarily transform the coefMatrix to valDerivOrder
    UInt tmpSolDeriv = solDerivOrder_;
    ComputeCoefficients(valDerivOrder,curTStepSize_);
    //we just take the corresponding line of the V-Matrix and multiply the glm entiries
    UInt row = numStages_ * (maxDerivOrder_+1) + tmpSolDeriv;
    UInt col = numStages_;
    for(UInt i=0;i<sizeGLMVec_;i++){
      Double tmpVal = 0.0;
      glm[i]->GetEntry(eqnNumber-1,tmpVal);
      retVal += tmpVal * schemeCoefs_[row][col+i];
    }
    //now we divide value by the special entry and we are good to go
    retVal += value * schemeCoefs_[row][col-1];
    //restore the coefficient matrix
    ComputeCoefficients(tmpSolDeriv,curTStepSize_);
    return retVal;
  }
}



//================================================================
//TRAPEZOIDAL SCHEME
//================================================================

Trapezoidal::Trapezoidal(Double gamma)
              : GLMScheme(){
  gamma_ = gamma;

  maxDerivOrder_ = 1;
  numStages_ = 1;
  numOldSols_ = 1;
  numSol1stDerivs_ = 1;
  numSol2ndDerivs_ = 0;
  sizeGLMVec_ = numOldSols_ + numSol1stDerivs_;

  lastStageIsSolution_ = true;
  usePredictors_ = true;

  //prepare coefficient matrix
  UInt numCols = numStages_ + ((maxDerivOrder_+1) * numOldSols_);
  UInt numRows = (maxDerivOrder_+1)*(numStages_) + sizeGLMVec_;
  schemeCoefs_.Resize(numRows,numCols);
  schemeCoefs_.Init();
}

void Trapezoidal::ComputeCoefficients(UInt solDerivOrder,Double deltaT){
  curTStepSize_ = deltaT;
  solDerivOrder_ = solDerivOrder;

  switch(solDerivOrder){
  case 1:
    solDerivOrder_ = 1;
    schemeCoefs_[0][0] = gamma_ * curTStepSize_;        //a_11
    schemeCoefs_[0][1] = -1;                            //U_11
    schemeCoefs_[0][2] = (- 1.0 + gamma_) * curTStepSize_;  //U_12

    schemeCoefs_[1][0] = 1;                             //a_21
    schemeCoefs_[1][1] = 0;                             //U_21
    schemeCoefs_[1][2] = 0;                             //U_22

    schemeCoefs_[2][0] = gamma_ * curTStepSize_;        //b_11
    schemeCoefs_[2][1] = 1;                             //V_11
    schemeCoefs_[2][2] = (1.0-gamma_) * curTStepSize_;  //V_12

    schemeCoefs_[3][0] = 1.0;                           //b_21
    schemeCoefs_[3][1] = 0;                             //V_21
    schemeCoefs_[3][2] = 0;                             //V_22
    break;
  case 0:
    solDerivOrder_ = 0;
    schemeCoefs_[0][0] = 1;
    schemeCoefs_[0][1] = 0;
    schemeCoefs_[0][2] = 0;
    schemeCoefs_[1][0] =    1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[1][1] =   1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[1][2] =  (1.0-gamma_) / gamma_;
    schemeCoefs_[2][0] = 1;
    schemeCoefs_[2][1] = 0;
    schemeCoefs_[2][2] = 0;
    schemeCoefs_[3][0] =    1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[3][1] =  - 1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[3][2] =  - (1.0-gamma_) / gamma_;
    break;
  }
}

//================================================================
//NEWMARK SCHEME
//================================================================

Newmark::Newmark(Double gamma,Double beta, Double alpha)
         : GLMScheme() {
  gamma_ = gamma;
  beta_ = beta;
  alpha_ = alpha;

  maxDerivOrder_ = 2;
  numStages_ = 1;
  numOldSols_ = 1;
  numSol1stDerivs_ = 1;
  numSol2ndDerivs_ = 1;
  sizeGLMVec_ = numOldSols_ + numSol1stDerivs_ + numSol2ndDerivs_;


  if(alpha == 0.0){
    usePredictors_ = true;
    lastStageIsSolution_ = true;
  }else {
    //alpha method is only implemented for effective stiffness right now
    usePredictors_ = false;
    lastStageIsSolution_ = false;
    //redefine beta and gamma accorin=ding to hughes
    beta_ = (1-alpha_)*(1-alpha_)/4;
    gamma_ = (1-2*alpha_)/2;
  }

  //prepare coefficient matrix
  UInt numCols = numStages_ + ((maxDerivOrder_+1) * numOldSols_);
  UInt numRows = (maxDerivOrder_+1)*(numStages_) + sizeGLMVec_;
  //THis is the effective mass tableau
  schemeCoefs_.Resize(numRows,numCols);
  schemeCoefs_.Init();
}

void Newmark::ComputeCoefficients(UInt solDerivOrder,Double deltaT){
  curTStepSize_ = deltaT;
  solDerivOrder_ = solDerivOrder;

  if(solDerivOrder != 0 && alpha_ != 0.0){
    EXCEPTION("Alpha-Method timestepping currently only supports effective stiffness formulation");
  }

  switch(solDerivOrder){
  case 0:
    solDerivOrder_ = 0;
    //zero order part
    schemeCoefs_[0][0] = 1.0+alpha_;
    schemeCoefs_[0][1] = alpha_;
    schemeCoefs_[0][2] = 0.0;
    schemeCoefs_[0][3] = 0.0;
    //fist order part
    schemeCoefs_[1][0] =  (1+alpha_) * gamma_ / (beta_ * curTStepSize_);
    schemeCoefs_[1][1] =  (1+alpha_) * gamma_ / (beta_ * curTStepSize_);
    schemeCoefs_[1][2] = ((1+alpha_) * gamma_ / beta_) - 1.0;
    schemeCoefs_[1][3] =  (1+alpha_) * curTStepSize_ * ((gamma_/beta_) - 2.0) * 0.5;
    //second order part
    schemeCoefs_[2][0] = 1.0 / ( beta_ * curTStepSize_ * curTStepSize_);
    schemeCoefs_[2][1] = 1.0 / ( beta_ * curTStepSize_ * curTStepSize_);
    schemeCoefs_[2][2] = 1.0 / (beta_ * curTStepSize_);
    schemeCoefs_[2][3] = (0.5/beta_ - 1) ;

    //UPDATE PART zero order
    schemeCoefs_[3][0] = 1.0;
    schemeCoefs_[3][1] = 0.0;
    schemeCoefs_[3][2] = 0.0;
    schemeCoefs_[3][3] = 0.0;
    //UPDATE PART first order
    schemeCoefs_[4][0] = 1.0 * gamma_ / (beta_ * curTStepSize_);
    schemeCoefs_[4][1] = -1.0 * gamma_ / (beta_ * curTStepSize_);
    schemeCoefs_[4][2] = (-1.0 * gamma_ / beta_) + 1.0;
    schemeCoefs_[4][3] = curTStepSize_ * (2.0 - (gamma_ / beta_)) * 0.5;
    //UPDATE PART second order
    schemeCoefs_[5][0] = 1.0 / (beta_ * curTStepSize_ * curTStepSize_);
    schemeCoefs_[5][1] = -1.0 / (beta_ * curTStepSize_ * curTStepSize_);
    schemeCoefs_[5][2] = -1.0 / (beta_ * curTStepSize_);
    schemeCoefs_[5][3] = (1.0 - 0.5/beta_);
    break;
  case 1 :
    solDerivOrder_ = 1;
    //zero order part
    schemeCoefs_[0][0] = curTStepSize_ * beta_ / gamma_;
    schemeCoefs_[0][1] = -1.0;
    schemeCoefs_[0][2] = curTStepSize_ * ( (beta_ / gamma_) - 1.0);
    schemeCoefs_[0][3] = ((beta_ / gamma_) - 0.5) * curTStepSize_ * curTStepSize_;
    //fist order part
    schemeCoefs_[1][0] = 1.0;
    schemeCoefs_[1][1] = 0.0;
    schemeCoefs_[1][2] = 0.0;
    schemeCoefs_[1][3] = 0.0;
    //second order part
    schemeCoefs_[2][0] = 1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[2][1] = 0.0;
    schemeCoefs_[2][2] = 1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[2][3] = (1.0 / gamma_) - 1.0;

    //UPDATE PART zero order
    schemeCoefs_[3][0] = curTStepSize_ * beta_ / gamma_;
    schemeCoefs_[3][1] = 1.0;
    schemeCoefs_[3][2] = curTStepSize_ * (1.0 - (beta_ / gamma_));
    schemeCoefs_[3][3] = (0.5 - (beta_ / gamma_)) * curTStepSize_ * curTStepSize_;
    //UPDATE PART first order
    schemeCoefs_[4][0] = 1.0;
    schemeCoefs_[4][1] = 0.0;
    schemeCoefs_[4][2] = 0.0;
    schemeCoefs_[4][3] = 0.0;
    //UPDATE PART second order
    schemeCoefs_[5][0] = 1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[5][1] = 0.0;
    schemeCoefs_[5][2] = -1.0 / (gamma_ * curTStepSize_);
    schemeCoefs_[5][3] = 1.0 - (1.0 / gamma_);
    break;
  case 2 :
    solDerivOrder_ = 2;
    //zero order part
    schemeCoefs_[0][0] = beta_ * curTStepSize_ * curTStepSize_;
    schemeCoefs_[0][1] = -1.0;
    schemeCoefs_[0][2] = -1.0 * curTStepSize_;
    schemeCoefs_[0][3] = (beta_ - 0.5) * curTStepSize_ * curTStepSize_;
    //fist order part
    schemeCoefs_[1][0] = gamma_ * curTStepSize_;
    schemeCoefs_[1][1] = 0.0;
    schemeCoefs_[1][2] = -1.0;
    schemeCoefs_[1][3] = (gamma_-1.0) * curTStepSize_;
    //second order part
    schemeCoefs_[2][0] = 1.0;
    schemeCoefs_[2][1] = 0.0;
    schemeCoefs_[2][2] = 0.0;
    schemeCoefs_[2][3] = 0.0;

    //zero order part
    schemeCoefs_[3][0] = beta_ * curTStepSize_ * curTStepSize_;
    schemeCoefs_[3][1] = 1.0;
    schemeCoefs_[3][2] = 1.0 * curTStepSize_;
    schemeCoefs_[3][3] = (0.5 - beta_) * (1+alpha_) * curTStepSize_ * curTStepSize_;
    //fist order part
    schemeCoefs_[4][0] = gamma_ * curTStepSize_;
    schemeCoefs_[4][1] = 0.0;
    schemeCoefs_[4][2] = 1.0;
    schemeCoefs_[4][3] = (1.0 - gamma_) * curTStepSize_;
    //second order part
    schemeCoefs_[5][0] = 1.0;
    schemeCoefs_[5][1] = 0.0;
    schemeCoefs_[5][2] = 0.0;
    schemeCoefs_[5][3] = 0.0;
    break;

  }
}

void Newmark::PrepareStage(UInt i,Double aTime, Domain* domain){
 domain->GetMathParser()->SetValue( MathParser_GLOB_HANDLER,
                                    "t", aTime+(alpha_*curTStepSize_) );
}


//================================================================
//BDF2 SCHEME
//================================================================

Bdf2::Bdf2()
      : GLMScheme(){

  maxDerivOrder_ = 1;
  numStages_ = 1;
  numOldSols_ = 2;
  numSol1stDerivs_ = 1;
  numSol2ndDerivs_ = 0;
  sizeGLMVec_ = numOldSols_ + numSol1stDerivs_;
  lastStageIsSolution_ = false;
  usePredictors_ = false;

  //prepare coefficient matrix
  UInt numCols = numStages_ + ((maxDerivOrder_+1) * numOldSols_);
  UInt numRows = (maxDerivOrder_+1)*(numStages_) + sizeGLMVec_;
  schemeCoefs_.Resize(numRows,numCols);
  schemeCoefs_.Init();
}

void Bdf2::ComputeCoefficients(UInt solDerivOrder,Double deltaT){  
  // --------------------------------------------------------------------------------
  // Adaptive Timestepping:
  // Stores Previus Timsteps for Errorestimation
  // The Scheme has been updated so, constant aswell as adaptive timesteps can be used.
  // w_ -> 1 for constant timstep, the coefficiats then revert back to previus values.
  // --------------------------------------------------------------------------------
  if(!initialized_)
  {
    dtCurrent_ = deltaT;
    dtPrev1_ = deltaT;
    dtPrev2_ = deltaT;
    initialized_ = true;
    coefChanged_ = true;
  }else
  {
    Double oldDt = dtCurrent_;
    dtPrev2_ = dtPrev1_;
    dtPrev1_ = dtCurrent_;
    dtCurrent_ = deltaT;
    coefChanged_ = (dtCurrent_ != oldDt);
  }

  curTStepSize_ = deltaT;
  solDerivOrder_ = solDerivOrder;
   w_ = dtPrev1_/dtCurrent_; // For constant Timestep w_=1 -> reverts to non adaptive BDF2
  double a0 = (1.0 + 2.0*w_) / (1.0 + w_);

  switch(solDerivOrder){
  case 1:
    solDerivOrder_ = 1;
    schemeCoefs_[0][0] = dtCurrent_ / a0;          // 2*curTStepSize_/3;
    schemeCoefs_[0][1] = -(1.0 + w_) / a0;         // -4/3;
    schemeCoefs_[0][2] = w_*w_ / ((1.0 + 2.0*w_));  // 1/3;
    schemeCoefs_[0][3] = 0;
    schemeCoefs_[1][0] = 1;
    schemeCoefs_[1][1] = 0;
    schemeCoefs_[1][2] = 0;
    schemeCoefs_[1][3] = 0;
    schemeCoefs_[2][0] = dtCurrent_ / a0;             // 2*curTStepSize_/3;
    schemeCoefs_[2][1] = (1.0 + w_) / a0;            // 4/3;
    schemeCoefs_[2][2] = -w_*w_ / (1.0 + 2.0*w_);    // -1/3;
    schemeCoefs_[2][3] = 0;
    schemeCoefs_[3][0] = 0;
    schemeCoefs_[3][1] = 0;
    schemeCoefs_[3][2] = 1;
    schemeCoefs_[3][3] = 0;
    schemeCoefs_[4][0] = 1;
    schemeCoefs_[4][1] = 0;
    schemeCoefs_[4][2] = 0;
    schemeCoefs_[4][3] = 0;
    break;
  case 0:
    solDerivOrder_ = 0;
    schemeCoefs_[0][0] = 1;
    schemeCoefs_[0][1] = 0;
    schemeCoefs_[0][2] = 0;
    schemeCoefs_[0][3] = 0;
    schemeCoefs_[1][0] = (1.0 + 2.0 *w_)/ ((1.0 + w_)*dtCurrent_); // 3/(2*curTStepSize_);
    schemeCoefs_[1][1] =  (1.0 +w_) / dtCurrent_;                 //  2/(curTStepSize_);
    schemeCoefs_[1][2] = -w_*w_ /((1.0+w_) * dtCurrent_);         // -0.5/curTStepSize_;
    schemeCoefs_[1][3] = 0;
    schemeCoefs_[2][0] = 1;
    schemeCoefs_[2][1] = 0;
    schemeCoefs_[2][2] = 0;
    schemeCoefs_[2][3] = 0;
    schemeCoefs_[3][0] = 0;
    schemeCoefs_[3][1] = 1;
    schemeCoefs_[3][2] = 0;
    schemeCoefs_[3][3] = 0;
    schemeCoefs_[4][0] = (1.0 + 2.0 * w_) / ((1.0 +w_) * dtCurrent_); // 3/(2*curTStepSize_);
    schemeCoefs_[4][1] = -(1.0 + w_) / dtCurrent_;                     // -2/(curTStepSize_);
    schemeCoefs_[4][2] = w_*w_ / ((1.0 + w_)* dtCurrent_);            // (0.5/curTStepSize_);
    schemeCoefs_[4][3] = 0;
    break;
  }
}


//================================================================
// RUNGE KUTTA fourth ORDER SCHEME
//================================================================

RungeKutta4::RungeKutta4()
         : GLMScheme() {

  maxDerivOrder_ = 1;
  numStages_ = 4;
  numOldSols_ = 1;
  numSol1stDerivs_ = 1;
  numSol2ndDerivs_ = 0;
  sizeGLMVec_ = numOldSols_ + numSol1stDerivs_;

  lastStageIsSolution_ = false;
  usePredictors_ = false;

  //prepare coefficient matrix
  UInt numCols = numStages_ + ((maxDerivOrder_+1) * numOldSols_);
  UInt numRows = (maxDerivOrder_+1)*(numStages_) + sizeGLMVec_;
  //THis is the effective mass tableau
  schemeCoefs_.Resize(numRows,numCols);
  schemeCoefs_.Init();
}

Double RungeKutta4::TransformBC(const StdVector< SingleVector* > & glm,
                                 Double value,
                                 UInt valDerivOrder,
                                 Integer eqnNumber){
  //test if the value deriv order is 0 if so, we compute an approximated time derivative
  if(valDerivOrder == solDerivOrder_){
      return value;
    }else if (valDerivOrder == 0){
      std::cerr << "WARNING: in case of Explicit schemes you need to specify the time derivative of the dirichlet value" << std::endl;
      return value;
      //Double retVal = 0.0;
      ////temporarily transform the coefMatrix to valDerivOrder
      //UInt tmpSolDeriv = solDerivOrder_;
      //ComputeCoefficients(valDerivOrder,curTStepSize_);
      ////we just take the corresponding line of the V-Matrix and multiply the glm entiries
      //UInt row = numStages_ * (maxDerivOrder_+1) + tmpSolDeriv;
      //UInt col = numStages_;
      //for(UInt i=0;i<sizeGLMVec_;i++){
      //  Double tmpVal = 0.0;
      //  glm[i]->GetEntry(eqnNumber-1,tmpVal);
      //  retVal += tmpVal * schemeCoefs_[row][col+i];
      //}
      ////now we multiply value by the special entry and we are good to go
      //retVal += value * schemeCoefs_[row][col-1];
      ////restore the coefficient matrix
      //ComputeCoefficients(tmpSolDeriv,curTStepSize_);
      //return retVal;
    }
  return -1;
}

void RungeKutta4::ComputeCoefficients(UInt solDerivOrder,Double deltaT){
  curTStepSize_ = deltaT;
  solDerivOrder_ = solDerivOrder;

  switch(solDerivOrder){
  case 0:
    EXCEPTION("RK4 scheme cannot solve for solution itself")
    break;
  case 1 :
    solDerivOrder_ = 1;

    //FIRST STAGE
    //zero order part
    schemeCoefs_[0][0] = 0.0;
    schemeCoefs_[0][1] = 0.0;
    schemeCoefs_[0][2] = 0.0;
    schemeCoefs_[0][3] = 0.0;
    schemeCoefs_[0][4] = -1.0;
    schemeCoefs_[0][5] = 0.0;
    //fist order part
    schemeCoefs_[1][0] = 1.0;
    schemeCoefs_[1][1] = 0.0;
    schemeCoefs_[1][2] = 0.0;
    schemeCoefs_[1][3] = 0.0;
    schemeCoefs_[1][4] = 0.0;
    schemeCoefs_[1][5] = 0.0;

    //SECOND STAGE
    //zero order part
    schemeCoefs_[2][0] = -curTStepSize_/2.0;
    schemeCoefs_[2][1] = 0.0;
    schemeCoefs_[2][2] = 0.0;
    schemeCoefs_[2][3] = 0.0;
    schemeCoefs_[2][4] = -1.0;
    schemeCoefs_[2][5] = 0.0;
    //fist order part
    schemeCoefs_[3][0] = 0.0;
    schemeCoefs_[3][1] = 1.0;
    schemeCoefs_[3][2] = 0.0;
    schemeCoefs_[3][3] = 0.0;
    schemeCoefs_[3][4] = 0.0;
    schemeCoefs_[3][5] = 0.0;

    //THIRD STAGE
    //zero order part
    schemeCoefs_[4][0] = 0.0;
    schemeCoefs_[4][1] = -curTStepSize_/2.0;
    schemeCoefs_[4][2] = 0.0;
    schemeCoefs_[4][3] = 0.0;
    schemeCoefs_[4][4] = -1.0;
    schemeCoefs_[4][5] = 0.0;
    //fist order part
    schemeCoefs_[5][0] = 0.0;
    schemeCoefs_[5][1] = 0.0;
    schemeCoefs_[5][2] = 1.0;
    schemeCoefs_[5][3] = 0.0;
    schemeCoefs_[5][4] = 0.0;
    schemeCoefs_[5][5] = 0.0;

    //FOURTH STAGE
    // zero order part
    schemeCoefs_[6][0] = 0.0;
    schemeCoefs_[6][1] = 0.0;
    schemeCoefs_[6][2] = -1.0*curTStepSize_;
    schemeCoefs_[6][3] = 0.0;
    schemeCoefs_[6][4] = -1.0;
    schemeCoefs_[6][5] = 0.0;
    // first order part
    schemeCoefs_[7][0] = 0.0;
    schemeCoefs_[7][1] = 0.0;
    schemeCoefs_[7][2] = 0.0;
    schemeCoefs_[7][3] = 1.0;
    schemeCoefs_[7][4] = 0.0;
    schemeCoefs_[7][5] = 0.0;


    //UPDATE PART zero order
    schemeCoefs_[8][0] = (1.0/6.0)*curTStepSize_;
    schemeCoefs_[8][1] = (1.0/3.0)*curTStepSize_;
    schemeCoefs_[8][2] = (1.0/3.0)*curTStepSize_;
    schemeCoefs_[8][3] = (1.0/6.0)*curTStepSize_;
    schemeCoefs_[8][4] = 1.0;
    schemeCoefs_[8][5] = 0.0;

    //UPDATE PART first order
    schemeCoefs_[9][0] = 1.0/6.0;
    schemeCoefs_[9][1] = 1.0/3.0;
    schemeCoefs_[9][2] = 1.0/3.0;
    schemeCoefs_[9][3] = 1.0/6.0;
    schemeCoefs_[9][4] = 0.0;
    schemeCoefs_[9][5] = 0.0;
    break;
  default:
    EXCEPTION("Specified time derivative order not supported by RK4");
    break;

  }
}

void RungeKutta4::PrepareStage(UInt i,Double aTime){
  //obtain current time
  switch(i){
    case 0:
    case 2:
      break;
    case 3:
      //set current time to t+dt
      domain->GetMathParser()->SetValue( MathParser_GLOB_HANDLER,
                                         "t", aTime+curTStepSize_ );
      break;
    case 1:
      //set current time to t+dt/2
      domain->GetMathParser()->SetValue( MathParser_GLOB_HANDLER,
                                         "t", aTime+(0.5*curTStepSize_) );
      break;
    default:
      EXCEPTION("RK4 Called with invalid stage number!");
      break;
   }
};


}
