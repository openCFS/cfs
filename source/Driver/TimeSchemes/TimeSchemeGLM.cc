/*
 * TimeSchemeGLM.cc
 *
 *  Created on: Jan 22, 2012
 *      Author: ahueppe
 */

#include <fstream>

#include "TimeSchemeGLM.hh"
#include "GLMSchemeLib.hh"
#include "MatVec/Vector.hh"

#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField{

  DECLARE_LOG(timeschemeglm)
   DEFINE_LOG(timeschemeglm, "timescheme.glm")


  TimeSchemeGLM::TimeSchemeGLM(SchemeType type, UInt solDerivOrder) :
         avoidUpdateIdx_(-1) {

    InitGLMs();

    curScheme_ = availSchemes[type];
    curType_ = type;

    curScheme_->solDerivOrder_ = solDerivOrder;
    solOrder_ = solDerivOrder;
  }

  TimeSchemeGLM::~TimeSchemeGLM(){
    for(UInt i=1;i<curScheme_->sizeGLMVec_;i++){
       if(avoidFreeingIdx_.find(i)==avoidFreeingIdx_.end())
         delete glmVector_[i];
    }

    for(UInt i=1;i<curScheme_->numStages_;i++){
      if((Integer)i!=avoidUpdateIdx_){
        delete stageVector_[i];
      }
    }

    delete curScheme_;
    curScheme_ = NULL;

    std::map<SchemeType, GLMScheme*>::iterator it = availSchemes.begin();
    while(it != availSchemes.end()){
      delete it->second;
      it->second = NULL;
    }
  }

  void TimeSchemeGLM::Init(SingleVector* solVec,Double dt){

    curScheme_->ComputeCoefficients(curScheme_->solDerivOrder_,dt);

    //now init GLM vector
    //this now highly depends on the used scheme
    //the first component of the glm vector always points to
    //the fe function vector by definition
    glmVector_.Resize(curScheme_->sizeGLMVec_);
    glmVector_[0] = solVec;
    for(UInt i=1;i<curScheme_->sizeGLMVec_;i++){
      glmVector_[i] = new Vector<Double>();
      glmVector_[i]->Resize(solVec->GetSize());
      glmVector_[i]->Init();
    }

    //now we go for the stage vector
    //if the scheme has the lastStageIsSolution_ flag set,
    //the corresponding stage vector points directly to the
    //glm Vector entry and the update step for this entry can be avoided
    stageVector_.Resize(curScheme_->numStages_);
    if(curScheme_->lastStageIsSolution_){
      if(solOrder_==0){
        avoidUpdateIdx_ = 0;
      }else if(solOrder_ == 1){
        avoidUpdateIdx_ = curScheme_->numOldSols_;
      }else if(solOrder_ == 2){
        avoidUpdateIdx_ = curScheme_->numOldSols_ + curScheme_->numSol1stDerivs_ ;
      }
      stageVector_[curScheme_->numStages_-1] = glmVector_[avoidUpdateIdx_];
    }
    for(UInt i=0;i<curScheme_->numStages_-1;i++){
        stageVector_[i] = new Vector<Double>();
        stageVector_[i]->Resize(solVec->GetSize());
        stageVector_[i]->Init();
    }

    //for some schemes we can reuse computed information in predictors
    //we create them anyway either we use them classically as predictors
    //or we use them as a temoprary storage for the update step
    UInt sPred = curScheme_->sizeGLMVec_;
    if(curScheme_->usePredictors_){
      predictors_.Resize(sPred);
      predictorCalculated_.Resize(sPred,false);
      for(UInt i=0;i<sPred;i++){
        predictors_[i] = new Vector<Double>();
        predictors_[i]->Resize(solVec->GetSize());
        predictors_[i]->Init();
      }
    }
  }

  void TimeSchemeGLM::ComputeStageRHS(UInt actStage, Integer derivId,
                                      SingleVector* rhsVec, Integer subIdx){

    //if the derivative id is equal to solution we do not need to do anything
    //as the corresponding line in the GLM is equal to zero
    if(derivId == (Integer)curScheme_->solDerivOrder_ ||
       derivId > (Integer)curScheme_->maxDerivOrder_ ||
       derivId < 0){
      rhsVec->Init();
      return;
    }
    rhsVec->Init();
    UInt dId = (UInt) derivId;
    //Calculate coefficient matrix row index
    UInt cRow = actStage * (curScheme_->maxDerivOrder_+1) + dId;

    if(curScheme_->usePredictors_ && predictorCalculated_[dId]){
      rhsVec->Add((*predictors_[dId]));
    }else{
      //update for old solutions
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        UInt col = curScheme_->numStages_+i;
        Double coef = curScheme_->schemeCoefs_[cRow][col];
        if(coef !=0){
          SingleVector * curVec = glmVector_[i];
          rhsVec->Add(coef,(*curVec));
        }
      }
      if(curScheme_->usePredictors_){
        predictors_[dId]->Add(*rhsVec);
        predictorCalculated_[dId] = true;
      }
    }

    //now loop over each column, scale the GLM vector Entry by the factor and add it to RHS
    for(UInt i=0;i<actStage;i++){
      Double coef = curScheme_->schemeCoefs_[cRow][i];
      if(coef !=0){
        rhsVec->Add(coef,(*stageVector_[i]));
      }
    }
  }

  void TimeSchemeGLM::FinishStep(){
    //update for old solutions
    if(curScheme_->usePredictors_){
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        if((Integer)i != avoidUpdateIdx_){
           glmVector_[i]->Init();
           SingleVector * curPVec = predictors_[i];
           glmVector_[i]->Add(-1.0,*curPVec);
           curPVec = NULL;

           for(UInt actS = 0; actS < curScheme_->numStages_; actS++){
             SingleVector * curSVec = stageVector_[actS];

             UInt row = curScheme_->numStages_ * (curScheme_->maxDerivOrder_+1);
             Double coef = curScheme_->schemeCoefs_[row+i][actS];
             glmVector_[i]->Add(coef,*curSVec);
             curSVec = NULL;
           }
        }
      }
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        predictors_[i]->Init();
        predictorCalculated_[i] = false;
      }
    }else{
      Exception("The general case is not implemented yet");
    }
  }

  void TimeSchemeGLM::AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap
                                          , std::map<FEMatrixType,Double> & matFactors){

    //The correct matrix factors are always at the same position
    std::map<FEMatrixType,Integer>::const_iterator mIt = matMap.begin();
    matFactors.clear();
    for(mIt = matMap.begin(); mIt != matMap.end() ; mIt++){
      if(mIt->second >= 0)
        matFactors[mIt->first] = curScheme_->schemeCoefs_[stage*curScheme_->maxDerivOrder_+ mIt->second][stage];
      else
        matFactors[mIt->first] = 0;
    }
    //
#ifndef NDEBUG
    for(mIt = matMap.begin(); mIt != matMap.end() ; mIt++){
      LOG_DBG(timeschemeglm) << "Matrix " << feMatrixType.ToString(mIt->first) << " value = " << matFactors[mIt->first] << std::endl;
    }
#endif
  }

  SingleVector* TimeSchemeGLM::GetTimeDerivative(UInt order){
    SingleVector* derivVec = NULL;
    if(order > curScheme_->maxDerivOrder_){
      Exception("Request for invalid time derivative order.");
    }else{
      switch(order){
      case 0:
        derivVec = glmVector_[0];
        break;
      case 1:
        if(curScheme_->numSol1stDerivs_ > 0){
          derivVec = glmVector_[curScheme_->numOldSols_];
        }else{
          //TODO: Calculate it, TO BE DONE
        }
        break;
      case 2:
        if(curScheme_->numSol2ndDerivs_ > 0){
          derivVec = glmVector_[curScheme_->numOldSols_+curScheme_->numSol1stDerivs_];
        }else{
          //TODO: Calculate it, TO BE DONE
        }
        break;
      }
    }
    return derivVec;
  }

  void TimeSchemeGLM::SetTimeDerivVector(UInt order,SingleVector * coefVector){
    if(order > curScheme_->maxDerivOrder_){
          Exception("Cannot Set Vector for order higher than scheme deriv order");
    }else{
      switch(order){
        case 1:
          if(curScheme_->numSol1stDerivs_ > 0){
            //delete the current vector and associate the new one
            //we take care if we solve for this solution order...
            delete glmVector_[curScheme_->numOldSols_];
            glmVector_[curScheme_->numOldSols_] = coefVector;
            avoidFreeingIdx_.insert(curScheme_->numOldSols_);
            if(solOrder_ == 1 && avoidUpdateIdx_ == (Integer)curScheme_->numOldSols_){
              stageVector_[curScheme_->numStages_-1] = coefVector;
            }
          }else{
            EXCEPTION("The coefVector cannot be set for the chosen time scheme");
          }
          break;
        case 2:
          if(curScheme_->numSol2ndDerivs_ > 0){
            //delete the current vector and associate the new one
            //we take care if we solve for this solution order...
            delete glmVector_[curScheme_->numOldSols_+curScheme_->numSol1stDerivs_];
            glmVector_[curScheme_->numOldSols_+curScheme_->numSol1stDerivs_] = coefVector;
            avoidFreeingIdx_.insert(curScheme_->numOldSols_+curScheme_->numSol1stDerivs_);
            if(solOrder_ == 1 && avoidUpdateIdx_ == (Integer)(curScheme_->numOldSols_+curScheme_->numSol1stDerivs_)){
              stageVector_[curScheme_->numStages_-1] = coefVector;
            }
          }else{
            EXCEPTION("The coefVector cannot be set for the chosen time scheme");
          }
          break;
        default:
          EXCEPTION("The specified time derivative order invalid and the vector cannot be set!");
          break;
      }
    }
  }

  void TimeSchemeGLM::InitGLMs(){
    //initialize GLMs
    availSchemes[TRAPEZOIDAL] = new Trapezoidal(0.513);

    availSchemes[NEWMARK] = new Newmark(0.5,0.25);
  }

}
