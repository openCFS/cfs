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

   DEFINE_LOG(timeschemeglm, "timescheme.glm")


  TimeSchemeGLM::TimeSchemeGLM(GLMScheme::SchemeType type, UInt solDerivOrder, NonLinType nlType) :
         avoidUpdateIdx_(-1) {

    InitGLMs();

    curScheme_ = availSchemes[type];
    curType_ = type;

    curScheme_->solDerivOrder_ = solDerivOrder;
    solOrder_ = solDerivOrder;
    nLinType_ = nlType;
  }

  TimeSchemeGLM::TimeSchemeGLM(GLMScheme* scheme, UInt solDerivOrder, NonLinType nlType) :
           avoidUpdateIdx_(-1) {

      InitGLMs();

      curScheme_ = scheme;
      curType_ = scheme->GetType();

      curScheme_->solDerivOrder_ = solDerivOrder;
      solOrder_ = solDerivOrder;
      nLinType_ = nlType;
  }

  // Copy constructor
  TimeSchemeGLM::TimeSchemeGLM(const TimeSchemeGLM & ts) {
    InitGLMs();

    curScheme_ = ts.curScheme_;
    curType_ = curScheme_->GetType();

    curScheme_->solDerivOrder_ = ts.solOrder_;
    solOrder_ = ts.solOrder_;
    // was not copied before
    nLinType_ = ts.nLinType_;
  }
  
  TimeSchemeGLM::~TimeSchemeGLM(){
    for(UInt i=1;i<curScheme_->sizeGLMVec_;i++){
       if(avoidFreeingIdx_.find(i)==avoidFreeingIdx_.end())
         if( i < glmVector_.GetSize())
           delete glmVector_[i];
    }
    glmVector_.Clear();

    for(UInt i=1;i<curScheme_->sizeGLMVec_;i++){
       if(avoidFreeingIdx_.find(i)==avoidFreeingIdx_.end())
         if( i < initialIterGlmVector_.GetSize())
           delete initialIterGlmVector_[i];
    }
    initialIterGlmVector_.Clear();

    for(UInt i=1;i<curScheme_->numStages_;i++){
      if((Integer)i!=avoidUpdateIdx_){
        if( i < stageVector_.GetSize())
          delete stageVector_[i];
      }
    }
    stageVector_.Clear();

    std::map<GLMScheme::SchemeType, GLMScheme*>::iterator it = availSchemes.begin();
    while(it != availSchemes.end()){
      delete it->second;
      it->second = NULL;
      ++it;
    }
    availSchemes.clear();
    
    // Note: as the "curScheme_" pointer was taken initially from the 
    // availSchemes map, it was already deleted in the previous statement
    // and we must NOT delete it again.
    //delete curScheme_;
    curScheme_ = NULL;
    
    for( UInt i = 0; i < predictors_.GetSize(); ++i ) {
      delete predictors_[i];
    }
    predictors_.Clear();
    
  }

  void TimeSchemeGLM::Init(SingleVector* solVec,Double dt){   //RD: Where does the GLM get the dt from, changing it to Mathlib constand should be possible and easiest solution 
    curScheme_->adaptiveBDF2 = false;
    curScheme_->ComputeCoefficients(curScheme_->solDerivOrder_,dt);
    prevPrevSol_ = nullptr;

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

    // now we do the same for the initial copy of the glmVec
    initialIterGlmVector_.Resize(curScheme_->sizeGLMVec_);
    for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
      initialIterGlmVector_[i] = new Vector<Double>();
      initialIterGlmVector_[i]->Resize(solVec->GetSize());
      initialIterGlmVector_[i]->Init();
    }

    // first loop, we don't need to reset the glmVector just yet
    resetGlmVector_ = false;

    //now we go for the stage vector
    //if the scheme has the lastStageIsSolution_ flag set,
    //the corresponding stage vector points directly to the
    //glm Vector entry and the update step for this entry can be avoided
    stageVector_.Resize(curScheme_->numStages_);
    UInt upperBound = curScheme_->numStages_-1;
    if(curScheme_->lastStageIsSolution_){
      if(solOrder_==0){
        avoidUpdateIdx_ = 0;
      }else if(solOrder_ == 1){
        avoidUpdateIdx_ = curScheme_->numOldSols_;
      }else if(solOrder_ == 2){
        avoidUpdateIdx_ = curScheme_->numOldSols_ + curScheme_->numSol1stDerivs_ ;
      }
      stageVector_[curScheme_->numStages_-1] = glmVector_[avoidUpdateIdx_];
    }else{
      //in case of a general scheme, we need to initialize the full stage vector
      upperBound++;
    }
    for(UInt i=0;i<upperBound;i++){
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

  void TimeSchemeGLM::ModifyInit(bool extrapolateStatic){
    extrapolateStatic_ = extrapolateStatic;

//    std::string fname = "glmModifyInit_" + std::to_string(extrapolateStatic) + ".txt";
//    std::fstream myfile(fname,  std::ios::out);
//    myfile << "This is the GLM Vector before the update:" << std::endl;
//    for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
//      myfile << "Index " << i << std::endl;
//      myfile << glmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
//      myfile << "Finish GLM Vector" << std::endl;
//    }
//    myfile << std::endl;
//
//    myfile << "This is the initial GLM Vector of this time step before the update" << std::endl;
//    for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
//      myfile << "Index " << i << std::endl;
//      myfile << initialIterGlmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
//      myfile << "Finish initial GLM Vector" << std::endl;
//    }

    if ( extrapolateStatic_ ) {
      // fill all primary values with the solVec
      SingleVector * glmVec0 = glmVector_[0];
      for(UInt i=1;i<curScheme_->numOldSols_;i++){
        glmVector_[i]->Add(1.0,*glmVec0);
        glmVec0 = NULL;
      }
    }

//    myfile << std::endl;
//    myfile << std::endl;
//
//    myfile << "This is the GLM Vector after the update:" << std::endl;
//    for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
//      myfile << "Index " << i << std::endl;
//      myfile << glmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
//      myfile << "Finish GLM Vector" << std::endl;
//    }
//    myfile << std::endl;
//
//    myfile << "This is the initial GLM Vector of this time step after the update" << std::endl;
//    for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
//      myfile << "Index " << i << std::endl;
//      myfile << initialIterGlmVector_[i]->ToString(TS_NONZEROS,"\n") << std::endl;
//      myfile << "Finish initial GLM Vector" << std::endl;
//    }
//    myfile << std::endl;
//
//    myfile << std::endl;
//    myfile.close();
  }
  
  void TimeSchemeGLM::BeginStep( bool updatePredictor, bool storeInitialIterGlmVector ) {

    //-------------------------------------------------------------
    // Adaptive Timestepping Initialation:
    // Sets adaptivetimestepping Flag in GLMSCchemeLib, and checks
    // if only PDE`s for which adaptive Timestepping is Implemented
    // are used.
    //-------------------------------------------------------------
     if(domain_ != nullptr && mathparser_ == nullptr)
     {
        mathparser_ = domain_->GetMathParser();
     }

    if(domain_ != nullptr && curScheme_->adaptiveBDF2 != true)  // 
    {
      Double adaptive = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "adaptiveEnabeled");
      if(adaptive == 1)
      {
        if(curScheme_->GetType() == 3)
        {
          curScheme_->adaptiveBDF2 = true;
          adaptiveStepCount_ = 0;
        }else
        {
          EXCEPTION("Adaptive Time Stepping only implmented for, Smooth PDE and AcousticMixedPDE.");
        }
      }
    }
    //-------------------------------------------------------------

    //-------------------------------------------------------------
    // Adaptive Timestepping
    // Updates timestep, according to MathParser value and recalculates
    // the BDF2 Coefficients.
    if(curScheme_->adaptiveBDF2)
    {
      double dt = mathparser_->GetExprVars(MathParser::GLOB_HANDLER,"dt");
      curScheme_->ComputeCoefficients(curScheme_->solDerivOrder_, dt);
    }
    //-------------------------------------------------------------
  
    //update for old solutions
    if(curScheme_->usePredictors_){
      if( updatePredictor ) {
        for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
          predictors_[i]->Init();
          predictorCalculated_[i] = false;
        }
      }
    }

    if (storeInitialIterGlmVector) {
      // if we have the first iteration we store a copy of the initial glmVec for this iteration in order to be able to reset it later
      initialIterGlmVector_.Resize(curScheme_->sizeGLMVec_);
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){

        StdVector< SingleVector* > tmpGLM;
        tmpGLM.Resize(curScheme_->sizeGLMVec_);
        for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
          tmpGLM[i] = new Vector<Double>();
          tmpGLM[i]->Resize(glmVector_[i]->GetSize());
          tmpGLM[i]->Init();
        }
        for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
          SingleVector * curSVec = glmVector_[i];
          tmpGLM[i]->Add(1,*curSVec);
          curSVec = NULL;
          initialIterGlmVector_[i]->operator =((*tmpGLM[i]));
          //free memory
          delete tmpGLM[i];
        }
        tmpGLM.Clear();
        //free the memory
      }
    }
    if ( resetGlmVector_ ) {
      // reset the glmVec to the initial copy if we did not converge
      // if we get to this point, we are already in the next iteration, so it is save to assume that we did not converge
      ProcessGlmVec(false);
      resetGlmVector_ = false;
    }
  }

	void TimeSchemeGLM::UpdateStageRHSWithVector(UInt actStage, Integer derivId, SingleVector* rhsVec,
                                            SingleVector* UpdateVector, Double factor, bool forceReset){
		
		// reset rhs vec first; might be useful if only the update is wanted 
		if(forceReset){
			rhsVec->Init();
		}
		
		if(derivId > (Integer)curScheme_->maxDerivOrder_ || derivId < 0){
      return;
    }
		
		UInt dId = (UInt) derivId;
    //Calculate coefficient matrix row index
    UInt cRow = actStage * (curScheme_->maxDerivOrder_+1) + dId;
		
		UInt col = curScheme_->numStages_;
		//std::cout << "current order " << col << "  0,1 = sol itself; 2 = first time deriv; 3 = second time deriv" << std::endl;
		Double coef = curScheme_->schemeCoefs_[cRow][col];
		//std::cout << "scaling factor: " << coef << std::endl;
		if(coef !=0){
			////This is currently a HACK!!
			//if ( curType_ == GLMScheme::BDF2 && col > 0 )
			//  coef = curScheme_->schemeCoefs_[cRow][col-1];

	 //   std::cout << "NL:  cRow: " << cRow << "  col: " << col << "  Value: " << coef << std::endl;
			SingleVector * curVec = UpdateVector;
			rhsVec->Add(coef * factor,(*curVec));
		}
		
	}
	
  void TimeSchemeGLM::ComputeStageRHS(UInt actStage, Integer derivId,
                                      SingleVector* rhsVec, Integer subIdx, bool skipIncremental){

    //std::cout << "Compute StageRHS" << std::endl;
    //std::cout << "Requested derivID = " << derivId << " (0 > mult with stiffness K; 1 > mult with D; 2 > mult with M)" << std::endl;
    //if the derivative id is equal to solution we do not need to do anything
    //as the corresponding line in the GLM is equal to zero, but not! always
    if(derivId > (Integer)curScheme_->maxDerivOrder_ ||
       derivId < 0){
      rhsVec->Init();
      return;
    }
    rhsVec->Init();
    UInt dId = (UInt) derivId;
    //Calculate coefficient matrix row index
    UInt cRow = actStage * (curScheme_->maxDerivOrder_+1) + dId;
    
		/*
		 * NOTE: the predictors that are computed here, do not coincide with
		 *			 \tilde{u} and \tilde{\dot{u}} from the text book
		 *			instead, they are obtained by setting in \tilde{u} and \tilde{\dot{u}}
		 *			into the RHS formulation and then rearrange the corresponding terms
		 *			by matrices; 
		 *			> those predictors correspond to a matrix each (stiffness, damping,
		 *				mass, etc)
		 */
			
    if(curScheme_->usePredictors_ && predictorCalculated_[dId]){
      //std::cout << "Update RHS with predictors - Reuse predictors" << std::endl;
      rhsVec->Add((*predictors_[dId]));
    }else{
      //update for old solutions
      //std::cout << "Update RHS with predictors - Compute predictors" << std::endl;
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        UInt col = curScheme_->numStages_+i;

        //std::cout << "current order " << col << "  0,1 = sol itself; 2 = first time deriv; 3 = second time deriv" << std::endl;
        
        Double coef = curScheme_->schemeCoefs_[cRow][col];
        //std::cout << "scaling factor: " << coef << std::endl;
        if(coef !=0){
          //std::cout << "cRow: " << cRow << "  col: " << col << "  Value: " << coef << std::endl;
          SingleVector * curVec = glmVector_[i];
          rhsVec->Add(coef,(*curVec));
        }
      }
      if(curScheme_->usePredictors_){
        predictors_[dId]->Add(*rhsVec);
        predictorCalculated_[dId] = true;
      }
    }
    //in case of non-linear PDEs in incremental formulation, we add the stage vector
    //this is done in a simple hack and works for NEWMARK and TRAPEZOIDAL in Effective STIFFNESS ONLY!
    //NOTES: Due to the current implementation in which the AssembleNonLinRhs already computes a nonLinear-stiffness
    // matrix to the right hand side this limitation is fundamental. For future implementations one might need to think
    // of alternative ways to accomplish the non-linear solution scheme
    // furthermore, we assume, that the stageVector always holds the solution at the current non-linear iteration
    // NOT the increment. Here it becomes apparent why solveStep and Timescheme are no longer separated as initially intended
    if((nLinType_ == INCREMENTAL)&&(skipIncremental == false)){
			bool forceReset = false; 
			Double fac = -1.0;
			LOG_DBG(timeschemeglm) << "Update with stageVector_[actStage] = " << stageVector_[actStage]->ToString() << std::endl;
			UpdateStageRHSWithVector(actStage, derivId, rhsVec, stageVector_[actStage], fac, forceReset);
			
//      //std::cout << "Incremental formulation: Update RHS with -(K_eff-K)*curSol" << std::endl;
//      //std::cout << "(K_eff-K) = damping (C) and mass (M) contributions to effective stiffness matrix" << std::endl;
//      UInt col = curScheme_->numStages_;
//      //std::cout << "current order " << col << "  0,1 = sol itself; 2 = first time deriv; 3 = second time deriv" << std::endl;
//      Double coef = curScheme_->schemeCoefs_[cRow][col];
//      //std::cout << "scaling factor: " << coef << std::endl;
//      if(coef !=0){
//        ////This is currently a HACK!!
//        //if ( curType_ == GLMScheme::BDF2 && col > 0 )
//        //  coef = curScheme_->schemeCoefs_[cRow][col-1];
//
//     //   std::cout << "NL:  cRow: " << cRow << "  col: " << col << "  Value: " << coef << std::endl;
//        SingleVector * curVec = stageVector_[actStage];
//        rhsVec->Add(coef * -1.0,(*curVec));
//      }
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
    //--------------------------------------------------------------------
    // Adaptive Timestepping:
    //---------------------------------------------------------------------
    
    if (curScheme_->adaptiveBDF2) {
      if (adaptiveStepCount_ >= 2) {
        // LTE and step-size control only once enough history exists
        LTELocalErrorEstimation();

        bool accepted = ComputeAdaptiveStepSize();

        if (!accepted) {
          std::cout << " [Adaptive] Step REJECTED: LocalError= " << curScheme_->local_error_
                    << " > tol, retrying with dt= "
                    << mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "dt") << "\n";
          ResetGlmVector();
          reset_dt();
          return;
        }else
        {
          //Saves ResetValues for when a step is rerun
          curScheme_->prev_dtCurrent_ = curScheme_->dtCurrent_;  
          curScheme_->prev_dtPrev1_= curScheme_->dtPrev1_;    
          curScheme_->prev_dtPrev2_ = curScheme_->dtPrev2_;    
        }
      }
      // First 2 steps: run with user-defined dt from XML, no adaptation

      // step accepted: save y_n as y_{n-1} for next LTE computation
      if (prevPrevSol_ == nullptr) {
        prevPrevSol_ = new Vector<Double>();
        prevPrevSol_->Resize(glmVector_[1]->GetSize());
      }
      prevPrevSol_->operator=(*glmVector_[1]);
      adaptiveStepCount_++;
    }
  
    
    //just hack for flow and BDF2
    //bool usePredictorsOK = true;
    //if (curType_ == GLMScheme::BDF2 && nLinType_ == INCREMENTAL)
    //  usePredictorsOK = false;

    //update for old solutions
    //if(curScheme_->usePredictors_ && usePredictorsOK ) { //&& nLinType_ != INCREMENTAL){
    if(curScheme_->usePredictors_) { //&& nLinType_ != INCREMENTAL){
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
        }else if(nLinType_!=NONE){
          glmVector_[i] = stageVector_[0];
        }
      }

    }else{
      //first we create a temporary GLM Vector for the solution
      //which will overwrite the solution later
      StdVector< SingleVector* > tmpGLM;
      tmpGLM.Resize(curScheme_->sizeGLMVec_);
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        tmpGLM[i] = new Vector<Double>();
        tmpGLM[i]->Resize(glmVector_[i]->GetSize());
        tmpGLM[i]->Init();
      }
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        //loop over all stages
        for(UInt actS = 0; actS < curScheme_->numStages_; actS++){
          SingleVector * curSVec = stageVector_[actS];
          UInt row = curScheme_->numStages_ * (curScheme_->maxDerivOrder_+1);
          Double coef = curScheme_->schemeCoefs_[row+i][actS];
          if(coef != 0){
            tmpGLM[i]->Add(coef,*curSVec);
          }
          curSVec = NULL;
        }
        //loop over all old solutions
        for(UInt actSol = 0; actSol < curScheme_->sizeGLMVec_; actSol++){
          SingleVector * curSVec = glmVector_[actSol];
          UInt row = curScheme_->numStages_ * (curScheme_->maxDerivOrder_+1);
          Double coef = curScheme_->schemeCoefs_[row+i][curScheme_->numStages_+actSol];
          if(coef != 0){
            tmpGLM[i]->Add(coef,*curSVec);
          }
          curSVec = NULL;
        }
      }
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        glmVector_[i]->operator =((*tmpGLM[i]));
        //free memory
        delete tmpGLM[i];
      }
      tmpGLM.Clear();
      //free the memory

    }
  }

  void TimeSchemeGLM::ProcessGlmVec(bool converged){
    if (!converged) {
      // we did not converge and have to reset the glmVec
      for(UInt i=0;i<curScheme_->sizeGLMVec_;i++){
        glmVector_[i]->operator =((*initialIterGlmVector_[i]));
      }
    }
  }

  void TimeSchemeGLM::AddMatFactors(UInt stage, const std::map<FEMatrixType,Integer> & matMap,
                                        std::map<FEMatrixType,Double> & matFactors){

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
    availSchemes[GLMScheme::TRAPEZOIDAL] = new Trapezoidal(0.513);
    //availSchemes[GLMScheme::TRAPEZOIDAL] = new Trapezoidal(1.0);

    availSchemes[GLMScheme::NEWMARK] = new Newmark(0.5,0.25);
    availSchemes[GLMScheme::BDF2] = new Bdf2();
    availSchemes[GLMScheme::RK4] = new RungeKutta4();
  }

  void TimeSchemeGLM::LTELocalErrorEstimation()
  {
    // Error BDF2
    
    Double h2 = curScheme_->dtCurrent_;  // h_{n+2}
    Double h1 = curScheme_->dtPrev1_;    // h_{n+1}
    Double h0 = curScheme_->dtPrev2_;    // h_n

    Double prefactor = (h1 + h2) / 6.0;
    Double c1 = 1.0 / h2;
    Double c2 = (1.0 + h2 / h1) / h1;
    Double c3 = h2 / (h1 * h0);

    Double ErrorScheme = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "ERROR_Scheme");
    
    double l2_norm = 0.0;

    Double maxLTE = 0.0;
    UInt n = stageVector_[0]->GetSize();
    Double sum = 0.0;
    for (UInt j = 0; j < n; j++) {
        Double yNp2, yNp1, yN, yNm1;
        stageVector_[0]->GetEntry(j, yNp2);
        glmVector_[0]->GetEntry(j, yNp1);
        glmVector_[1]->GetEntry(j, yN);
        prevPrevSol_->GetEntry(j, yNm1);

        Double lte = std::abs(prefactor * (
              c1 * (yNp2 - yNp1)
            - c2 * (yNp1 - yN)
            + c3 * (yN  - yNm1)
        ));
        if(ErrorScheme == 2){sum = sum + std::pow(lte,2);}
        if (lte > maxLTE){ maxLTE = lte; }
    }

    if(ErrorScheme == 2)
    { // Normalized (Eucleadian) Error 
      l2_norm = std::sqrt(sum/n);
      curScheme_->local_error_ = l2_norm;
      mathparser_->SetValue(MathParser::GLOB_HANDLER, "MAX_LOCAL_ERROR", l2_norm);
    }else
    {
      curScheme_->local_error_ = maxLTE;
      mathparser_->SetValue(MathParser::GLOB_HANDLER, "MAX_LOCAL_ERROR", maxLTE);
    }
  }

  bool TimeSchemeGLM::ComputeAdaptiveStepSize()
  {
    Double Rtol  = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "adaptiveTol");
    Double est   = curScheme_->local_error_;
    Double h     = curScheme_->dtCurrent_;

    const Double z_U      = 0.1;
    const Double z_S      = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "adaptiveSigma") +1; // zs has to be between 1 and 2
    const Double F_U      = 10.0;
    const Double maxRatio = 1.0 + std::sqrt(2.0);  // BDF2 stability limit (growth only)

    Double h_next;
    bool   accepted;

    if (est == 0.0) {
      h_next   = F_U * h;
      accepted = true;
    } else {
      Double z = z_S * std::pow((est / Rtol), (1.0 / 3.0));
      if (z <= z_U) {
        h_next = F_U * h;  accepted = true;
      } else if (z <= z_S) {
        h_next = h / z;    accepted = true;
      } else {
        h_next = h / z;    accepted = false;
      }
    }

    std::cout << "[DEBUG] h=" << h << " est=" << est 
          << " Rtol=" << Rtol << "\n";


    // BDF2 stability: ratio h_next/h must not exceed 1+sqrt(2).
    // Shrinking is always stable for BDF2, so only the growth direction is capped.
    if (h_next / h > maxRatio)
      h_next = h * maxRatio;

    // Apply user-defined bounds after the stability cap
    Double dtMin = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "adaptiveDtMin");
    Double dtMax = mathparser_->GetExprVars(MathParser::GLOB_HANDLER, "adaptiveDtMax");
    Double h_next_unclamped = h_next;  // save before bounds to detect genuine dtMin floor hit
    h_next = std::max(h_next, dtMin);
    h_next = std::min(h_next, dtMax);

    // Force-accept only when the formula truly wanted to go below dtMin —
    // i.e. the bounds prevented further shrinking, not just any rejected step.
    if (!accepted && h_next_unclamped < dtMin)
    {
      accepted = true;
      mathparser_->SetValue(MathParser::GLOB_HANDLER, "toleranceNotReachable", 1.0);
    }
    
    mathparser_->SetValue(MathParser::GLOB_HANDLER, "dt",           h_next);
    mathparser_->SetValue(MathParser::GLOB_HANDLER, "stepRejected", accepted ? 0.0 : 1.0);

    return accepted;
  }

  void TimeSchemeGLM::reset_dt()
  {
    curScheme_->dtCurrent_ = curScheme_->prev_dtCurrent_;  
    curScheme_->dtPrev1_ =curScheme_->prev_dtPrev1_;    
    curScheme_->dtPrev2_ = curScheme_->prev_dtPrev2_;    
  }

}
