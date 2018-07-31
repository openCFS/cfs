#include <fstream>
#include <iostream>
#include <string>

#include "SolveStepHyst.hh"
#include "Driver/Assemble.hh"
#include "PDE/StdPDE.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Results/BaseResults.hh"
#include "Utils/EvalIntegrals/BiotSavart.hh"
#include "Utils/Timer.hh"
#include "Driver/SingleDriver.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/algsys/SolStrategy.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

//#include "MatVec/SingleVector.hh"

namespace CoupledField {
  
  // declare logging stream
  DECLARE_LOG(solvehyst)
  DEFINE_LOG(solvehyst, "solvehyst")
  
  SolveStepHyst::SolveStepHyst(StdPDE & apde)
  :StdSolveStep(apde)
  {
		trans_ = true;
		flagsInitialized_ = false;
    materialTensorsHystDepenedent_ = false;    
	}
	
  //! Destructor
  SolveStepHyst::~SolveStepHyst() {
    //logFile_.close();
  }
  
	void SolveStepHyst::ReadNonLinDataHyst(){
    
		PtrParamNode hystNode = solStrat_->GetHystNode();
    
    // Check, if any nonlinear node was found
    Integer deltaForm;
    if( !hystNode ) {
      WARN("Taking default parameters for nonlinear data" );
      towardsLastIteration_ = true;
      deltaMatCoupling_ = false;
      deltaMatStrain_ = false;
      nonLinMethod_ = "HYST_deltaMat_TS";
      deltaForm = -1;
      lineSearch_ = "None";
      evalDepth_ = 2;
      measurePerformance_ = false;
      testInversion_ = 0;
      incStopCrit_ = 1e-4;
      residualStopCrit_ = 1e-4;
      failBackCrit_ = -1.0;
      nonLinMaxIter_ = 35;
      
      allowReset_ = true;
      minNumberItTillReset_ = 5;
      minResDecreaseTillReset_ = 1e-5;
      nonLinMethodAfterReset_ = "HYST_deltaMat_IT";
      towardsLastIterationReset_ = true;
      deltaMatCouplingReset_ = false;
      deltaMatStrainReset_ = false;
      lineSearchAfterReset_ = "HYST_minResidual";
      useRelativeNormForRes_ = true;
      useRelativeNormForInc_ = true;
    } else {
      /*
       * New: after introducing a new hysteresis node, we have to read in parameter from that
       * node instead of from the nonlin node
       * > thus we also have to read in linesearch, method, tolerance criteria and so on
       */
      
      // solution method
      // type of line search
      towardsLastIteration_ = true;
      deltaMatCoupling_ = false;
      deltaMatStrain_ = false;
      nonLinMethod_ = "HYST_deltaMat_IT";
      
      if( hystNode->Has("method") ) {
        PtrParamNode methodNode = hystNode->Get("method");
        
        if(methodNode->Has("debug")){
          nonLinMethod_ = "HYST_debug";
        } else if(methodNode->Has("fixpoint")){
          if(methodNode->Get("fixpoint")->Has("towardsLastIteration")){
            methodNode->Get("fixpoint")->GetValue( "towardsLastIteration", towardsLastIteration_, ParamNode::PASS );
          }
          if(towardsLastIteration_){
            nonLinMethod_ = "HYST_fixPoint_IT";
          } else {
            nonLinMethod_ = "HYST_fixPoint_TS";
          }
        } else if(methodNode->Has("deltaMat")){
          if(methodNode->Get("deltaMat")->Has("towardsLastIteration")){
            //            std::cout << "towardsLastIteration found" << std::endl;
            methodNode->Get("deltaMat")->GetValue( "towardsLastIteration", towardsLastIteration_, ParamNode::PASS );
          }
          if(methodNode->Get("deltaMat")->Has("deltaMatCoupling")){
            //            std::cout << "deltaMatCoupling found" << std::endl;
            methodNode->Get("deltaMat")->GetValue( "deltaMatCoupling", deltaMatCoupling_, ParamNode::PASS );
          }
          if(methodNode->Get("deltaMat")->Has("deltaMatStrain")){
            //            std::cout << "deltaMatStrain found" << std::endl;
            methodNode->Get("deltaMat")->GetValue( "deltaMatStrain", deltaMatStrain_, ParamNode::PASS );
          }
          
          if(towardsLastIteration_){
            nonLinMethod_ = "HYST_deltaMat_IT";
          } else {
            nonLinMethod_ = "HYST_deltaMat_TS";
          }
        }
      }
      
      //      std::cout << "Read in parameter: " << std::endl;
      //      std::cout << "nonLinMethod_ = " << nonLinMethod_ << std::endl;
      //      std::cout << "towardsLastIteration = " << towardsLastIteration_ << std::endl;
      //      std::cout << "deltaMatCoupling = " << deltaMatCoupling_ << std::endl;
      //      std::cout << "deltaMatStrain = " << deltaMatStrain_ << std::endl;
      //      
      // perform logging?
      hystNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );
      
      // type of line search
      if( hystNode->Has("lineSearch") ) {
        hystNode->Get( "lineSearch")->GetValue( "type", lineSearch_,ParamNode::PASS );
      }
      
      // incremental stopping criterion
      hystNode->Get("incStopCrit")->GetValue("value",incStopCrit_, ParamNode::PASS);
      hystNode->Get("incStopCrit")->GetValue("relative",useRelativeNormForInc_, ParamNode::PASS);
      //      hystNode->GetValue( "incStopCrit", incStopCrit_, ParamNode::PASS );
      
      // residual stopping criterion
      hystNode->Get("resStopCrit")->GetValue("value",residualStopCrit_, ParamNode::PASS);
      hystNode->Get("resStopCrit")->GetValue("relative",useRelativeNormForRes_, ParamNode::PASS);
      //      hystNode->GetValue( "resStopCrit", residualStopCrit_, ParamNode::PASS );
      
      //TODO: add defaultValue to xmlSchema!!!!!
      
      //      nonLinNode->GetValue( "evalVersion", evalVersion_, ParamNode::INSERT );
      //   //   std::cout << "evalVersion: " << evalVersion_ << std::endl;
      //      nonLinNode->GetValue( "forceResidualReevaluation", forceReevaluation_, ParamNode::INSERT );
      //   //   std::cout << "forceReevaluation_: " << forceReevaluation_ << std::endl;
      
      // maximal number of NL-iterations
      hystNode->GetValue( "maxNumIters", nonLinMaxIter_, ParamNode::PASS );
      
      // abort if max number of iterations is reached?
      hystNode->GetValue("abortOnMaxNumIters",abortOnMaxIter_,ParamNode::INSERT);
      
      //      LOG_TRACE(stdsolvestep) << "Nonlinear convergence criteria were read:";
      //      LOG_DBG3(stdsolvestep) << "\tincremental Stopping Criterion: " << incStopCrit_;
      //      LOG_DBG3(stdsolvestep) << "\tresidual Stopping Criterion: " << residualStopCrit_;
      //      LOG_DBG3(stdsolvestep) << "\tmaxNumIters: " << nonLinMaxIter_;
      if( hystNode->Has("stopOnLimit") ) {
        std::string solutionString;
        // quantity
        solutionString = hystNode->Get( "stopOnLimit")
                ->Get( "quantity" )->As<std::string>();
        solutionLimit_ = SolutionTypeEnum.Parse(solutionString);
        if (feFunctions_.find(solutionLimit_) == feFunctions_.end() ) 
          EXCEPTION("ERROR: Solution type '" << solutionString << "' is not part of PDE '" << pdename_ << 
                  "' and cannot serve as stopping limit criterion");
        // minimum value
        hystNode->Get( "stopOnLimit")
        ->GetValue( "min", minValidValue_, ParamNode::PASS );
        // maximum value
        hystNode->Get( "stopOnLimit")
        ->GetValue( "max", maxValidValue_, ParamNode::PASS );
        // region
        std::string solutionRegion;
        solutionRegion = hystNode->Get( "stopOnLimit")
                ->Get( "region" )->As<std::string>();
        solutionLimitReg_ = ptgrid_->GetRegion().Parse(solutionRegion);
        if (solutionLimitReg_ != ALL_REGIONS) {
          if( subdoms_.Find(solutionLimitReg_) == -1 ) 
            EXCEPTION("ERROR: Region '" << solutionRegion <<"' is not part of PDE '" << pdename_ << 
                    "' and cannot serve as stopping limit region");
          
          // don't know how to implement region checking as we only have the solution as vector (glmvector_[0])
          EXCEPTION("checking limits in a region is not implemented (and I don't know how) (SE). " 
                  << "Just remove the 'region=XXX' tag");
        }
        //        LOG_DBG3(stdsolvestep) << "\tStop on Limit: " << solutionString;
        //        LOG_DBG3(stdsolvestep) << "\t\tmin:     " << minValidValue_;
        //        LOG_DBG3(stdsolvestep) << "\t\tmax:     " << maxValidValue_;
        //        LOG_DBG3(stdsolvestep) << "\t\tin Region:     " << solutionRegion;
        //        LOG_DBG3(stdsolvestep) << "\t\tRegionID:     " << solutionLimitReg_;
      }
      
      
      if( hystNode->Has("HYST_evalDepth") ) {
        hystNode->GetValue( "HYST_evalDepth", evalDepth_, ParamNode::INSERT );
      } else {
        evalDepth_ = 2; // extended
      }
      
      if( hystNode->Has("HYST_measurePerformance") ) {
        hystNode->GetValue( "HYST_measurePerformance", measurePerformance_, ParamNode::INSERT );
      } else {
        measurePerformance_ = 0;
      }
      
      if( hystNode->Has("HYST_failBackTol") ) {
        hystNode->GetValue( "HYST_failBackTol", failBackCrit_, ParamNode::INSERT );
      } else {
        // negative norm not possible > failback will not trigger
        failBackCrit_ = -1.0;
      }
      
      allowReset_ = true;
      minNumberItTillReset_ = 5;
      minResDecreaseTillReset_ = 1e-5;
      nonLinMethodAfterReset_ = "HYST_deltaMat_IT";
      lineSearchAfterReset_ = "HYST_minResidual";
      
      if( hystNode->Has("HYST_reset") ){
        PtrParamNode resetNode = hystNode->Get("HYST_reset");
        if(resetNode->Has("allowReset")){
          resetNode->GetValue( "allowReset", allowReset_, ParamNode::PASS );
        }
        if(resetNode->Has("minNumberItTillReset")){
          resetNode->GetValue( "minNumberItTillReset", minNumberItTillReset_, ParamNode::PASS );
        }
        if(resetNode->Has("minResDecreaseTillReset")){
          resetNode->GetValue("minResDecreaseTillReset",minResDecreaseTillReset_, ParamNode::PASS);
        }
        if(resetNode->Has("nonLinMethodAfterReset")){
          
          PtrParamNode methodNode = resetNode->Get("nonLinMethodAfterReset");
          
          if(methodNode->Has("debug")){
            nonLinMethodAfterReset_ = "HYST_debug";
          } else if(methodNode->Has("fixpoint")){
            if(methodNode->Get("fixpoint")->Has("towardsLastIteration")){
              methodNode->Get("fixpoint")->GetValue( "towardsLastIteration", towardsLastIterationReset_, ParamNode::PASS );
            }
            if(towardsLastIterationReset_){
              nonLinMethodAfterReset_ = "HYST_fixPoint_IT";
            } else {
              nonLinMethodAfterReset_ = "HYST_fixPoint_TS";
            }
          } else if(methodNode->Has("deltaMat")){
            if(methodNode->Get("deltaMat")->Has("towardsLastIteration")){
              //              std::cout << "towardsLastIteration found" << std::endl;
              methodNode->Get("deltaMat")->GetValue( "towardsLastIteration", towardsLastIterationReset_, ParamNode::PASS );
            }
            if(methodNode->Get("deltaMat")->Has("deltaMatCoupling")){
              //              std::cout << "deltaMatCoupling found" << std::endl;
              methodNode->Get("deltaMat")->GetValue( "deltaMatCoupling", deltaMatCouplingReset_, ParamNode::PASS );
            }
            if(methodNode->Get("deltaMat")->Has("deltaMatStrain")){
              //              std::cout << "deltaMatStrain found" << std::endl;
              methodNode->Get("deltaMat")->GetValue( "deltaMatStrain", deltaMatStrainReset_, ParamNode::PASS );
            }
            if(towardsLastIterationReset_){
              nonLinMethodAfterReset_ = "HYST_deltaMat_IT";
            } else {
              nonLinMethodAfterReset_ = "HYST_deltaMat_TS";
            }
          }
        }
        //          resetNode->GetValue( "nonLinMethodAfterReset", nonLinMethodAfterReset_, ParamNode::PASS );
        //        }
        if(resetNode->Has("lineSearchAfterReset")){
          resetNode->Get("lineSearchAfterReset")->GetValue( "type", lineSearchAfterReset_,ParamNode::PASS );
        }
      } 
      
			testNode_ = NULL;
      if( hystNode->Has("HYST_testOperator") ) {
				testInversion_ = 1;
        testNode_ = hystNode->Get("HYST_testOperator");
      } else {
        testInversion_ = 0;
      }
      
      /*
       * Important: evaluationDepth has to be set first as it triggers the
       * initializstion of the storage
       */
      PDE_.SetFlagInCoefFncHyst("evaluationDepth",evalDepth_);
      
      // currently we do not distinguish between -1 and +1 when checking
      // for this flag in CoefFncHyst
      // we check however, if flag != 0
      // if 0 > no deltaMat will be computed!
      if(nonLinMethod_ == "HYST_deltaMat_TS"){
        deltaForm = -1;
      } else if(nonLinMethod_ == "HYST_deltaMat_IT"){
        deltaForm = 1;
      } else {
        deltaForm = 0;
      }
      PDE_.SetFlagInCoefFncHyst("deltaForm",deltaForm);
      
      PDE_.SetFlagInCoefFncHyst("measurePerformance",measurePerformance_);
		}
    
    LOG_DBG(solvehyst) << "The following parameter were retrieved from .xml file:";
    LOG_DBG(solvehyst) << "nonLinMethod_: " << nonLinMethod_ << "( ==> deltaForm = " << deltaForm << ")";
    LOG_DBG(solvehyst) << "lineSearch_: " << lineSearch_;
    LOG_DBG(solvehyst) << "evalDepth_: " << evalDepth_;
    LOG_DBG(solvehyst) << "measurePerformance_: " << measurePerformance_;
    LOG_DBG(solvehyst) << "testInversion_: " << testInversion_;
    LOG_DBG(solvehyst) << "incStopCrit_: " << incStopCrit_;
    LOG_DBG(solvehyst) << "residualStopCrit_: " << residualStopCrit_;
    LOG_DBG(solvehyst) << "nonLinMaxIter_: " << nonLinMaxIter_;    
	}
	
	void SolveStepHyst::InitTimeStepping(){
		LOG_TRACE(solvehyst) << "SolveStepHyst::InitTimeStepping";
		// call InitTimeStepping of base class first
		StdSolveStep::InitTimeStepping();
		
		// Now initialize additional storage	
    LinRhsVec_.Resize(feFunctions_.size());
    NonLinRhsVec_.Resize(feFunctions_.size());
		predictorCorrectorUpdate_.Resize(feFunctions_.size());
		
		resVec_.Resize(feFunctions_.size());
		stageRHSUpdate_.Resize(feFunctions_.size());
		
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    UInt rhsSize = 0;
    
    //reserve memory for the rhs
    for(fncIt = feFunctions_.begin(); fncIt != feFunctions_.end();++fncIt ){
      rhsSize = fncIt->second->GetSingleVector()->GetSize();
      FeFctIdType id = fncIt->second->GetFctId();
      
			LinRhsVec_.SetSubVector(new Vector<Double>(),id);
      LinRhsVec_.GetPointer(id)->Resize(rhsSize);
			
			NonLinRhsVec_.SetSubVector(new Vector<Double>(),id);
      NonLinRhsVec_.GetPointer(id)->Resize(rhsSize);
      
			predictorCorrectorUpdate_.SetSubVector(new Vector<Double>(),id);
      predictorCorrectorUpdate_.GetPointer(id)->Resize(rhsSize);
      
			resVec_.SetSubVector(new Vector<Double>(),id);
      resVec_.GetPointer(id)->Resize(rhsSize);
			
			stageRHSUpdate_.SetSubVector(new Vector<Double>(),id);
      stageRHSUpdate_.GetPointer(id)->Resize(rhsSize);
    }
    
		LinRhsVec_.Init();
    NonLinRhsVec_.Init();
		predictorCorrectorUpdate_.Init();
		
		rhsVec_.Init();
		resVec_.Init();
		stageRHSUpdate_.Init();
    
    /*
     * Time levelflags
     * (-10 is no valid option for later usage and will trigger an initial setup
     *  in each case)
     */
    timeLevelMaterialCurrent_ = -10;
    timeLevelDeltaMatPolCurrent_ = -10;
    timeLevelDeltaMatStrainCurrent_ = -10;
    timeLevelDeltaMatCouplingCurrent_ = -10;
    //	  timeLevelDeltaMatCurrent_ = -10;
    timeLevelRHSPolCurrent_ = -10;
    timeLevelRHSStrainCurrent_ = -10;
    timeLevelRHSCouplingCurrent_ = -10;
    //    timeLevelRHSHystCurrent_ = -10;
    
    debugOutput_ = false;
  }
	
  void SolveStepHyst::SolveStepTrans() {
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
      if( fncIt->second->GetTimeScheme()->isIncremental() == false){
        //        WARN("SolveStepHyst: Incremental stepping will be enforced");
        fncIt->second->GetTimeScheme()->forceIncremental();
      }
    } 		
    
    StepTransNonLinHysteresis();
  }  
  
	Double SolveStepHyst::SetupRES(SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly, UInt iterationCounter){
		/*
     * Compute residual for current iteration k+1
     * using full evaluation, i.e.
     * 
     * res_k+1 = rhs_k+1 (including contribution from time-stepping)
     *				 - K_k+1 (std stiffness, no delta)*u_k+1
     * 
     */	
		/*
     * timeLevelMaterial
     * -2: material tensors (eps, e, h, ...) are independent of hysteresis
     * -1: material tensors are to be evaluated with hyst output from last ts (n)
     *  0: material tensors are to be evaluated with current hyst output (k)
     *  1: material tensors are to be evaluated with hyst output from previous it (k-1)
     */
		Integer timeLevelMaterial; 
		/*
     * timeLevelDeltaMat
     * -2: no delta matrix formulation
     * -1: delta matrix formulation with delta towards last ts (n) 
     *			e.g. deltaP = P_k - P_n
     *  0: --- invalid
     *  1: delta matrix formulation with delta towards previous it (k-1)
     *			e.g. deltaP = P_k - P_k-1
     */
    //		Integer timeLevelDeltaMat; 
    Integer timeLevelDeltaMatPol;
    Integer timeLevelDeltaMatStrain;
    Integer timeLevelDeltaMatCoupling;
    
		/*
     * timeLevelRHSHyst
     * -2: hyst operator will not be added to rhs
     * -1: output of hyst operator from last ts (n) will be added to rhs
     *  0: current output of hyst operator (k) will be added to rhs
     *  1: output of hyst operator from last it (k-1) will be added to rhs
     */
    //		Integer timeLevelRHSHyst;
    Integer timeLevelRHSPol;
    Integer timeLevelRHSStrain;
    Integer timeLevelRHSCoupling;
		/*
     * solutionForUpdate
     *  solutionForUpdate = solVec_ = u_k
     *			> update rhs with K_eff * u_k
     *			> incremental stepping towards last iteration
     *			> deltaU_k+1 = u_k+1 - u_k
     *	solutionForUpdate = solVecLastTS_ = u_n
     *			> update rhs with K_eff * u_n
     *			> incremental stepping towards last ts
     *			> deltaU_k+1 = u_k+1 - u_n
     */
		//SBM_Vector solutionForUpdate;
    
		/*
     * Currently only full evaluation implemented, i.e.
     * all evaluations use the current output of the hyst operator
     * but no delta-matrix
     */
		//solutionForUpdate = solVec_;
		timeLevelMaterial = 0;
    
    //		timeLevelDeltaMat = -2; > no deltaMat for residual computation > always full evaluation
    timeLevelDeltaMatPol = -2;
    timeLevelDeltaMatStrain = -2;
    timeLevelDeltaMatCoupling = -2;
    
    //		timeLevelRHSHyst = 0;
    timeLevelRHSPol = 0;
    timeLevelRHSStrain = 0;
    timeLevelRHSCoupling = 0;
    
		return SetupRESorRHS(timeLevelMaterial, timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling,
            timeLevelRHSPol, timeLevelRHSStrain, timeLevelRHSCoupling, 
            solutionForUpdate, solutionForMatrixAssembly, solutionForRHSAssembly, resVec_);
	}
	
	Double SolveStepHyst::SetupRHS(SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly, UInt iterationCounter){
		/*
     * During timestep (ts) n+1 and iteration (it) k+1
     * compute rhs vector rhs_k for INCREMENTAL stepping
     * 
     * (K_k + C + M) deltaU_k+1 = rhs_k
     *
     * Note: solution u_k is the current solution as u_k+1 has yet to be computed
     *			 therefore, k stands represents for the last known value (current value)
     *			 k-1 for the previous iteration and n for the last time step
     * 
     * For computation of rhs_k the following flags and vectors need to be set: 
     */
		/*
     * timeLevelMaterial
     * -2: material tensors (eps, e, h, ...) are independent of hysteresis
     * -1: material tensors are to be evaluated with hyst output from last ts (n)
     *  0: material tensors are to be evaluated with current hyst output (k)
     *  1: material tensors are to be evaluated with hyst output from previous it (k-1)
     */
		Integer timeLevelMaterial; 
		/*
     * timeLevelDeltaMat
     * -2: no delta matrix formulation
     * -1: delta matrix formulation with delta towards last ts (n) 
     *			e.g. deltaP = P_k - P_n
     *  0: --- invalid
     *  1: delta matrix formulation with delta towards previous it (k-1)
     *			e.g. deltaP = P_k - P_k-1
     */
    //		Integer timeLevelDeltaMat; 
    Integer timeLevelDeltaMatPol;
    Integer timeLevelDeltaMatStrain;
    Integer timeLevelDeltaMatCoupling;
		/*
     * timeLevelRHSHyst
     * -2: hyst operator will not be added to rhs
     * -1: output of hyst operator from last ts (n) will be added to rhs
     *  0: current output of hyst operator (k) will be added to rhs
     *  1: output of hyst operator from last it (k-1) will be added to rhs
     */
    //		Integer timeLevelRHSHyst;
    Integer timeLevelRHSPol;
    Integer timeLevelRHSStrain;
    Integer timeLevelRHSCoupling;
		/*
     * solutionForUpdate
     *  solutionForUpdate = solVec_ = u_k
     *			> update rhs with K_eff * u_k
     *			> incremental stepping towards last iteration
     *			> deltaU_k+1 = u_k+1 - u_k
     *	solutionForUpdate = solVecLastTS_ = u_n
     *			> update rhs with K_eff * u_n
     *			> incremental stepping towards last ts
     *			> deltaU_k+1 = u_k+1 - u_n
     * 
     * > note: solVecForUpdate_ will be set in StepTransNonLinHysteresis 
     *         according to the comment above
     */
		//SBM_Vector solutionForUpdate = solVecForUpdate_;
    
		/*   
     * Distinguish the following cases:
     */		
		if(nonLinMethod_ == "HYST_debug"){
			/*
       * 0) debug mode = fixpoint iteration without considering hysteresis at all;
       *		incremental stepping towards last iteration k:
       *			deltaU_k+1 = u_k+1 - u_k
       *			rhs_k = rhs_lin + rhs_nonlin_k + rhs_timestepping(u_k) 
       *							- K_k (std stiffness, no delta)*u_k
       */
			//solutionForUpdate = solVec_;
			timeLevelMaterial = -2;
      //			timeLevelDeltaMat = -2;
      //			timeLevelRHSHyst = -2;
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      timeLevelRHSPol = -2;
      timeLevelRHSStrain = -2;
      timeLevelRHSCoupling = -2;
      
		} else if(nonLinMethod_ == "HYST_fixPoint_IT"){
			/*
       * 1) fixpoint iteration towards last iteration k:
       *			deltaU_k+1 = u_k+1 - u_k
       *			rhs_k = rhs_lin + rhs_nonlin_k + rhs_timestepping(u_k)
       *							- K_k (std stiffness, no delta)*u_k
       *							+ rhs_hyst_k
       */
			//solutionForUpdate = solVec_;
			timeLevelMaterial = 0;
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      //			timeLevelDeltaMat = -2;
      //			timeLevelRHSHyst = 0; 
      timeLevelRHSPol = 0;
      timeLevelRHSStrain = 0;
      timeLevelRHSCoupling = 0;
		} else if(nonLinMethod_ == "HYST_fixPoint_TS"){
			/*
       * 2) fixpoint iteration towards last timestep n:
       *			deltaU_k+1 = u_k+1 - u^n
       *			rhs_k = rhs_lin + rhs_nonlin_k + rhs_timestepping(u_n)
       *							- K_k (std stiffness, no delta)*u^n
       *							+ rhs_hyst_n
       */
			//solutionForUpdate = solVecLastTS_;
			timeLevelMaterial = 0;
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      //			timeLevelDeltaMat = -2;
      //			timeLevelRHSHyst = -1;
      timeLevelRHSPol = -1;
      timeLevelRHSStrain = -1;
      timeLevelRHSCoupling = -1;
		} else if(nonLinMethod_ == "HYST_deltaMat_IT"){
			/*
       * 3) deltaMat computation towards last iteration k:
       *			deltaU_k+1 = u_k+1 - u_k
       *			rhs_k = rhs_lin + rhs_nonlin_k + rhs_timestepping(u_k)
       *							- K_k (std stiffness, no delta)*u_k
       *							+ rhs_hyst_k-1
       */
			//solutionForUpdate = solVec_;
			timeLevelMaterial = 0;
      
      // no delta mat for rhs!
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      //			timeLevelDeltaMat = -2;
      
      int timeLevelRHSHyst = 0;//1; actually k-1 (i.e. timeLevelRHSHyst = 1) would be 
      // correct according to formulas; however, according to formula we would get
      // for the case that deltaU_k+1 approx deltaU_k
      //  deltaP_k/deltaU_k * deltaU_k+1 approx deltaP_k
      //  which is added to rhs P_k-1, i.e. we get P_k on rhs 
      //  for correct solution we need P_k+1!
      // > set P_k on rhs and hope that deltaP_k+1 approx deltaP_k for deltaU_k+1 approx deltaU_k
      // > if you set this value to 1 and use the previous it value, nonlin solution will converge
      //    very slowly (if at all) and each second iteration will result in a near-zero increment
      //    (> reason: for second iteration, the same rhs will be used as for first (due to using k-1
      //                which is -1,0 during first 2 iterations, i.e. rhs_nl = 0)
      //                > no solution update can be found (as long as no other nl terms exist)
      //                > during third iteration, we will get new solution, but as the first two resulted
      //                    in the same solution, so will the stored hyst values be the same
      //                > game repeats ...)
      timeLevelRHSPol = timeLevelRHSHyst;
      
      if(deltaMatStrain_){
        timeLevelRHSStrain = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelRHSCoupling = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSCoupling = 0;
      }
      //			
		} else if(nonLinMethod_ == "HYST_deltaMat_TS"){
			/*
       * 4) deltaMat computation towards last timestep n:
       *			deltaU_k+1 = u_k+1 - u^n
       *			rhs_k = rhs_lin + rhs_nonlin_k + rhs_timestepping(u_n)
       *							- K_k (std stiffness, no delta)*u_n
       *							+ rhs_hyst_n
       * 
       */
			//solutionForUpdate = solVecLastTS_;
			timeLevelMaterial = 0;
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      //			timeLevelDeltaMat = -2;
      //			timeLevelRHSHyst = -1;
      int timeLevelRHSHyst = -1; // lastTS
      timeLevelRHSPol = timeLevelRHSHyst;
      
      if(deltaMatStrain_){
        timeLevelRHSStrain = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelRHSCoupling = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSCoupling = 0;
      }
      
		} else if(nonLinMethod_ == "HYST_deltaMat_mix"){
      /*
       * 5) HYST_deltaMat_mix
       *    1. iteration: estimateSlope, no deltaMat, hystLastTS on rhs
       *    2. iteration: deltaMat_TS, hystLastTS on rhs
       *    3+ iteration: deltaMat_IT, hystLastIt on rhs (actual k-1, see comment above)
       */
      timeLevelMaterial = 0;
      timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
      //			timeLevelDeltaMat = -2;
      int timeLevelRHSHyst = 0;
      if(iterationCounter < 2){
        
        //			timeLevelRHSHyst = -1;
        timeLevelRHSHyst = -1; // lastTS
      } else {
        //        timeLevelRHSHyst = 1; // not 0
        timeLevelRHSHyst = 1; // lastIT
      }
      
      timeLevelRHSPol = timeLevelRHSHyst;
      
      if(deltaMatStrain_){
        timeLevelRHSStrain = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelRHSCoupling = timeLevelRHSHyst;
      } else {
        // no deltaMat for strains > take current value
        timeLevelRHSCoupling = 0;
      }
      
    } else {
      std::stringstream except;
      except << "Nonlinear method " << nonLinMethod_ << " not implemented for hysteresis.";
			EXCEPTION(except.str())
		}
    
    Double normOfSystemOnly = SetupRESorRHS(timeLevelMaterial, timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling,
            timeLevelRHSPol, timeLevelRHSStrain, timeLevelRHSCoupling, 
            solutionForUpdate, solutionForMatrixAssembly, solutionForRHSAssembly, rhsVec_);
    //		SetupRESorRHS(timeLevelMaterial, timeLevelDeltaMat, timeLevelRHSHyst, solutionForUpdate, solutionForMatrixAssembly, solutionForRHSAssembly,  rhsVec_);
		algsys_->InitRHS(rhsVec_);
    return normOfSystemOnly;
	}
	
  void SolveStepHyst::AssembleSystem(SBM_Vector& solutionForMatrixAssembly, UInt iterationCounter){
    /*
     *  Assemble system according to nonLinMethod_
     *  > delta or fixpoint; activated or not  
     */
    
    /*
     * timeLevelMaterial
     * -2: material tensors (eps, e, h, ...) are independent of hysteresis
     * -1: material tensors are to be evaluated with hyst output from last ts (n)
     *  0: material tensors are to be evaluated with current hyst output (k)
     *  1: material tensors are to be evaluated with hyst output from previous it (k-1)
     */
		Integer timeLevelMaterial; 
		/*
     * timeLevelDeltaMat
     * -2: no delta matrix formulation
     * -1: delta matrix formulation with delta towards last ts (n) 
     *			e.g. deltaP = P_k - P_n
     *  0: --- invalid
     *  1: delta matrix formulation with delta towards previous it (k-1)
     *			e.g. deltaP = P_k - P_k-1
     */
		Integer timeLevelDeltaMat; 
    Integer timeLevelDeltaMatPol;
    Integer timeLevelDeltaMatStrain;
    Integer timeLevelDeltaMatCoupling;
    
    //    std::cout << "nonLinMethod_: " << nonLinMethod_ << std::endl;
    if(nonLinMethod_ == "HYST_debug"){
			/*
       * 0) debug mode = fixpoint iteration without considering hysteresis at all;
       *    > linear system has to be assembled
       */
			timeLevelMaterial = -2;
			timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
		} else if( (nonLinMethod_ == "HYST_fixPoint_IT")||(nonLinMethod_ == "HYST_fixPoint_TS") ){
			/*
       * 1,2) fixpoint iteration towards last iteration k or last ts n:
       *    > nonlinear system (mat tensors depend on hysteresis) but no delta
       */
			timeLevelMaterial = 0;
			timeLevelDeltaMatPol = -2;
      timeLevelDeltaMatStrain = -2;
      timeLevelDeltaMatCoupling = -2;
		} else if(nonLinMethod_ == "HYST_deltaMat_IT"){
			/*
       * 3) deltaMat computation towards last iteration k:
       *			deltaU_k+1 = u_k+1 - u_k
       *    > delta mat to be evaluated with respect to last iteration (k-1)
       */
			timeLevelMaterial = 0;
			timeLevelDeltaMat = 1;
      
      timeLevelDeltaMatPol = timeLevelDeltaMat;
      if(deltaMatStrain_){
        timeLevelDeltaMatStrain = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelDeltaMatCoupling = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatCoupling = 0;
      }
      
		} else if(nonLinMethod_ == "HYST_deltaMat_TS"){
			/*
       * 4) deltaMat computation towards last timestep n:
       *			deltaU_k+1 = u_k+1 - u^n
       *    > delta mat to be evaluated with respect to last timestep (n)
       */
			timeLevelMaterial = 0;
      //			timeLevelDeltaMat = -1;
      timeLevelDeltaMat = -1;
      
      timeLevelDeltaMatPol = timeLevelDeltaMat;
      if(deltaMatStrain_){
        timeLevelDeltaMatStrain = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelDeltaMatCoupling = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatCoupling = 0;
      }
		} else if(nonLinMethod_ == "HYST_deltaMat_mix"){
      /*
       * 5) HYST_deltaMat_mix
       *    1. iteration: estimateSlope, no deltaMat, hystLastTS on rhs
       *    2. iteration: deltaMat_TS, hystLastTS on rhs
       *    3+ iteration: deltaMat_IT, hystLastIt on rhs (actual k-1, see comment above)
       */
      timeLevelMaterial = 0;
      if(iterationCounter == 0){
        // timeLevelDeltaMat has to be -1 or 1; otherwise we do not call getDeltaMat which 
        // will return the estimated slope (does not work so well ...)
        timeLevelDeltaMat = -1;
        //Double steppingLength = 1e-18;
        //Double scaling = 1;
        //PDE_.EstimateCurrentSlopeForHysteresis(steppingLength, scaling);
      } else if(iterationCounter == 1){
        timeLevelDeltaMat = -1;
      } else {
        timeLevelDeltaMat = 1;
      }
      
      timeLevelDeltaMatPol = timeLevelDeltaMat;
      if(deltaMatStrain_){
        timeLevelDeltaMatStrain = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatStrain = 0;
      }
      
      if(deltaMatCoupling_){
        timeLevelDeltaMatCoupling = timeLevelDeltaMat;
      } else {
        // no deltaMat for strains > take current value
        timeLevelDeltaMatCoupling = 0;
      }
      
    } else {
      std::stringstream except;
      except << "Nonlinear method " << nonLinMethod_ << " not implemented for hysteresis.";
			EXCEPTION(except.str())
		}
    
    //    AssembleSystem(timeLevelMaterial, timeLevelDeltaMat, solutionForMatrixAssembly);
    AssembleSystem(timeLevelMaterial, timeLevelDeltaMatPol,timeLevelDeltaMatStrain,timeLevelDeltaMatCoupling , solutionForMatrixAssembly);
  }
  
  void SolveStepHyst::AssembleSystem(Integer timeLevelMaterial, Integer timeLevelDeltaMatPol, Integer timeLevelDeltaMatStrain, Integer timeLevelDeltaMatCoupling, 
          SBM_Vector& solutionForMatrixAssembly){
    //    std::cout << "#### AssembleSystem - 1" << std::endl;
    
    // TODO: remove timeLevelDeltaMat
    int timeLevelDeltaMat = timeLevelDeltaMatPol;
		if(debugOutput_){
			LOG_DBG(solvehyst) << "AssembleSystem";
			LOG_DBG(solvehyst) << "timeLevelMaterial (-2 = deactivated; -1 = lastTS; 0 = curIT; +1 = lastIT): " << timeLevelMaterial; 
      //	  LOG_DBG(solvehyst) << "timeLevelDeltaMat (-2 = deactivated; -1 = lastTS; 0 = INVALID; +1 = lastIT): " << timeLevelDeltaMat; 
			LOG_DBG(solvehyst) << "solutionForAssembly: " << solutionForMatrixAssembly.ToString();
			LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString();
		}
		/*
     *  Assemble system according to flags timeLevelDeltaMat, timeLevelMaterial
     *  > this is used for res and rhs assembly as here the system might be different
     *    from the one that actually has to be solved
     */
    SBM_Vector diff;
    diff = lastUpdateVecForMatrices_;
    diff.Add(-1.0,solutionForMatrixAssembly);
    //    std::cout << "#### AssembleSystem - 2" << std::endl;
    bool solutionChanged = (diff.NormL2() > 1e-16);
    bool deltaMatStatusChanged = ( (timeLevelDeltaMatPol != timeLevelDeltaMatPolCurrent_) ||
    (timeLevelDeltaMatStrain != timeLevelDeltaMatStrainCurrent_) ||
    (timeLevelDeltaMatCoupling != timeLevelDeltaMatCouplingCurrent_) );
    
    bool materialStatusChanged = (timeLevelMaterial != timeLevelMaterialCurrent_);
    //    std::cout << "#### AssembleSystem - 2.1" << std::endl;
    // system depends on solution, if
    // a) the material tensors depend on hysteresis (e.g. the coupling tensors)
    // b) we have a deltaMatrix formulation (timeLevelDeltaMat == -1 or +1)
    bool systemSolutionDependent = ( materialTensorsHystDepenedent_ || (abs(timeLevelDeltaMat) == 1) );
		//if(debugOutput_){
    LOG_DBG(solvehyst) << "Assemble matrices - solutionChanged? " << solutionChanged;
    LOG_DBG(solvehyst) << "Assemble matrices - deltaMatStatusChanged? " << deltaMatStatusChanged;
    LOG_DBG(solvehyst) << "Assemble matrices - materialStatusChanged? " << materialStatusChanged;
    LOG_DBG(solvehyst) << "Assemble matrices - systemSolutionDependent? " << systemSolutionDependent;
		//}
    //      std::cout << "#### AssembleSystem - 2.2" << std::endl;
		// we need to reassemble, if
    // a) we switched deltaMat status, i.e. from non-delta to delta form
    // b) if material status changed, i.e. from non-dependent to depenedent (should normally not occur)
    // c) if solution has changed and system is solution dependent
    if( deltaMatStatusChanged || materialStatusChanged || (systemSolutionDependent && solutionChanged) ){
      LOG_DBG(solvehyst) << "Assemble matrices - Try to perform reassembly ";
      //    bool systemRequiresReassembly = (diff.NormL2() < 1e-16); // '<' should have been '>' 
      //    
      //    if( (timeLevelDeltaMat != timeLevelDeltaMatCurrent_) ||
      //				(timeLevelMaterial != timeLevelMaterialCurrent_) ||
      //				(systemRequiresReassembly == true) ) {
      //			LOG_TRACE(solvehyst) << "SolveStepHyst::Assemble matrices";
      //      std::cout << "#### AssembleSystem - 2.3" << std::endl;
      // evalPurpose 1 > assemble
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
			PDE_.SetFlagInCoefFncHyst("SetTimeLevelMaterial",timeLevelMaterial);
      //			PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMat",timeLevelDeltaMat);
      
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatPol",timeLevelDeltaMatPol);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatStrain",timeLevelDeltaMatStrain);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatCoupling",timeLevelDeltaMatCoupling);
      
      // important: allow reevaluation of hyst operator > might have changed
      //      std::cout << "#### AssembleSystem - 2.4" << std::endl;
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);
      //			std::cout << "#### AssembleSystem - 3" << std::endl;
      debugOutput_ = false;
			if(debugOutput_){
				// NEVER assign directly! SBM_Vector diff = solVec_ will assign a reference!
				// computing diff-soltutionForAssembly will modify solVec_!!!!
				// instead we need two steps SBM_Vector = diff; diff = solVec_ > copy value!
				SBM_Vector diff;
				diff = solVec_;
				diff.Add(-1.0,solutionForMatrixAssembly);
				
				LOG_DBG(solvehyst) << "SolVec: " << solVec_.ToString();
				LOG_DBG(solvehyst) << "diff: " << diff.ToString();
			}	
			
			SBM_Vector solBackup;
			solBackup = solVec_;
      //      std::cout << "#### AssembleSystem - 3.1" << std::endl;
			solVec_ = solutionForMatrixAssembly;
			// AssembleMatrices accesses solVec_ for assembly > set solVec_ accordingly
      //      std::cout << "#### AssembleSystem - 3.2" << std::endl;
      assemble_->AssembleMatrices(false); 
			//algsys_->GetRHSVal( NonLinrhsVec_ );
			// reset solVec_
      //      std::cout << "#### AssembleSystem - 3.3" << std::endl;
			solVec_ = solBackup;
      //      std::cout << "#### AssembleSystem - 4" << std::endl;
			timeLevelMaterialCurrent_ = timeLevelMaterial;
			timeLevelDeltaMatPolCurrent_ = timeLevelDeltaMatPol;
      timeLevelDeltaMatStrainCurrent_ = timeLevelDeltaMatStrain;
      timeLevelDeltaMatCouplingCurrent_ = timeLevelDeltaMatCoupling;
      lastUpdateVecForMatrices_ = solutionForMatrixAssembly;
		}
    //      std::cout << "#### AssembleSystem - END" << std::endl;
  }
  
	Double SolveStepHyst::SetupRESorRHS(Integer timeLevelMaterial, Integer timeLevelDeltaMatPol, Integer timeLevelDeltaMatStrain, Integer timeLevelDeltaMatCoupling, 
          Integer timeLevelRHSPol, Integer timeLevelRHSStrain, Integer timeLevelRHSCoupling, 
          SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly, 
          SBM_Vector& returnVector){
		
    // TODO: remove this flag
    //    int timeLevelRHSHyst = timeLevelRHSPol;
    int timeLevelDeltaMat = timeLevelDeltaMatPol;
    
		if(debugOutput_){
			LOG_DBG(solvehyst) << "SetupRESorRHS";
			LOG_DBG(solvehyst) << "timeLevelMaterial (-2 = deactivated; -1 = lastTS; 0 = curIT; +1 = lastIT): " << timeLevelMaterial; 
			LOG_DBG(solvehyst) << "timeLevelDeltaMat (-2 = deactivated; -1 = lastTS; 0 = INVALID; +1 = lastIT): " << timeLevelDeltaMat; 
			LOG_DBG(solvehyst) << "solutionForUpdate: " << solutionForUpdate.ToString();
			LOG_DBG(solvehyst) << "solutionForMatrixAssembly: " << solutionForMatrixAssembly.ToString();
      LOG_DBG(solvehyst) << "solutionForRHSAssembly: " << solutionForRHSAssembly.ToString();
			LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString();
		}
    //    std::cout << "#### SetupRESorRHS::IT - 1" << std::endl;
		/*
     * 0. Init returnVector
     */
		returnVector.Init();
		
		/*
     * 1. Assemble system if required
     * > here we have to use the version with prescribed timelevels as those
     *   might differ from the ones that are required by nonLinMethod
     */
    //    std::cout << "#### SetupRESorRHS::IT - 1-2" << std::endl;
    AssembleSystem(timeLevelMaterial, timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling, solutionForMatrixAssembly);
    //    std::cout << "#### SetupRESorRHS::IT - 2" << std::endl;
		/*
     * 2. Setup Linear part and predictor corrector update
     */
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Assembly of linear RHS/RES required? "<< LinRHSRequiresAssembly_;
		}
		if(LinRHSRequiresAssembly_){
      //			LOG_DBG(solvehyst) << "SolveStepHyst::Assemble linear RHS/RES";
      
			LinRhsVec_.Init();
      algsys_->InitRHS(LinRhsVec_);
      assemble_->AssembleLinRHS();
      algsys_->GetRHSVal( LinRhsVec_ );
      
			LinRHSRequiresAssembly_ = false;
		}
		returnVector.Add(1.0,LinRhsVec_);
    //    std::cout << "#### SetupRESorRHS::IT - 3" << std::endl;
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Return vector - linearRHS: " << returnVector.ToString();
		}
		
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Compute predictor/corrector update required? "<< PredictorCorrectorRequiresAssembly_;
		}		
		if(PredictorCorrectorRequiresAssembly_){
      //		LOG_DBG(solvehyst) << "SolveStepHyst::Compute predictor/corrector update";
      
			predictorCorrectorUpdate_.Init();
      algsys_->InitRHS(predictorCorrectorUpdate_);
			
			if(trans_){
				// remember: (currently) we only support single stage timestepping!
				UInt stage = 0;
				
				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
				std::map<FEMatrixType,Integer>::iterator matIt;
				
				// update RHS according to time stepping
				// here we only want to get the update due to predictor corrector, but
				// not the update for incremental stepping (this has to be done separately
				// as the value for update differs)
				bool skipIncremental = true;
				
				for(matIt = matrices.begin();matIt != matrices.end();matIt++){
					if(matIt->second < 0)
						continue;
          
					for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
						FeFctIdType fctId = fncIt->second->GetFctId();
						fncIt->second->GetTimeScheme()->ComputeStageRHS(stage,matIt->second,stageRHS_.GetPointer(fctId),-1,skipIncremental);
					}
          
					algsys_->UpdateRHS(matIt->first,stageRHS_,true);
				}	
			}
			algsys_->GetRHSVal( predictorCorrectorUpdate_ );
			
			PredictorCorrectorRequiresAssembly_ = false;
		}
		returnVector.Add(1.0,predictorCorrectorUpdate_);
    //    std::cout << "#### SetupRESorRHS::IT - 4" << std::endl;
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Return vector - linearRHS+predictorCorrector: " << returnVector.ToString();
		}
    
		
		/*
     * 3. Setup nonlin RHS with hysteresis
     * > nonlin RHS is solution dependent 
     * > if current solution differes from previously used one >  reassemble
     */
    SBM_Vector diff;
    diff = lastUpdateVecForNonLinRHS_;
    diff.Add(-1.0,solutionForRHSAssembly);
    bool NonLinRHSRequiresReassembly = (diff.NormL2() > 1e-16);
		
		//if(debugOutput_){
    LOG_DBG(solvehyst) << "NonLinRHSRequiresReassembly?" << NonLinRHSRequiresReassembly;
    LOG_DBG(solvehyst) << "CurrentSolution for RHS assembly: " << solutionForRHSAssembly.ToString();
    LOG_DBG(solvehyst) << "LastUsedSolution for assembly: " << lastUpdateVecForNonLinRHS_.ToString();
    LOG_DBG(solvehyst) << "Diff input: " << diff.ToString();
    LOG_DBG(solvehyst) << "NonLinRhsVec_ pre recomputation: " << NonLinRhsVec_.ToString();
		//}		
    diff = NonLinRhsVec_;
    //    std::cout << "#### SetupRESorRHS::IT - 5" << std::endl;
    /*
     * NOTE: assembleNonLinRHS has its own check if solution did change or not
     *      > in case that the solution vector is the same, no reassembly is done
     *      > time level might be different though (e.g. 1 of iteration k+1 == 0 of iteration k!)
     *          > in that case, we enter this case, set NonLinRhsVec_ to 0 and skip assembly
     *          > rhs vector wrong!
     *      > do no longer check for timeLevel!
     *      > same issue for assemble matrices?!
     * > Edit: we may not skip the timeLevelcheck; we rather have to force a reassembly!
     *        > reason: solution might be the same as explained above; however, the
     *          hysteresis rhs part might require reassembly! 
     *          remember: at timestep k+1 using deltaMat_IT we need hyst_k-1 on rhs!
     *                    > this would have been detected by the timelevel comparison
     *                      but that lead to a missing reassembly due to same solution
     *                      (note: the solution is not required to obtain hyst_k-1 as
     *                             hyst_k-1 is stored in the hyst operator and just needs
     *                             to be loaded)
     *                 
     */
    bool forceHystReassembly = ( (timeLevelRHSPol != timeLevelRHSPolCurrent_) ||
    (timeLevelRHSStrain != timeLevelRHSStrainCurrent_) ||
    (timeLevelRHSCoupling != timeLevelRHSCouplingCurrent_) );
    
    //    LOG_DBG(solvehyst) << "Requested Timelevel of RHSHyst: "<<timeLevelRHSHyst;
    //    LOG_DBG(solvehyst) << "Last Used Timelevel of RHSHyst: "<<timeLevelRHSHystCurrent_;
    //		if(timeLevelRHSHyst != timeLevelRHSHystCurrent_) {
    //      forceHystReassembly = true;
    //      LOG_DBG(solvehyst) << "Force reassembly as we need different hysteresis state.";
    //    }
    if( NonLinRHSRequiresReassembly || forceHystReassembly ){
      //	LOG_DBG(solvehyst) << "SolveStepHyst::Assemble non-linear RHS/RES";
      
			NonLinRhsVec_.Init();
			algsys_->InitRHS(NonLinRhsVec_); //load storage into algsys
      // evalPurpose 1 > assemble
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
      //			PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSHyst",timeLevelRHSHyst);	
      
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSPol",timeLevelRHSPol);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSStrain",timeLevelRHSStrain);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSCoupling",timeLevelRHSCoupling);
      
      // important: allow reevaluation of hyst operator > might have changed
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);
			
			SBM_Vector solBackup;
			solBackup = solVec_;
			solVec_ = solutionForRHSAssembly;
			// AssembleNonLinRHS accesses solVec_ for assembly > set solVec_ accordingly
      
      LOG_DBG(solvehyst) << "Try to assemble nonlin rhs";
      assemble_->AssembleNonLinRHS(); 
			algsys_->GetRHSVal( NonLinRhsVec_ );
			// reset solVec_
			solVec_ = solBackup;
      
      timeLevelRHSPolCurrent_ = timeLevelRHSPol;
      timeLevelRHSStrainCurrent_ = timeLevelRHSStrain;
      timeLevelRHSCouplingCurrent_ = timeLevelRHSCoupling;
      //			timeLevelRHSHystCurrent_ = timeLevelRHSHyst;
      lastUpdateVecForNonLinRHS_ = solutionForRHSAssembly;
		}
		returnVector.Add(1.0,NonLinRhsVec_);
    //		std::cout << "#### SetupRESorRHS::IT - 6" << std::endl;
    LOG_DBG(solvehyst) << "NonLinRhsVec_ after recomputation: " << NonLinRhsVec_.ToString();
    diff.Add(-1.0,NonLinRhsVec_);
    LOG_DBG(solvehyst) << "Diff NonLinRhsVec_: " << diff.ToString();
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Return vector - linearRHS,predictorCorrector,nonLinRHS: " << returnVector.ToString();
		}
    
		/*
     * 4. Update res/rhs according to update solution vector
     */
    SBM_Vector systemUpdate;
    systemUpdate = rhsVec_; // take any value to get correct size
    systemUpdate.Init();
    // note: compute KeffTimesSolution does not assemble system! it takes the current system
    // as assembled above
    ComputeKeffTimesSolution(solutionForUpdate, systemUpdate);
    Double normOfSystemUpdate = systemUpdate.NormL2();
    //    std::cout << "#### SetupRESorRHS::IT - 7" << std::endl;
    returnVector.Add(-1.0,systemUpdate);
    //    
    //		// load returnVector into algsys; then update for incremental formulation
    //		algsys_->InitRHS(returnVector);
    //
    //		// Stiffness part
    //		//LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u";
    //		//LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 1: K";
    //		solutionForUpdate.ScalarMult(-1.0);
    //		algsys_->UpdateRHS(STIFFNESS,solutionForUpdate,true);
    //		if(!trans_){
    //			algsys_->UpdateRHS(STIFFNESS_UPDATE,solutionForUpdate,true);
    //		}
    //		solutionForUpdate.ScalarMult(-1.0);
    //
    //		algsys_->GetRHSVal( returnVector );
    //		
    //		if(debugOutput_){
    //			LOG_DBG(solvehyst) << "Return vector - linearRHS,predictorCorrector,nonLinRHS,solUpdate(stiff): " << returnVector.ToString();
    //		}
    //		
    //		// Damping, Mass part (with corresponding scaling > done via timestepping)
    //		if(trans_){
    //			LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 2: C,M";
    //			std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    //			std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    //			std::map<FEMatrixType,Integer>::iterator matIt;
    //			
    //			bool forceReset = true;
    //			UInt stage = 0;
    //			for(matIt = matrices.begin();matIt != matrices.end();matIt++){
    //				if(matIt->second < 0)
    //					continue;
    //
    //				// write scaled solution vector (scaled by timestepping factors) to stageRHSUpdate_
    //				for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
    //					FeFctIdType fctId = fncIt->second->GetFctId();
    //					fncIt->second->GetTimeScheme()->UpdateStageRHSWithVector(stage,matIt->second,stageRHSUpdate_.GetPointer(fctId),
    //																						solutionForUpdate.GetPointer(fctId), -1.0, forceReset);
    //
    //				}
    //				// now update loaded rhs with matrix*scaled solution vector
    //				algsys_->UpdateRHS(matIt->first,stageRHSUpdate_,true);
    //			}
    //		}
    //		/*
    //		 * 5. Retrieve vector
    //		 */
    //		algsys_->GetRHSVal( returnVector );
    //		
		if(debugOutput_){
			LOG_DBG(solvehyst) << "Return vector - linearRHS,predictorCorrector,nonLinRHS,solUpdate(full): " << returnVector.ToString();
		}
    //		std::cout << "#### SetupRESorRHS::END" << std::endl;
    //	LOG_DBG(solvehyst) << "returnVector: " << returnVector.ToString(); 
    return normOfSystemUpdate;
	}
  
  Double SolveStepHyst::LineSearchHyst(SBM_Vector& currentSol, SBM_Vector& solIncrement,UInt iterationCounter)  {
    //
    LOG_DBG(solvehyst) << "Perform linesearch with method " << lineSearch_;
    SBM_Vector solBak;
    solBak = solVec_;
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
    
    if(lineSearch_ == "HYST_quadApprox"){
      /*
       * Approximate resiudal by quadratic equation 
       * > take minimum of approximation
       * 
       * res(eta) = a0 + a1*eta + a2*eta^2
       * dres/deta = a1 + 2a2*eta
       * > etaExtremum = -a1/2a2
       * > etaExtremum is min if 2a2 > 0
       * 
       * res(0) = a0
       * res(0.5) = a0 + a1/2 + a2/4
       * res(1) = a0 + a1 + a2
       * 
       * > a0 = res(0)
       * > a1 = -3 res(0) + 4 res(0.5) - res(1)
       * > a2 = 2 res(0) - 4 res(0.5) + 2 res(1)
       * 
       */
      Double a1,a2;
      Double res0,res05,res1;
      
      SBM_Vector tmpSol;
      tmpSol = currentSol;
      tmpSol.Add(0,solIncrement);
      
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      
      // SetupRes will write the residual directly to resVec_
      // It will also reassemble the system! > we have to call this function 
      // again after the linesearch to get the curret system set up!
      SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
      res0 = resVec_.NormL2();
      
      tmpSol = currentSol;
      tmpSol.Add(0.5,solIncrement);
      solVec_ = tmpSol;
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      
      SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
      res05 = resVec_.NormL2();
      
      tmpSol = currentSol;
      tmpSol.Add(1.0,solIncrement);
      solVec_ = tmpSol;
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      
      SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
      res1 = resVec_.NormL2();
      
      // restore solution vector
      solVec_ = solBak;
      
      //      a0 = res0;
      a1 = -3*res0 + 4*res05 - res1;
      a2 = 2*res0 - 4*res05 +2*res1;
      
      //      UInt solCase = 0;
      Double etaSol = 1;
      if(a2 > 0){
        // minimum
        //        solCase = 1;
        etaSol = -a1/(2*a2);
      } else {
        // either a2 == 0 or a2 < 0 i.e. extremum would be a maximum
        // take minimum of res05 and res1
        if(res05 < res1){
          //          solCase = 2;
          etaSol = 0.5;
        } else {
          //          solCase = 3;
          etaSol = 1.0;
        }
      }
      
      //      std::cout << "Solution case: " << solCase << std::endl;
      //      std::cout << "Found eta: " << etaSol << std::endl;
      return etaSol;
    } else if( ( lineSearch_ == "HYST_minResidual") || (lineSearch_ == "HYST_minResidual_largerValues") ){
      /*
       * Search for eta by trial and error;
       * take eta that leads to smallest residual norm
       */
      const UInt nrEtas = 5;
      Double eta[nrEtas];
      if(lineSearch_ == "HYST_minResidual"){
        eta[0] = 1.0;
        eta[1] = 0.9;
        eta[2] = 0.7;
        eta[3] = 0.3;
        eta[4] = 0.1;
        //        Double eta[nrEtas] = {1.0,0.9,0.7,0.3,0.1};
      } else {
        eta[0] = 1.0;
        eta[1] = 0.9;
        eta[2] = 0.75;
        eta[3] = 0.5;
        eta[4] = 0.25;
      }
      
      Double bestNorm = 1e16;
      Double bestEta = 1.0;
      Double resNorm = 0.0;
      SBM_Vector tmpSol;
      for(UInt i = 0; i < nrEtas; i++){
        //  LOG_DBG(solvehyst) << "Eta: " << eta[i];
        tmpSol = currentSol;
        tmpSol.Add(eta[i],solIncrement);
        
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        solVec_ = solBak;
        
        // SetupRes will write the residual directly to resVec_
        // It will also reassemble the system! > we have to call this function 
        // again after the linesearch to get the curret system set up!
        SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
        resNorm = resVec_.NormL2();
        //  LOG_DBG(solvehyst) << "resNorm: " << resNorm;
        if(resNorm < bestNorm){
          bestNorm = resNorm;
          bestEta = eta[i];
        }
      }
      return bestEta;
      
    } else if( (lineSearch_ == "HYST_overRelaxation") || (lineSearch_ == "HYST_bisect") ){
      /*
       * Project residual onto solIncrement; try to 
       * reduce the value of this projection to 0
       */
      Double etaUp;
      if((lineSearch_ == "HYST_overRelaxation")){
        // allow to step further than increment 
        etaUp = 2.0;
      } else {
        etaUp = 1.0;
      }
      //Double etaLow = 0.001;
      Double etaLow = 1e-5;
      Double etaMid = (etaUp+etaLow)/2;
      Double projLow,projUp;
      //Double resNormUp,resNormLow;
      
      SBM_Vector tmpSol;
      tmpSol = currentSol;
      tmpSol.Add(etaUp,solIncrement);
      
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      solVec_ = solBak;
      
      // SetupRes will write the residual directly to resVec_
      // It will also reassemble the system! > we have to call this function 
      // again after the linesearch to get the curret system set up!
      SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
      resVec_.Inner(solIncrement,projUp);
      //resNormUp = resVec_.NormL2();
      
      tmpSol = currentSol;
      tmpSol.Add(etaLow,solIncrement);
      
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      solVec_ = solBak;
      
      SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
      resVec_.Inner(solIncrement,projLow);
      //resNormLow = resVec_.NormL2();
      UInt numSteps = 6;
      if(projLow*projUp >= 0){
        /*
         * both projections have the same sign or are 0
         * > take the eta that leads to smaller projection (or residual?!)
         * > projection works better
         */
        //if(abs(resNormLow) < abs(resNormUp)){
        //        if(abs(projLow) < abs(projUp)){
        //          return etaLow;
        //        } else {
        //          return etaUp;
        //        }
        // reduce lower eta and try to find change in sign
        for(UInt i = 0; i < numSteps; i++){
          etaLow = etaLow/2.0;
          tmpSol = currentSol;
          tmpSol.Add(etaLow,solIncrement);
          
          solVec_ = tmpSol;
          // set solution vector to tmp solution, then evaluate hyst operators
          // do not reevaluate matrices for inversion > flag = 0
          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
          solVec_ = solBak;
          
          SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
          resVec_.Inner(solIncrement,projLow);
          
          if(projLow*projUp < 0){
            return etaLow;
          }       
        }
        
        // could not find changing sign; take eta with smallest projection
        if(abs(projLow) < abs(projUp)){
          return etaLow;
        } else {
          return etaUp;
        }
        
      } else {
        /*
         * Here we start the actual bisection algorithm
         * limit the number iterations as it might take too many of them to
         * actually reduce the projection close to 0
         */
        
        Double tol = 1e-7;
        Double projMid;
        for(UInt i = 0; i < numSteps; i++){
          etaMid = (etaUp + etaLow)/2.0;
          
          tmpSol = currentSol;
          tmpSol.Add(etaMid,solIncrement);
          
          solVec_ = tmpSol;
          // set solution vector to tmp solution, then evaluate hyst operators
          // do not reevaluate matrices for inversion > flag = 0
          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
          solVec_ = solBak;
          
          SetupRES(tmpSol,tmpSol,tmpSol,iterationCounter);
          resVec_.Inner(solIncrement,projMid);
          
          if(abs(projMid) < tol){
            // solution good enough
            // take etaMid
            return etaMid;
          }
          if( projMid*projLow >= 0 ){
            // mid and low have same sign > search between mid and up
            etaLow = etaMid;
            projLow = projMid;
          } else {
            etaUp = etaMid;
            projUp = projMid;
          }
        }
      }
      return etaMid;
    } else {
      std::stringstream warnmsg;
      warnmsg << "Linesearch method " << lineSearch_ << " not implemented for hysteresis. ";
      warnmsg << "No linesearch is performed.";
      WARN(warnmsg.str());
      return 1.0;
    }
    
  }
  
  void SolveStepHyst::ComputeKeffTimesSolution(SBM_Vector& solutionForUpdate, SBM_Vector& returnVector){
    
    /*
     * Function: Take CURRENT system (no reassembly done here!), multiply it with solution
     *           and store the result in returnVector
     */
    
    Double factor = 1.0;
    returnVector.Init();
    algsys_->InitRHS(returnVector);
    
    // Stiffness part
    //LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u";
    //LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 1: K";
    solutionForUpdate.ScalarMult(factor);
    algsys_->UpdateRHS(STIFFNESS,solutionForUpdate,true);
    if(!trans_){
      algsys_->UpdateRHS(STIFFNESS_UPDATE,solutionForUpdate,true);
    }
    solutionForUpdate.ScalarMult(factor);
    
    algsys_->GetRHSVal( returnVector );
    
    // Damping, Mass part (with corresponding scaling > done via timestepping)
    if(trans_){
      //LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 2: C,M";
      std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
      std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
      std::map<FEMatrixType,Integer>::iterator matIt;
      
      bool forceReset = true;
      UInt stage = 0;
      for(matIt = matrices.begin();matIt != matrices.end();matIt++){
        if(matIt->second < 0)
          continue;
        
        // write scaled solution vector (scaled by timestepping factors) to stageRHSUpdate_
        for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
          FeFctIdType fctId = fncIt->second->GetFctId();
          fncIt->second->GetTimeScheme()->UpdateStageRHSWithVector(stage,matIt->second,stageRHSUpdate_.GetPointer(fctId),
                  solutionForUpdate.GetPointer(fctId), factor, forceReset);
          
        }
        // now update loaded rhs with matrix*scaled solution vector
        algsys_->UpdateRHS(matIt->first,stageRHSUpdate_,true);
      }
    }
    /*
     * Retrieve vector
     */
    algsys_->GetRHSVal( returnVector );
    
  }
  
  void SolveStepHyst::StepTransNonLinHysteresis() {
		static UInt timestepCnt = 0;
		timestepCnt++;
		LOG_TRACE(solvehyst) << "##################################### " << timestepCnt;
		LOG_TRACE(solvehyst) << "#### SolveStepHyst::Begin of timestep " << timestepCnt;
		LOG_DBG(solvehyst) << "SolVec_ at start of timestep: " << solVec_.ToString(); 
		LOG_DBG(solvehyst) << "resVec_ at start of timestep: " << resVec_.ToString(); 
    //		std::cout << "#### SolveStepHyst::Start of timestep" << std::endl;
		/*
     * Initial checks
     */
		if(!flagsInitialized_){
			ReadNonLinDataHyst();
			flagsInitialized_ = true;
      materialTensorsHystDepenedent_ = PDE_.MaterialTensorsHystDependent();
		}
		
		if(testInversion_ == 1){
			PDE_.TestInversionOfHystOperator(testNode_);
			testInversion_ = 0;
		}
    
		
		//obtain the number of stages
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 ){
      /*
       * maybe it is quite easy to enable multiple stages, maybe it does not work;
       * feel free to find out
       */
      EXCEPTION("SolveStepHyst: just one stage time-stepping allowed");
    }
    UInt stage = 0;
    
    /*
     * Setup and initialize solution vectors and rhs
     */
    SBM_Vector solInc(BaseMatrix::DOUBLE);
		solInc.Init();
    
    SBM_Vector gatheredSolInc(BaseMatrix::DOUBLE);
		gatheredSolInc.Init();
    
    SBM_Vector bestSol(BaseMatrix::DOUBLE);
		bestSol.Init();
    Double bestResNorm = 1e15;
    //SBM_Vector actSol(BaseMatrix::DOUBLE);
    //actSol = solVec_;
    
    // solVec_ currently holds the solution from the previous timestep
    // create backup for later usage
    solVecLastIT_ = solVec_;
    solVecLastTS_ = solVec_;
    lastUpdateVecForNonLinRHS_ = solVecLastTS_;
    lastUpdateVecForMatrices_ = solVecLastTS_;
    //LOG_DBG(solvehyst) << "lastUpdateVecForNonLinRHS_: " << lastUpdateVecForNonLinRHS_.ToString() ;
    //LOG_DBG(solvehyst) << "lastUpdateVecForMatrices_: " << lastUpdateVecForMatrices_.ToString() ;
    
    /*
     * Compute norm of linear rhs as reference value for relative res. error
     */
    rhsVec_.Init();
    algsys_->InitRHS(rhsVec_);
    
    // setup right hand side
    //    Double loadFactor = 1.0;
    //    Double RhsLinL2Norm = SetLinRHS(loadFactor);
    Double normOfSystem; // norm of K_eff(u_k+1)u_k+1
    /*
     * Init time stepping
     */
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
    
    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }
    
    /*
     * Important: Although stageSol seems not to be necessary, the TimeScheme will rely on it
     * It HAS to store the total solution (not the increment!) -> see TimeSchemeGLM.cc
     * Important2: Incremental stepping is only working (correctly) for EffectiveStiffness formulation
     */
    SBM_Vector stageSol;
    stageSol.Resize(feFunctions_.size());
		
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      FeFctIdType fctId = fncIt->second->GetFctId();
      stageSol.SetSubVector(fncIt->second->GetTimeScheme()->GetStageVector(stage),fctId);
      fncIt->second->GetTimeScheme()->InitStage(stage,actTime_,PDE_.GetDomain());
    }
    stageSol.SetOwnership(false);
    
    /*
     * During the later computations, we will add all incremental solutions to
     * stageSol; unlike the case of StepTransNonLin, we add the results to the
     * current solution instead of only gathering the increments.
     * Note that we initialize stageSol with actSol in StepTransNonLin, too, but
     * there, we will overwrite stageSol with 0.0 during the first iteration
     */
    // stageSol = actSol;
    // solVec_  = actSol;
    stageSol = solVec_;
		
    /*
     * disable output of switching and rotation state as BMP images (even if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",false);
    
    /*
     * Setup flags
     */
    
    static UInt overallIterationCounter = 0;
    UInt iterationCounter = 0;
    UInt iterationCounterLog = 0;
    bool performOneMoreStep;
    Double solL2Norm;
    Double residualErr;
    Double incrementalErr;
    Double residualErrRel;  
    Double incrementalErrRel; 
    
    // keep track of residual; if it is not decreasing over several iterations > reset (see below)
    UInt numResets = 0;
    std::deque<Double> resNormHistory(minNumberItTillReset_, 50000);
    
    /*
     * Reset flags
     * Note: 
     * - system and nonlinRHS are solution depended > need to be evaluated each 
     *   iteration
     * - linRHS is solution independent and predictor corrector uses only the
     *   values from last ts; BUT: predictor-corrector in effective stiffness
     *   formulation requires damping and mass matrices from the current iteration,
     *   i.e. if those are solution dependent, the predictor-corrector parts need
     *   reassembly, too
     *   For hysteresis: neither C, nor M are solution dependend (only K_elec and K_mag are)
     *    > no reevaluation needed
     */
    LinRHSRequiresAssembly_ = true;
    PredictorCorrectorRequiresAssembly_ = true;
    debugOutput_ = false;
		
    // DEBUG -- FORCE delta last ts
    //nonLinMethod_ = "HYST_deltaMat_TS";
    //nonLinMethod_ = "HYST_deltaMat_mix";
    
    bool steppingTowardsLastTS = false;
    if( (nonLinMethod_ == "HYST_fixPoint_TS")||(nonLinMethod_ == "HYST_deltaMat_TS") ){
      steppingTowardsLastTS = true;
    }
    
    // after reset, we step towards zero;
    // if steppingTowardsLastTS = true, we start from 0 in each new iteration
    // else we only start from 0 once
    bool steppingTowardsZero = false;
    
    // currently, we allow only incremental stepping, i.e. we compute only
    // the increment to the last ts solution; in that sense, the dirichlet
    // boundary conditions have to be the difference of current values and the
    // values from the last ts
    bool deltaIDBC;
    bool setIDBC;
    bool lastTS = false;
    std::string nonLinBackup = nonLinMethod_;
    std::string lineSearchBackup = lineSearch_;
    bool towardsLastIterationBackup = towardsLastIteration_;
    bool deltaMatCouplingBackup = deltaMatCoupling_;
    bool deltaMatStrainBackup = deltaMatStrain_;
    
    /*
     *  Begin iterations
     */
    Double lastTSSatAvg, lastItSatAvg, curItSatAvg, oppositeDirAsTSAvg, oppositeDirAsItAvg;
    
    /*
     * Before doing the first solve step > estimate slope
     * Background:
     *  the first iteration is especially crucial for the solution process as it
     *  has to be a full step (i.e. etaLinesearch = 1.0) in order to satisfy the boundary conditions
     * 
     *  this means, however, if the initial slope is too small (or high) the first step might go
     *  so bad (and cannot be scaled down) such that the solution process never can recover again
     * 
     *  to avoid this, we estimate the slope at the current point; this slope of course is not
     *  perfect either, but it may avoid to dramatic steppings
     * 
     */
    // note: steppingLength will be multiplied with xSaturated during calculation
    // Double steppingLength = 1e-18;
    // Double scaling = 1;
    // Estimated slope will only be used if timeLevelDeltaMat = -1 or 1, i.e. even for FP iteration, we have
    // to set this to -1 or 1 for the first iteration if we want to use the estimated slope
    //PDE_.EstimateCurrentSlopeForHysteresis(steppingLength, scaling);
    
    //LOG_TRACE(solvehyst) << "Set previous hyst values for deltaIT prior to first iteration to current value";
    // first deltaMat will be 0, though
    //PDE_.SetPreviousHystVals(lastTS);
    
    // the following variable indicates towards which previous state the deltaMatrix shall be computed
    // (only in case of deltaMatrix and deltaMat towards Iteration of course)
    UInt numIterationsTillUpdateOfLastITValues = 1;
    
    /*
     * For magnetic PDE, we have the problem, that the matrix that is used
     * to compute H (i.e. mu/nu) might depend on H/M itself
     * > compute mu/nu for inversion only once each iteration
     */
    // allow setting of matrices
    PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
    // EvaluateHystOperators; flag = 3 > set Matrix for Inverion and rotation states for possible vector extension of scalar model
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",3);
    
    do {
      //      std::cout << "#### SolveStepHyst::Start of IT" << std::endl;
      
      if(debugOutput_){
        PDE_.CheckSaturationOfHystOperators(lastTSSatAvg, lastItSatAvg, curItSatAvg, oppositeDirAsTSAvg, oppositeDirAsItAvg);
        
        LOG_DBG(solvehyst) << "Saturation information: ";
        LOG_DBG(solvehyst) << lastTSSatAvg*100.0 << " percent of elements were saturated during last timestep";
        LOG_DBG(solvehyst) << lastItSatAvg*100.0 << " percent of elements were saturated during last iteration";
        LOG_DBG(solvehyst) << curItSatAvg*100.0 << " percent of elements are currently saturated";
        LOG_DBG(solvehyst) << oppositeDirAsTSAvg*100.0 << " percent of elements changed direction of saturation by more than 90 degree compared to last timestep";
        LOG_DBG(solvehyst) << oppositeDirAsItAvg*100.0 << " percent of elements changed direction of saturation by more than 90 degree compared to last iteration";
      }
      
      /*
       * Step 0: Increase counter 
       */
      overallIterationCounter++;
      iterationCounter++;
      iterationCounterLog++;
      LOG_TRACE(solvehyst) << "## Iteration: " << iterationCounter ;
      
      /*
       * Step 0.2: Reset solution UPDATE vector (if necessary)
       * (this vector is the one that gets subtracted from the rhs, i.e. the
       * one the stepping is going to)
       * > NOTE: do not overwrite solVec_ here, as the system still needs to be
       *         assembled using the current solution!
       * 
       */
      //      if( (steppingTowardsLastTS)&&(!steppingTowardsZero) ){
      //        // increment has to be towards last ts > solVecForUpdate_ has to hold this
      //        // value at begin of iteration
      //        solVecForUpdate_ = solVecLastTS_;
      //      } else if( (steppingTowardsLastTS)&&(steppingTowardsZero) ){
      //        // increment has to be towards 0 > solVecForUpdate_ has to be reset to 0
      //        solVecForUpdate_.Init();
      //      } else {
      //        // stepping towards last iteration > solVecForUpdate_ = solVec_
      //        solVecForUpdate_ = solVecLastIT_;
      //      }
      
      if(steppingTowardsZero){
        solVecForUpdate_.Init();
      } else if(steppingTowardsLastTS){
        solVecForUpdate_ = solVecLastTS_;
      } else {
        // DOES NOT WORK IF COMMENTED OUT, reasons not clear yet
        //if((iterationCounterLog-1)%numIterationsTillUpdateOfLastITValues == 0){
        solVecForUpdate_ = solVec_; // solVecLastIT_; < solVecLastIT is actually the solution from two iterations ago
        // i.e. during iteration k+1, solVec = sol_k and solVecLastIT = sol_k-1
        //}
      }
      
			if(debugOutput_){
				LOG_DBG(solvehyst) << "Solution vectors at start of iteration: ";
				LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
				LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
				LOG_DBG(solvehyst) << "solVecLastIT_: " << solVecLastIT_.ToString(); 
				LOG_DBG(solvehyst) << "solVecForUpdate_: " << solVecForUpdate_.ToString(); 
			}
      //      std::cout << "#### SolveStepHyst::IT - 1" << std::endl;
      //      
      //      LOG_DBG(solvehyst) << "Try following nonlin method: " << nonLinMethod_;
      //      LOG_DBG(solvehyst) << "Flags: ";
      //      LOG_DBG(solvehyst) << "steppingTowardsLastTS: " << steppingTowardsLastTS;
      //      LOG_DBG(solvehyst) << "steppingTowardsZero: " << steppingTowardsZero;
      
      /*
       * Step 1: setup rhs
       */
      //debugOutput_ = true;
      LOG_TRACE(solvehyst) << "SetupRHS";
      SetupRHS(solVecForUpdate_,solVec_,solVec_,iterationCounter);
      //      std::cout << "#### SolveStepHyst::IT - 1-2" << std::endl;
      //debugOutput_ = false;
      LOG_DBG(solvehyst) << "RHS to solve for: : " << rhsVec_.ToString(); 
      
      /*
       * Step 2: (re-)assemble matrices according to nonLinMethod
       * Note: During SetupRHS() the matrices already get assembled, but 
       *        most probably with different flags > reassembly needed
       */
      // LOG_DBG(solvehyst) << "AssembleSystem";
			AssembleSystem(solVec_, iterationCounter);
      //      std::cout << "#### SolveStepHyst::IT - 2" << std::endl;
      /*
       * Step 3: setup system, including boundary conditions, solver, etc.
       */
      matrix_factor_.clear();
      algsys_->InitMatrix(SYSTEM);
      for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
        FeFctIdType fctId = fncIt->second->GetFctId();
        fncIt->second->GetTimeScheme()->AddMatFactors(0,matrices,matrix_factor_[fctId]);
        algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
      }
      
      PDE_.SetBCs();
      algsys_->BuildInDirichlet();
      algsys_->SetupPrecond();
      algsys_->SetupSolver();
      
      /*
       * In the following cases, we must set DBCs
       * a) Incremental stepping towards last iteration
       *     a1) for iterationCounter == 1
       *        > for later iterations, the contribution from IDBC is already con-
       *          tained in solVec
       *     a2) after reset (i.e. setting the solution vector back to 0 and start from there
       * b) Incremental stepping towards last timestep
       *     b1) for EACH iteration, as we reset the solution to solVecLastTS_ again
       *     b2) after reset
       * 
       */
      if(steppingTowardsLastTS || steppingTowardsZero || ( iterationCounter == 1 && couplingIter_ == 0 ) ){
        setIDBC = true;
      } else {
        setIDBC = false;
      }
      
      if(steppingTowardsZero){
        // if we step towards zero, we do not subtract the old idbc values
        deltaIDBC = false;
      } else {
        deltaIDBC = true;
      }
      
      if( (steppingTowardsZero)&&(!steppingTowardsLastTS) ){
        // if stepping towards 0 was actived due to a reset, check if
        // we have stepping towards last ts or towards last it
        // > case a) towards last it > reset flag again
        // > case b) towards last ts > for all further iterations go back to 0
        steppingTowardsZero = false;
      }
      
      /*
       * STORE CURRENT VALUES FOR NEXT ITERATION
       * > this has to be done AFTER setting up the system
       * > otherwise if we compute the delta between current and old iteration
       *    this delta would be 0 all the time
       * > only needed for computation of incremental error
       */
      solVecLastIT_ = solVec_;
      //      std::cout << "#### SolveStepHyst::IT - 3" << std::endl;
      lastTS = false;
      // Note: by calling SetPreviousHystVals with flag lastTS = false,
      //        coefFncHyst will evaluate the hysteresis operator with
      //        the current state of solVec_; for this evaluation, the memory of
      //        the hysteresis operator will be locked, i.e. the evaluation will
      //        lead to reversible changes
      
      // NEW: set previous hyst values only, if solution did change significantly
      //      otherwise delta in the denominator gets very small
      
      
      
      //      OLD VERSION: SetPreviousHystValues evaluated hyst operators and stored results
      //      accrodingly; this meant, that this step had to be done after assembly where the
      //      old values were still required and before the new solution was found
      //              > NEW VERSION: take the preevaluated values and copy them to the corresponding
      //              arrays; this has to be done before evaluating the final result for the iteration
      //      
      //      /*
      //       * for setup of RHS and deltaMat towards IT, the value stored in coefFctHyst
      //       * will be used; if we want to have a delta over multiple iterations, we just 
      //       * have to skip the setting of the previous hyst values; by this, the old it
      //       * values will not be overwritten
      //       */
      //      if((iterationCounterLog-1)%numIterationsTillUpdateOfLastITValues == 0){
      //        //if(gatheredSolInc.NormL2() > 1e-5){
      //          LOG_TRACE(solvehyst) << "Set previous hyst values";
      //          PDE_.SetPreviousHystVals(lastTS);
      ////          gatheredSolInc.Init();
      ////        } else {
      ////          LOG_TRACE(solvehyst) << "Do not set previous hyst values as gathered solInc is yet too small";
      ////        }
      //      }
      // } 
      
      /*
       * Step 4: solve system and retrieve the solution (-update)
       */
      //LOG_DBG(solvehyst) << "Solve";
      algsys_->Solve(setIDBC,deltaIDBC);
      //      std::cout << "#### SolveStepHyst::IT - 4" << std::endl;
      /*
       *  RETRIEVE SOLUTION INCREMENT
       *  - nonLinMethod_ = "HYST_debug","HYST_fixPoint_IT","HYST_deltaMat_IT" > increment towards last iteration
       *  - nonLinMethod_ = "HYST_fixPoint_TS","HYST_deltaMat_TS" > increment towards last ts
       *  - steppingTowardsZero > increment towards 0 > full step
       */
      algsys_->GetSolutionVal(solInc, setIDBC, deltaIDBC );
      ////      std::cout << "#### SolveStepHyst::IT - 5" << std::endl;
      //      /*
      //       * Short check of solution:
      //       *  > without reassembly anything, the found vector solInc should
      //       *    reduce the vector
      //       *      rhs - Keff*solInc to 0
      //       *  > check this
      //       * 
      //       * > DOES somehow change the result > why?
      //       */
      //      if(debugOutput_){
      //        LOG_TRACE(solvehyst) << "Check actual solution vector against system";
      //        SBM_Vector rhsBackup;
      //        SBM_Vector tmpVector;
      //        rhsBackup = rhsVec_;
      //        if( (nonLinMethod_ == "HYST_deltaMat_IT")||(nonLinMethod_ == "HYST_deltaMat_TS") ){
      //          LOG_TRACE(solvehyst) << " Current solution strategy used deltaMatrix, i.e. we solved ";
      //          LOG_TRACE(solvehyst) << "  K_eff_delta_k * deltaU_k+1 = rhs_k - rhs_hyst_k - K_eff_k * U_k ";
      //          LOG_TRACE(solvehyst) << "1. Check if:  rhs_k - rhs_hyst_k - K_eff_k * U_k - K_eff_delta_k * deltaU_k+1 approx 0 ";
      //          // note: according to math, we should have rhs_hyst_k-1 on rhs but see SetupRHS for more details why k is used
      //          tmpVector = rhsVec_; // get copy just to obtain correct size
      //          tmpVector.Init();
      //          // note: Keff is K_eff_delta_k at this point
      //          ComputeKeffTimesSolution(solInc, tmpVector);
      //          // note: rhsVec_ = rhs_k - rhs_hyst_k - K_eff_k * U_k at this point
      //          
      //          tmpVector.Add(-1.0,rhsBackup);
      //          LOG_TRACE(solvehyst) << " > result:  rhs_k - rhs_hyst_k - K_eff_k * U_k - K_eff_delta_k * deltaU_k+1 = " << tmpVector.ToString();
      //          LOG_TRACE(solvehyst) << " > result:  |rhs_k - rhs_hyst_k - K_eff_k * U_k - K_eff_delta_k * deltaU_k+1| = " << tmpVector.NormL2();
      //          
      //          // for our solution process we have hyst_k-1 on rhs as the difference to hyst_k is hidden in deltaMatrix;
      //          // if we entangle the deltaMatrix to the standard matrix without delta, we have to move that difference to hyst_k-1, i.e. turning it
      //          // to rhs_hyst_k
      //          LOG_TRACE(solvehyst) << "2. Check if:  rhs_k - rhs_hyst_k - K_eff_k * U_k+1 approx 0 ";
      //          tmpVector.Init();
      //         
      //          SBM_Vector solForAssemblyTMP;
      //          SBM_Vector solForRHSTMP;
      //          SBM_Vector solForUpdateTMP;
      //          // we want to assemble matrix and rhs with u_k as solution; at this point solVec_ is already overwritten
      //          // by system due to solution process
      //          // > take backup from start of iteration
      //          solForAssemblyTMP = solVecForUpdate_;
      //          solForRHSTMP = solVecForUpdate_;
      //          // to compute K_eff_k*U_k+1 we set solForUpdateTMP to U_k+1
      //          solForUpdateTMP = solVecForUpdate_;
      //          solForUpdateTMP.Add(1.0,solInc);
      //          
      //          // 0 > evaluate with given solution vector solutionForAssembly; -2 deactivate deltaMat
      //          Integer TL_Mat = 0;
      //          Integer TL_DeltaMat = -2;
      //          Integer TL_RHSHyst = 0;
      //          SetupRESorRHS(TL_Mat, TL_DeltaMat,  TL_RHSHyst, solForUpdateTMP, solForAssemblyTMP, solForRHSTMP, tmpVector);
      //          
      //          LOG_TRACE(solvehyst) << " > result:  rhs_k - rhs_hyst_k - K_eff_k * U_k+1 = " << tmpVector.ToString();
      //          LOG_TRACE(solvehyst) << " > result:  |rhs_k - rhs_hyst_k - K_eff_k * U_k+1| = " << tmpVector.NormL2();
      //          
      //          LOG_TRACE(solvehyst) << "3. Compute actual residual (full evaluation, no linesearch), i.e.:";
      //          LOG_TRACE(solvehyst) << "  res_k+1 = rhs_k+1 - rhs_hyst_k+1 - K_eff_k+1 * U_k+1";
      //          tmpVector.Init();
      //         
      //          solForUpdateTMP = solVecForUpdate_;
      //          solForUpdateTMP.Add(1.0,solInc);
      //          solForAssemblyTMP = solForUpdateTMP;
      //          // 0 > evaluate with given solution vector solutionForAssembly; -2 deactivate deltaMat
      //          TL_Mat = 0;
      //          TL_DeltaMat = -2;
      //          TL_RHSHyst = 0;
      //          SetupRESorRHS(TL_Mat, TL_DeltaMat, TL_RHSHyst, solForUpdateTMP, solForAssemblyTMP, solForAssemblyTMP, tmpVector);
      //          
      //          LOG_TRACE(solvehyst) << " > result:  res_k+1 = rhs_k+1 - rhs_hyst_k+1 - K_eff_k+1 * U_k+1 = " << tmpVector.ToString();
      //          LOG_TRACE(solvehyst) << " > result:  |res_k+1| = |rhs_k+1 - rhs_hyst_k+1 - K_eff_k+1 * U_k+1| = " << tmpVector.NormL2();
      //        }
      //       rhsVec_ = rhsBackup;
      //      }
      //LOG_DBG(solvehyst) << "solInc (computed): " << solInc.ToString();
      
      /*
       * The solution vector will be the increment towards solVecForUpdate_,
       * i.e. the new solVec will be solVecForUpdate + eta*solInc
       */
      solVec_ = solVecForUpdate_;
      //      std::cout << "#### SolveStepHyst::IT - 6" << std::endl;
      /*
       * Step 5: check solution by evaluating residual and increment
       */
      Double residualL2Norm = 0.0;
      Double incrementL2Norm = 0.0;
      Double etaLinesearch;
      
			// unlock hyst operator (if locked)
			PDE_.SetFlagInCoefFncHyst("lockHystOperator",0);
			
      // before we evaluate the hyst operator with the actual solution, 
      // store its old value for next iteration
      // > do this BEFORE linesearch as operators get evaluated during LS, too
      if(lastTS){
        PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);
      } else {
        PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
      }
      
      //LOG_DBG(solvehyst) << "SetupRES/Perform LS";
      if( setIDBC || ( lineSearch_ == "none" ) ){
        /*
         * Linesearch is only allowed if the increment does not hold contributions
         * of the Dirichlet boundary as in that case, we have to step with length 1
         * to get the boundary terms correct!
         * > during first iteration for incremental towards last iteration
         * > during each iteration for incremental towards last ts or 0
         */
        etaLinesearch = 1.0;
      } else {
        etaLinesearch = LineSearchHyst(solVec_, solInc,iterationCounter);
      }
      solVec_.Add(etaLinesearch,solInc);
      
      /*
       * New solution has been found
       * > evaluate hyst operator
       */
      // allow setting of matrices
      // NOTE: matrices will be set with values from last iteration so make
      // sure that SetPreviousHystValues was called before
      PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
      // EvaluateHystOperators; flag = 1 > set Matrix for Inverion, too
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",1);
      
      //LOG_DBG(solvehyst) << "Updated solution: " << solVec_.ToString();
      // although LineSearchHyst will setup the residual, we have to
      // call this function once more with the actual best eta to ensure
      // that the system is assemble correctly!
      //debugOutput_ = true;
      LOG_TRACE(solvehyst) << "SetupRESIDUAL with found etaLinesearch: " << etaLinesearch;
      //LOG_DBG(solvehyst) << "OLD Nonlinear RHS (i.e. Hystpart) before setup: " << NonLinRhsVec_.ToString();
      normOfSystem = SetupRES(solVec_,solVec_,solVec_,iterationCounter);
      //debugOutput_ = false;
      //LOG_DBG(solvehyst) << "NEW Nonlinear RHS (i.e. Hystpart) after setup: " << NonLinRhsVec_.ToString();
      //LOG_DBG(solvehyst) << "Remaining Residual: " << resVec_.ToString();
      // update actSol to currently computed value
      // > is actSol needed at all?
			// actSol = solVec_;
      
      //      std::cout << "residualNorm: " << resVec_.NormL2() << std::endl;
      //      std::cout << "normOfSystem: " << normOfSystem << std::endl;
      //            
      residualL2Norm = resVec_.NormL2();
      //      if ( RhsLinL2Norm > 1.0 ){
      //        residualErr = residualL2Norm / RhsLinL2Norm;
      //      } else {
      residualErr = residualL2Norm;
      //      }
      
      if ( normOfSystem > 0.0 ){
        residualErrRel = residualErr/normOfSystem;
      } else {
        std::cout << "normOfSystem = 0!" << std::endl;
        residualErrRel = residualErr;
      }
      
      //      std::cout << "relError: " << residualErrRel << std::endl;
      
      /*
       * For incremental error checking, compute the actual update between
       *  the current and the previous iteration
       */
      solInc = solVec_;
      solInc.Add(-1.0,solVecLastIT_);
      
      gatheredSolInc.Add(1.0,solInc);
      
      //debugOutput_ = true;
			if(debugOutput_){
				LOG_DBG(solvehyst) << "solInc (towards lastIT): " << solInc.ToString();
      }
      //debugOutput_ = false;
      
			
      incrementL2Norm = solInc.NormL2();
      solL2Norm  = solVec_.NormL2();
      
      //      if ( solL2Norm > 1.0 ){
      //        incrementalErr = incrementL2Norm / solL2Norm;
      //      } else {
      incrementalErr = incrementL2Norm;
      //      }
      
      if ( solL2Norm > 0.0 ){
        incrementalErrRel = incrementalErr/solL2Norm;
      } else {
        //        std::cout << "solVec = 0!" << std::endl;
        incrementalErrRel = incrementalErr;
      }
      Double incrementalErrCheck;
      Double residualErrCheck;
      
      if(useRelativeNormForRes_){
        residualErrCheck = residualErrRel;
      } else {
        residualErrCheck = residualErr;
      }
      
      if(residualErrCheck < bestResNorm){
        bestResNorm = residualErrCheck;
        bestSol = solVec_;
      }
      
      resNormHistory.push_back(residualErrCheck);
      resNormHistory.pop_front();
      
      if(useRelativeNormForInc_){
        incrementalErrCheck = incrementalErrRel;
      } else {
        incrementalErrCheck = incrementalErr;
      }
           
      //      std::cout << "Use useRelativeNormForResetCheck_? " << useRelativeNormForResetCheck_ << std::endl;
      //      std::cout << "Use useRelativeNormForRes_? " << useRelativeNormForRes_ << std::endl;
      //      std::cout << "Use useRelativeNormForInc_? " << useRelativeNormForInc_ << std::endl;     
      //      
      performOneMoreStep =
              (incrementalErrCheck > incStopCrit_) || (residualErrCheck > residualStopCrit_);
      
      //      std::cout << "incrementalErr: " << incrementalErr << std::endl;
      //      std::cout << "incrementalErr < crit: " << (incrementalErr < incStopCrit_) << std::endl;
      //      std::cout << "incrementalErrRel: " << incrementalErrRel << std::endl;
      //      std::cout << "incrementalErrRel < crit: " << (incrementalErrRel < incStopCrit_) << std::endl;
      //      
      //      std::cout << "#### SolveStepHyst::IT - 8" << std::endl;
      
      // output of norms and data to info.xml
      if ( nonLinLogging_ == true ) {
        // get current step
        UInt actStep = PDE_.GetSolveStep()->GetActStep();
        
        if (PDE_.IsIterCoupled() == false) {
          couplingIter_ = 0;
        }
        
        WriteHystIterToInfoXML(pdename_,
                nonLinMethod_,
                lineSearch_,
                couplingIter_,
                actStep,
                iterationCounterLog,
                iterationCounter,
                residualErr, residualErrRel, useRelativeNormForRes_, 
                incrementalErr, incrementalErrRel, useRelativeNormForInc_,
                etaLinesearch, false, false);
        
        //        if (PDE_.IsIterCoupled()) {
        //          WriteHystIterToInfoXML(pdename_,nonLinMethod_, couplingIter_, actStep,iterationCounterLog,iterationCounter, residualErrCheck, incrementalErrCheck, etaLinesearch);
        //          WriteHystIterToInfoXML("relErrors","->", couplingIter_,actStep,iterationCounterLog,iterationCounter, residualErrRel, incrementalErrRel, etaLinesearch);
        //          WriteHystIterToInfoXML("absErrors","->", couplingIter_,actStep,iterationCounterLog,iterationCounter, residualErr, incrementalErr, etaLinesearch);
        //        } else {
        //          //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
        //          WriteHystIterToInfoXML(pdename_,nonLinMethod_, actStep,iterationCounterLog,iterationCounter, residualErrCheck, incrementalErrCheck, etaLinesearch);
        //          WriteHystIterToInfoXML("relErrors","->", actStep,iterationCounterLog,iterationCounter, residualErrRel, incrementalErrRel, etaLinesearch);
        //          WriteHystIterToInfoXML("absErrors","->", actStep,iterationCounterLog,iterationCounter, residualErr, incrementalErr, etaLinesearch);
        //        }
        // write norm to file
        logFile_ <<  iterationCounterLog << "\t"
                << residualErr << "\t"
                << incrementalErr << "\t"
                << etaLinesearch << std::endl;
      }
      
      // compare residual at front and back of deque; 
      // residual should have decreased (at least a bit); otherwise try reset
      if( (performOneMoreStep) && (resNormHistory.front()-minResDecreaseTillReset_ <= resNormHistory.back()) ){
				//debugOutput_ = true;
        if((numResets == 0)&&allowReset_){
					//continue;
          LOG_TRACE(solvehyst) << "Residual did not decrease by at least " << minResDecreaseTillReset_ << " using solution method " << nonLinMethod_; 
          nonLinMethod_ = nonLinMethodAfterReset_;
          towardsLastIteration_ = towardsLastIterationReset_;
          deltaMatCoupling_ = deltaMatCouplingReset_;
          deltaMatStrain_ = deltaMatStrainReset_;
          lineSearch_ = lineSearchAfterReset_;
          numIterationsTillUpdateOfLastITValues = 1;
          LOG_TRACE(solvehyst) << "Switch to  " << nonLinMethod_ << " with LS " << lineSearch_; 
          LOG_TRACE(solvehyst) << "Additionally, compute deltaMat over " << numIterationsTillUpdateOfLastITValues << " itartions"; 
          
          if((nonLinMethod_ == "HYST_deltaMat_IT")||(nonLinMethod_ == "HYST_fixPoint_IT")){
            steppingTowardsLastTS = false;
          } else {
            steppingTowardsLastTS = true;
          }
          
          // reset to remanence is not so good as tests have shown; 
          bool goToRemanence = false;
          if(goToRemanence){
            LOG_TRACE(solvehyst) << "Reset solVec to 0, put hystoperators to REMANENCE and try solution method HYST_deltaMat_IT"; 
            PDE_.SetFlagInCoefFncHyst("forceRemanence",1);
            PDE_.SetFlagInCoefFncHyst("lockHystOperator",1);
          } else {
            LOG_TRACE(solvehyst) << "Reset solVec to 0 and try solution method HYST_deltaMat_IT"; 
          }
          // important: we have to evaluate the hyst operator with the reset solution value
          // to be able to compute correct deltaMatrices
          // for this evaluation, we have to take care of:
          // a) correct storage 
          //      for stepping towards lastIT, we have to set the storage for lastIT
          //      for stepping towards lastTS, we have to set the storage for lastTS
          //  > this is adjusted by the flag lastTS (as used prior)
          // b) memory of hyst operator must not be overwritten (we do not want up with
          //    some minor loops due to iterations)
          //  > force memory lock (required if lastTS = true)	
          solVec_.Init();
          
          // sol vec was reset > evaluate hyst operator
          // + reset matrix for inversion
          //      > here with flag 2 > set solution to 0 first, then set matrix, then eval hyst operator
          // allow setting of matrices
          PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",2);
          
          
          lastTS = steppingTowardsLastTS;
          
          if(lastTS){
            PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);
          } else {
            PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
          }
          //bool forceMemoryLock = true;
          //PDE_.SetPreviousHystVals(lastTS,forceMemoryLock);
          
          //PDE_.EstimateCurrentSlopeForHysteresis(steppingLength, scaling);
          //LOG_DBG(solvehyst) << "Iteration " << iterationCounter << " of timestep " << timestepCnt;
          
          // reset deque
          for(UInt i = 0; i < resNormHistory.size(); i++){
            resNormHistory[i] = 50000;
          }
          steppingTowardsZero = true;
          numResets++;
          iterationCounterLog = 0;
          
          if ( nonLinLogging_ == true ) {
            // get current step
            UInt actStep = PDE_.GetSolveStep()->GetActStep();
            
            //            if (PDE_.IsIterCoupled()) {
            //              WriteHystIterToInfoXML(pdename_,"Reset", couplingIter_, actStep,0,0, -1.0, -1.0, 0);
            //            } else {
            //              //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
            //              WriteHystIterToInfoXML(pdename_,"Reset", actStep,0,0, -1.0, -1.0, 0);
            //            }
            
            WriteHystIterToInfoXML(pdename_,
                    nonLinMethod_,
                    lineSearch_,
                    couplingIter_,
                    actStep,
                    iterationCounterLog,
                    iterationCounter,
                    residualErr, residualErrRel, useRelativeNormForRes_, 
                    incrementalErr, incrementalErrRel, useRelativeNormForInc_,
                    etaLinesearch, true, false);
            
            // write norm to file
            logFile_ <<  iterationCounterLog << "\t"
                    << residualErr << "\t"
                    << incrementalErr << "\t"
                    << etaLinesearch << std::endl;
          }
          
        } else {
          // check for failback tolerance
          // > check best solution so far
          LOG_TRACE(solvehyst) << "Checking best found solution against failBackCriterion";
          if(bestResNorm < failBackCrit_){
            LOG_TRACE(solvehyst) << "Best found solution so far passes failBack > take it";
            solVec_ = bestSol;
            residualErr = bestResNorm;
            // output of norms and data to info.xml
            if ( nonLinLogging_ == true ) {
              // get current step
              UInt actStep = PDE_.GetSolveStep()->GetActStep();
              
              //              if (PDE_.IsIterCoupled()) {
              //                WriteHystIterToInfoXML(pdename_,"Failback", couplingIter_, actStep,0,0, residualErrCheck, -1.0, 0);
              //              } else {
              //                //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
              //                WriteHystIterToInfoXML(pdename_,"Failback", actStep,0,0, residualErrCheck, -1.0, 0);
              //              }
              
              WriteHystIterToInfoXML(pdename_,
                      nonLinMethod_,
                      lineSearch_,
                      couplingIter_,
                      actStep,
                      iterationCounterLog,
                      iterationCounter,
                      residualErr, residualErrRel, useRelativeNormForRes_, 
                      incrementalErr, incrementalErrRel, useRelativeNormForInc_,
                      etaLinesearch, false, true);
              
              // write norm to file
              logFile_ <<  iterationCounterLog << "\t"
                      << residualErr << "\t"
                      << incrementalErr << "\t"
                      << etaLinesearch << std::endl;
            }
            performOneMoreStep = false;
            break;
          } else {
            // compute new deltaMatrices using a percentage of the previous one too
            LOG_TRACE(solvehyst) << "failBackCriterion not satisfied (yet); try to improve by including old deltaMatrices in computation";
            UInt percentage = 0; // 35;
            // does not work so well yet
            PDE_.SetFlagInCoefFncHyst("includeOldDeltaInPercent",percentage);            
          }
        }
      } else if((numResets != 0)&&allowReset_){
        // after first reset, check every Nth step for failback tolerance
        // N is the number of steps to reset from inputfile
        if(iterationCounter%minNumberItTillReset_ == 0){
          if(bestResNorm < failBackCrit_){
            LOG_TRACE(solvehyst) << "Best found solution so far passes failBack > take it";
            solVec_ = bestSol;
            residualErr = bestResNorm;
            // output of norms and data to info.xml
            if ( nonLinLogging_ == true ) {
              // get current step
              UInt actStep = PDE_.GetSolveStep()->GetActStep();
              
              //              if (PDE_.IsIterCoupled()) {
              //                WriteHystIterToInfoXML(pdename_,"Failback", couplingIter_, actStep,0,0, residualErrCheck, -1.0, 0);
              //              } else {
              //                //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
              //                WriteHystIterToInfoXML(pdename_,"Failback", actStep,0,0, residualErrCheck, -1.0, 0);
              //              }
              
              WriteHystIterToInfoXML(pdename_,
                      nonLinMethod_,
                      lineSearch_,
                      couplingIter_,
                      actStep,
                      iterationCounterLog,
                      iterationCounter,
                      residualErr, residualErrRel, useRelativeNormForRes_, 
                      incrementalErr, incrementalErrRel, useRelativeNormForInc_,
                      etaLinesearch, false, true);
              
              // write norm to file
              logFile_ <<  iterationCounterLog << "\t"
                      << residualErr << "\t"
                      << incrementalErr << "\t"
                      << etaLinesearch << std::endl;
            }
            performOneMoreStep = false;
            break;
          }
        }
      }
      //      
      //      LOG_DBG(solvehyst) << "Solution vectors at end of iteration: ";
      //      LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
      //      LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
      //      LOG_DBG(solvehyst) << "solVecForUpdate_: " << solVecLastTS_.ToString(); 
      //  
      //      std::cout << "#### SolveStepHyst::End of IT" << std::endl;
      
    } while(performOneMoreStep && (iterationCounter < nonLinMaxIter_));
    
    abortOnMaxIter_ = true;
    if (performOneMoreStep && (iterationCounter >= nonLinMaxIter_) && abortOnMaxIter_) {
      EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
              << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
              << "' at iteration '" << iterationCounter
              << "'.\n ==> incremental error (abs): " << incrementalErr
              << "'.\n ==> incremental error (rel): " << incrementalErrRel
              << "\n ==> residual error (abs): " << residualErr
              << "\n ==> residual error (rel): " << residualErrRel);
    }
    
    // check solution
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator limitFeFctIt;
    limitFeFctIt = feFunctions_.find(solutionLimit_);
    //std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    UInt pos;
    if (limitFeFctIt != feFunctions_.end() ) {
      for(pos = 0,fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt,++pos){
        if (fncIt == limitFeFctIt) { // pos is now referring to the corresponding subVec[pos]
          //const SingleVector * subv = solVec_.GetPointer(pos);
          Vector<Double> & dsubVec = dynamic_cast<Vector<Double> & > (*(solVec_.GetPointer(pos)));
          for (UInt j=0; j < dsubVec.GetSize(); j++) {
            if (dsubVec[j] >= maxValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
                      "' is larger than the allowed maximum limit set in the XML: "
                      << maxValidValue_);
            }
            if (dsubVec[j] <= minValidValue_) {
              EXCEPTION("A value ('" << dsubVec[j] << "') in the solution of PDE '" << pdename_ <<
                      "' is smaller than the allowed minimum limit set in the XML: "
                      << minValidValue_);
            }
          }
        }
      }
    }
    
    //update stage
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }
    
    // recover solution strategy
    nonLinMethod_ = nonLinBackup;
    towardsLastIteration_ = towardsLastIterationBackup;
    deltaMatCoupling_ = deltaMatCouplingBackup;
    deltaMatStrain_ = deltaMatStrainBackup;
    lineSearch_ = lineSearchBackup;
    
    // reset percentage of old deltaMat participation
    UInt percentage = 0;
    PDE_.SetFlagInCoefFncHyst("includeOldDeltaInPercent",percentage);  
    
    /*
     * solution has been found > set values AND overwrite memory
     * (-1 > same as +1 but additionally overwrite hyst memory
     */
    PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",-1);
    
    /*
     * STORE CURRENT VALUES FOR NEXT Time step
     */
    lastTS = true;
    // Note: by calling SetPreviousHystVals with flag lastTS = true,
    //        coefFncHyst will evaluate the hysteresis operator with
    //        the final state of solVec_; for this evaluation, the memory of
    //        the hysteresis operator will be unlocked, i.e. the evaluation will
    //        lead to irreversible changes
    //PDE_.SetPreviousHystVals(lastTS);
    // do this storing AFTER the final value was computed
    if(lastTS){
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);
    } else {
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
    }
    
    
    /*
     * Store IDBC values from the current time step in form of its rhs representation
     * (i.e. the effect that it will have on the rhs)
     * these values are needed to compute the deltaIDBC values
     * (currently not used)
     */
    algsys_->SetOldDirichletValues();
    
    /*
     * allow to output switching and rotation state as BMP images (if flag
     * printOut is set to 1 in material file)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",true);
    
    /*
     * set evaluationPurpose to 4, i.e. output
     * this will lock the hysteresis operator again and will only evaluate
     * the hystOperator at the midpoint of each element regardless of the
     * evaluation depth
     * > for more infos see coefFncHyst
     */
    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",4);
    //    
    //    LOG_DBG(solvehyst) << "Solution vectors at end of time step: ";
    //    LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
    //    LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
    // 
    //    std::cout << "#### SolveStepHyst::End of timestep" << std::endl;
  }
  
  void SolveStepHyst::WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const std::string& linesearchMethod,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const UInt iterationCounterTotal,
          const Double residualErr, const Double residualErrRel, bool useRelativeRelErr, 
          const Double incrementalErr, const Double incrementalErrRel, bool useRelativeIncErr,
          double etaLineSearch, bool reset, bool failback){
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter->Get("pde")->SetValue(pdeName);
    PtrParamNode iteration;
    if(coupledIterStep > 0){
      PtrParamNode coupledIteration = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT)->Get("coupledIt",ParamNode::APPEND);
      iteration = coupledIteration->GetByVal("outerIt","value",coupledIterStep,ParamNode::INSERT);
    } else {
      //      iteration = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT)->Get("it",ParamNode::APPEND);
      iteration = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT);
    }
    PtrParamNode info = iteration->GetByVal("it","value",iterationCounterTotal,ParamNode::INSERT)->Get("info",ParamNode::APPEND);
    if(reset){
      info->Get("reset",ParamNode::APPEND);
    } else if(failback){
      info->Get("failback",ParamNode::APPEND);
    } else {
      info->Get("method")->SetValue(nonLinSolveMethod);
      info->Get("linesearch")->SetValue(linesearchMethod);
      if(linesearchMethod != "none"){
        info->Get("etaLS")->SetValue(etaLineSearch);
      }
      PtrParamNode increment = iteration->GetByVal("it","value",iterationCounterTotal,ParamNode::INSERT)->Get("increment",ParamNode::APPEND);
      increment->Get("incErr")->SetValue(incrementalErr);
      increment->Get("incErrRel")->SetValue(incrementalErrRel);
      increment->Get("incErrCriterion")->SetValue(incStopCrit_);
      if(useRelativeIncErr){
        increment->Get("useRelativeErr")->SetValue("True");
      } else {
        increment->Get("useRelativeErr")->SetValue("False");
      }
      PtrParamNode residual = iteration->GetByVal("it","value",iterationCounterTotal,ParamNode::INSERT)->Get("residual",ParamNode::APPEND);
      residual->Get("resErr")->SetValue(residualErr);
      residual->Get("resErrRel")->SetValue(residualErrRel);
      residual->Get("resErrCriterion")->SetValue(residualStopCrit_);
      if(useRelativeIncErr){
        residual->Get("useRelativeErr")->SetValue("True");
      } else {
        residual->Get("useRelativeErr")->SetValue("False");
      }
    }
  }
  
  void SolveStepHyst::WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const UInt solStep,
          const UInt iterationCounter,
          const UInt iterationCounterTotal,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT)
            ->Get("it",ParamNode::APPEND);
    iter->Get("pde")->SetValue(pdeName);
    iter->Get("NLMethod")->SetValue(nonLinSolveMethod);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("nrTotal")->SetValue(iterationCounterTotal);
    iter->Get("resErr")->SetValue(residualErr);
    iter->Get("incErr")->SetValue(incrementalErr);
    // SE: include PDE name and time step
    if(etaLineSearch)
      iter->Get("etaLS")->SetValue(etaLineSearch);
  }
  
  void SolveStepHyst::WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const UInt iterationCounterTotal,
          const Double residualErr, 
          const Double incrementalErr, double etaLineSearch)
  {
    
    PtrParamNode iter = PDE_.GetInfoNode()->Get("nonlinearConvergence");
    iter = iter->GetByVal("solStep","value",solStep,ParamNode::INSERT);
    iter = iter->GetByVal("couplingStep", "value", coupledIterStep, ParamNode::INSERT)
            ->Get("it",ParamNode::APPEND);
    iter->Get("pde")->SetValue(pdeName);
    iter->Get("NLMethod")->SetValue(nonLinSolveMethod);
    iter->Get("nr")->SetValue(iterationCounter);
    iter->Get("nrTotal")->SetValue(iterationCounterTotal);
    iter->Get("resErr")->SetValue(residualErr);
    iter->Get("incErr")->SetValue(incrementalErr);
    // SE: include PDE name and time step
    if(etaLineSearch)
      iter->Get("etaLS")->SetValue(etaLineSearch);
  }
  
  
} // end of namespace
