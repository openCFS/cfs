#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>

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
  }
	
  //! Destructor
  SolveStepHyst::~SolveStepHyst() {
    //logFile_.close();
  }
  
	void SolveStepHyst::ReadNonLinDataHyst(){

		PtrParamNode nonLinNode = solStrat_->GetNonLinNode();
    
    // Check, if any nonlinear node was found
    Integer deltaForm;
    if( !nonLinNode ) {
      WARN("Taking default parameters for nonlinear data" );
      nonLinMethod_ = "HYST_deltaMat_TS";
      deltaForm = -1;
      lineSearch_ = "None";
      evalDepth_ = 2;
      measurePerformance_ = false;
      testInvesion_ = false;
      incStopCrit_ = 1e-4;
      residualStopCrit_ = 1e-4;
      nonLinMaxIter_ = 35;
    }
    
		if(nonLinNode){
      /*
       * method and linesearch already get read in base class
       * > here we only need to read in hyst specific paramater
       */
      if(nonLinMethod_ == "fixPoint"){
        std::stringstream except;
        except << "Nonlinear solution method fixPoint not directly implemented for Hysteresis.";
        except << "Please select HYST_fixPoint_TS or HYST_fixPoint_IT instead.";
        EXCEPTION(except.str())    
      } else if(nonLinMethod_ == "newton"){
        std::stringstream except;
        except << "Newton method not applicable to Hysteresis.";
        except << "Try HYST_deltaMat_TS or HYST_deltaMat_IT instead.";
        EXCEPTION(except.str())    
      }
     
      if( nonLinNode->Has("HYST_evalDepth") ) {
        nonLinNode->GetValue( "HYST_evalDepth", evalDepth_, ParamNode::INSERT );
      } else {
        evalDepth_ = 2; // extended
      }
      
      if( nonLinNode->Has("HYST_measurePerformance") ) {
        nonLinNode->GetValue( "HYST_measurePerformance", measurePerformance_, ParamNode::INSERT );
      } else {
        measurePerformance_ = 0;
      }
      
      if( nonLinNode->Has("HYST_testInversion") ) {
        nonLinNode->GetValue( "HYST_testInversion", testInvesion_, ParamNode::INSERT );
      } else {
        testInvesion_ = 0;
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
      PDE_.SetFlagInCoefFncHyst("testInversion",testInvesion_);
		}
    
    LOG_DBG(solvehyst) << "The following parameter were retrieved from .xml file:";
    LOG_DBG(solvehyst) << "nonLinMethod_: " << nonLinMethod_ << "( ==> deltaForm = " << deltaForm << ")";
    LOG_DBG(solvehyst) << "lineSearch_: " << lineSearch_;
    LOG_DBG(solvehyst) << "evalDepth_: " << evalDepth_;
    LOG_DBG(solvehyst) << "measurePerformance_: " << measurePerformance_;
    LOG_DBG(solvehyst) << "testInvesion_: " << testInvesion_;
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
		
		RhsVec_.Resize(feFunctions_.size());
		ResVec_.Resize(feFunctions_.size());
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
				
			RhsVec_.SetSubVector(new Vector<Double>(),id);
      RhsVec_.GetPointer(id)->Resize(rhsSize);
			
			ResVec_.SetSubVector(new Vector<Double>(),id);
      ResVec_.GetPointer(id)->Resize(rhsSize);
			
			stageRHSUpdate_.SetSubVector(new Vector<Double>(),id);
      stageRHSUpdate_.GetPointer(id)->Resize(rhsSize);
    }
    
		LinRhsVec_.Init();
    NonLinRhsVec_.Init();
		predictorCorrectorUpdate_.Init();
		
		RhsVec_.Init();
		ResVec_.Init();
		stageRHSUpdate_.Init();
    
    /*
     * Time levelflags
     * (-10 is no valid option for later usage and will trigger an initial setup
     *  in each case)
     */
    timeLevelMaterialCurrent_ = -10;
	  timeLevelDeltaMatCurrent_ = -10;
    
    
  }
	
  void SolveStepHyst::SolveStepTrans() {
    
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
      if( fncIt->second->GetTimeScheme()->isIncremental() == false){
        WARN("SolveStepHyst: currently only incremental stepping implemented");
        fncIt->second->GetTimeScheme()->forceIncremental();
      }
    } 		
   
    StepTransNonLinHysteresis();
  }  
  
	void SolveStepHyst::SetupRES(SBM_Vector& solutionForUpdate){
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
		Integer timeLevelDeltaMat; 
		/*
		 * timeLevelRHSHyst
		 * -2: hyst operator will not be added to rhs
		 * -1: output of hyst operator from last ts (n) will be added to rhs
		 *  0: current output of hyst operator (k) will be added to rhs
		 *  1: output of hyst operator from last it (k-1) will be added to rhs
		 */
		Integer timeLevelRHSHyst;
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
		timeLevelDeltaMat = -2;
		timeLevelRHSHyst = 0;
		
		SetupRESorRHS(timeLevelMaterial, timeLevelDeltaMat, timeLevelRHSHyst, solutionForUpdate, ResVec_);
	}
	
	void SolveStepHyst::SetupRHS(SBM_Vector& solutionForUpdate){
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
		Integer timeLevelDeltaMat; 
		/*
		 * timeLevelRHSHyst
		 * -2: hyst operator will not be added to rhs
		 * -1: output of hyst operator from last ts (n) will be added to rhs
		 *  0: current output of hyst operator (k) will be added to rhs
		 *  1: output of hyst operator from last it (k-1) will be added to rhs
		 */
		Integer timeLevelRHSHyst;
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
			timeLevelDeltaMat = -2;
			timeLevelRHSHyst = -2;
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
			timeLevelDeltaMat = -2;
			timeLevelRHSHyst = 0; 
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
			timeLevelDeltaMat = -2;
			timeLevelRHSHyst = -1;
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
			timeLevelDeltaMat = -2;
			timeLevelRHSHyst = 1;
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
			timeLevelDeltaMat = -2;
			timeLevelRHSHyst = -1;
		} else {
      std::stringstream except;
      except << "Nonlinear method " << nonLinMethod_ << " not implemented for hysteresis.";
			EXCEPTION(except.str())
		}
		
		SetupRESorRHS(timeLevelMaterial, timeLevelDeltaMat, timeLevelRHSHyst, solutionForUpdate, RhsVec_);
		algsys_->InitRHS(RhsVec_);
	}
	
  void SolveStepHyst::AssembleSystem(){
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
    
    if(nonLinMethod_ == "HYST_debug"){
			/*
			 * 0) debug mode = fixpoint iteration without considering hysteresis at all;
       *    > linear system has to be assembled
			 */
			timeLevelMaterial = -2;
			timeLevelDeltaMat = -2;
		} else if( (nonLinMethod_ == "HYST_fixPoint_IT")||(nonLinMethod_ == "HYST_fixPoint_TS") ){
			/*
			 * 1,2) fixpoint iteration towards last iteration k or last ts n:
       *    > nonlinear system (mat tensors depend on hysteresis) but no delta
			 */
			timeLevelMaterial = 0;
			timeLevelDeltaMat = -2;
		} else if(nonLinMethod_ == "HYST_deltaMat_IT"){
			/*
			 * 3) deltaMat computation towards last iteration k:
			 *			deltaU_k+1 = u_k+1 - u_k
       *    > delta mat to be evaluated with respect to last iteration (k-1)
			 */
			timeLevelMaterial = 0;
			timeLevelDeltaMat = 1;
		} else if(nonLinMethod_ == "HYST_deltaMat_TS"){
			/*
			 * 4) deltaMat computation towards last timestep n:
			 *			deltaU_k+1 = u_k+1 - u^n
			 *    > delta mat to be evaluated with respect to last timestep (n)
			 */
			timeLevelMaterial = 0;
			timeLevelDeltaMat = -1;
		} else {
      std::stringstream except;
      except << "Nonlinear method " << nonLinMethod_ << " not implemented for hysteresis.";
			EXCEPTION(except.str())
		}
    
    AssembleSystem(timeLevelMaterial, timeLevelDeltaMat);
    
  }
  
  void SolveStepHyst::AssembleSystem(Integer timeLevelMaterial, Integer timeLevelDeltaMat){
    /*
     *  Assemble system according to flags timeLevelDeltaMat, timeLevelMaterial
     *  > this is used for res and rhs assembly as here the system might be different
     *    from the one that actually has to be solved
     */
    SBM_Vector diff;
    diff = lastUpdateVecForMatrices_;
    diff.Add(-1.0,solVec_);
    bool systemRequiresReassembly = (diff.NormL2() < 1e-16);
    
    if( (timeLevelDeltaMat != timeLevelDeltaMatCurrent_) ||
				(timeLevelMaterial != timeLevelMaterialCurrent_) ||
				(systemRequiresReassembly == true) ) {
			LOG_TRACE(solvehyst) << "SolveStepHyst::Assemble matrices";
      
      // evalPurpose 1 > assemble
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
			PDE_.SetFlagInCoefFncHyst("SetTimeLevelMaterial",timeLevelMaterial);
			PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMat",timeLevelDeltaMat);
      // important: allow reevaluation of hyst operator > might have changed
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);
			assemble_->AssembleMatrices(false);

			timeLevelMaterialCurrent_ = timeLevelMaterial;
			timeLevelDeltaMatCurrent_ = timeLevelDeltaMat;
      lastUpdateVecForMatrices_ = solVec_;
		}
  }
  
	void SolveStepHyst::SetupRESorRHS(Integer timeLevelMaterial, Integer timeLevelDeltaMat, 
		Integer timeLevelRHSHyst, SBM_Vector& solutionForUpdate, SBM_Vector& returnVector){
		
    LOG_DBG(solvehyst) << "SetupRESorRHS";
    LOG_DBG(solvehyst) << "timeLevelMaterial (-2 = deactivated; -1 = lastTS; 0 = curIT; +1 = lastIT): " << timeLevelMaterial; 
    LOG_DBG(solvehyst) << "timeLevelDeltaMat (-2 = deactivated; -1 = lastTS; 0 = INVALID; +1 = lastIT): " << timeLevelDeltaMat; 
    LOG_DBG(solvehyst) << "solutionForUpdate: " << solutionForUpdate.ToString(); 
    
		/*
		 * 0. Init returnVector
		 */
		returnVector.Init();
		
		/*
		 * 1. Assemble system if required
     * > here we have to use the version with prescribed timelevels as those
     *   might differ from the ones that are required by nonLinMethod
		 */
    AssembleSystem(timeLevelMaterial, timeLevelDeltaMat);
    
		/*
		 * 2. Setup Linear part and predictor corrector update
		 */
		if(LinRHSRequiresAssembly_){
			LOG_DBG(solvehyst) << "SolveStepHyst::Assemble linear RHS/RES";
		
			LinRhsVec_.Init();
      algsys_->InitRHS(LinRhsVec_);
      assemble_->AssembleLinRHS();
      algsys_->GetRHSVal( LinRhsVec_ );
						
			LinRHSRequiresAssembly_ = false;
		}
		returnVector.Add(1.0,LinRhsVec_);
		
		if(PredictorCorrectorRequiresAssembly_){
			LOG_DBG(solvehyst) << "SolveStepHyst::Compute predictor/corrector update";
		
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
		
		/*
		 * 3. Setup nonlin RHS with hysteresis
     * > nonlin RHS is solution dependent 
     * > if current solution differes from previously used one >  reassemble
		 */
    SBM_Vector diff;
    diff = lastUpdateVecForNonLinRHS_;
    diff.Add(-1.0,solVec_);
    bool NonLinRHSRequiresReassembly = (diff.NormL2() < 1e-16);
    
		if( (timeLevelRHSHyst != timeLevelRHSHystCurrent_) ||
				(NonLinRHSRequiresReassembly == true) ) {
			LOG_DBG(solvehyst) << "SolveStepHyst::Assemble non-linear RHS/RES";

			NonLinRhsVec_.Init();
			algsys_->InitRHS(NonLinRhsVec_); //load storage into algsys
      // evalPurpose 1 > assemble
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
			PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSHyst",timeLevelRHSHyst);	
      // important: allow reevaluation of hyst operator > might have changed
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);
			
      assemble_->AssembleNonLinRHS(); 
			algsys_->GetRHSVal( NonLinRhsVec_ );

			timeLevelRHSHystCurrent_ = timeLevelRHSHyst;
      lastUpdateVecForNonLinRHS_ = solVec_;
		}
		returnVector.Add(1.0,NonLinRhsVec_);
		
		/*
		 * 4. Update res/rhs according to update solution vector
		 */
		// load returnVector into algsys; then update for incremental formulation
		algsys_->InitRHS(returnVector);

		// Stiffness part
		LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u";
		LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 1: K";
		solutionForUpdate.ScalarMult(-1.0);
		algsys_->UpdateRHS(STIFFNESS,solutionForUpdate,true);
		if(!trans_){
			algsys_->UpdateRHS(STIFFNESS_UPDATE,solutionForUpdate,true);
		}
		solutionForUpdate.ScalarMult(-1.0);

		// Damping, Mass part (with corresponding scaling > done via timestepping)
		if(trans_){
			LOG_DBG(solvehyst) << "SolveStepHyst::Update with K_eff*u - Part 2: C,M";
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
																						solutionForUpdate.GetPointer(fctId), -1.0, forceReset);

				}
				// now update loaded rhs with matrix*scaled solution vector
				algsys_->UpdateRHS(matIt->first,stageRHSUpdate_,true);
			}
		}
		/*
		 * 5. Retrieve vector
		 */
		algsys_->GetRHSVal( returnVector );
		LOG_DBG(solvehyst) << "returnVector: " << returnVector.ToString(); 
	}
	  
  // TODO: UPDATE / UPGRADE for new implementation
  Double SolveStepHyst::LineSearchHyst(SBM_Vector& solIncrement, SBM_Vector& stageSol, Double& etaLineSearch, UInt evalVersion,
          UInt callingCnt, bool trans, bool performLineSearch,
          bool forceReevaluation, bool debugOutput, bool reset, UInt allowedSteps)  {
    
    static Timer* lineSearchTimer = new Timer();
    static UInt numCalls = 0;
    if(debugOutput){
      std::cout << "LineSearchHyst" << std::endl;
      lineSearchTimer->Start();
      numCalls++;
    }
    
    /*
     * backup nonlinear rhs, residual and solution vector
     */
    SBM_Vector solOld(BaseMatrix::DOUBLE);
    solOld = solVec_;
    
    SBM_Vector resOld(BaseMatrix::DOUBLE);
    resOld = currentResVec_;
    
    
    const UInt nrEtas = 13;
    //Double eta[nrEtas] = {1.0,0.5,0.2,0.1};
    Double eta[nrEtas] = {1.0,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.0025,0.001,0.0005,0.00025,0.0001};
    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01,0.00316,0.001,0.000316,0.0001,0.0000316,0.00001,0.0000316,0.00001};
    //Double eta[nrEtas] = {1.0,0.075,0.05,0.025,0.1};
    //Double eta[nrEtas] = {1.0,0.316,0.1,0.0316,0.01};
    //Double eta[nrEtas] = {1.0,0.75,0.5,0.25,0.1};
    bool assumeQuadratic = false;
    
    // initialize etaOpt or receive compiler warning
    Double etaOpt = 0.0;
    Double residualL2NormOpt = 1e15;
    
    if(performLineSearch == false){
      /*
       * only take first eta = 1
       * -> we do not have to lock hysteresis nor do we have to backup
       *    the solution or the residual as we directly take the first value
       */
      etaOpt = 1.0;
    } else {
      /*
       * lock hysteresis memory
       * NOTE: in the current evaluation, the locking is done automatically
       * in coefFnc hyst if the system is called during assembly;
       * only if SetPreviousHystValues is used with setLastTS = true, the
       * hysteresis memory is unlocked
       */
      
      /*
       * test different eta
       */
     // bool skipReassembly = true;
      UInt nrEtasToTest = std::min(allowedSteps,nrEtas);
      for( UInt i=0; i<nrEtasToTest; i++) {
        if(debugOutput){
          std::cout << "testing eta: " << eta[i] << std::endl;
        }
        
        //CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,eta[i], 0, callingCnt, evalVersion,
        //        trans, forceReevaluation, skipReassembly, debugOutput, reset);
        
        /*
         * calculate error norm as criterion
         */
        Double residualL2Norm = currentResVec_.NormL2();
        
        if(debugOutput){
          std::cout << "L2-norm of computed residual for eta = " << eta[i] << ": " << residualL2Norm << std::endl;
        }
        
        /*
         * restore old solution vector, old residual and nonlinear rhs
         * -> needed to perform a correct next call to CalcResidual...
         */
        solVec_ = solOld;
        currentResVec_ = resOld;
        
        if (residualL2Norm < residualL2NormOpt) {
          residualL2NormOpt = residualL2Norm;
          etaOpt = eta[i];
        } else {
          if(assumeQuadratic){
            /*
             * here we assume that the residual follows a quadratic curve
             * i.e. we only have one minimum which is the global minimum;
             * in that case, if the residual increases between two iterations,
             * we can assume that we passed the minimum already and stop
             */
            break;
          }
        }
      }
      
      /*
       * unlock hysteresis memory
       * -> no longer done here but first when we store the lastTS values
       */
      
    }
    
    /*
     * now we have found etaOpt, so we can call CalcResidualAndConfigSystemForHysteresis
     * one last time to setup the correct system for next computation
     */
    etaLineSearch = etaOpt;
    
    if(debugOutput){
      lineSearchTimer->Stop();
      std::cout << "Total time needed for linesearch: " << lineSearchTimer->GetCPUTime() << std::endl;
      std::cout << "Average time needed for linesearch: " << lineSearchTimer->GetCPUTime()/numCalls << std::endl;
      std::cout << "(optimal) eta: " << etaOpt << std::endl;
    }
    
    // calculate res and assemble system once more using the optimal eta
   // bool skipReassembly = false;
   // CalcResidualAndConfigSystemForHysteresis(solOld,solIncrement,stageSol,etaLineSearch, 0, callingCnt, evalVersion,
    //        trans, forceReevaluation, skipReassembly, debugOutput, reset);
    
    if(debugOutput){
      std::cout << "minimal residual: " << residualL2NormOpt << std::endl;
    }
    
    return currentResVec_.NormL2();
  }
    
  void SolveStepHyst::StepTransNonLinHysteresis() {
		static UInt timestepCnt = 0;
		timestepCnt++;
		LOG_TRACE(solvehyst) << "#### SolveStepHyst::Begin of timestep " << timestepCnt;
		LOG_DBG(solvehyst) << "SolVec_ at start of timestep: " << solVec_.ToString(); 
    /*
     * Initial checks
     */
		if(!flagsInitialized_){
			ReadNonLinDataHyst();
			flagsInitialized_ = true;	
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
    SBM_Vector actSol(BaseMatrix::DOUBLE);
    actSol = solVec_;
    
    // solVec_ currently holds the solution from the previous timestep
    // create backup for later usage
    solVecLastIT_ = solVec_;
    solVecLastTS_ = solVec_;
    lastUpdateVecForNonLinRHS_ = solVecLastTS_;
    lastUpdateVecForMatrices_ = solVecLastTS_;
    LOG_DBG(solvehyst) << "lastUpdateVecForNonLinRHS_: " << lastUpdateVecForNonLinRHS_.ToString() ;
    LOG_DBG(solvehyst) << "lastUpdateVecForMatrices_: " << lastUpdateVecForMatrices_.ToString() ;
    
    /*
     * Compute norm of linear rhs as reference value for relative res. error
     */
    RhsVec_.Init();
    algsys_->InitRHS(RhsVec_);
    
    // setup right hand side
    Double loadFactor = 1.0;
    Double RhsLinL2Norm = SetLinRHS(loadFactor);
    
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
    stageSol = actSol;
    solVec_  = actSol;
    
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
    bool performOneMoreStep;
    Double solL2Norm;
    Double residualErr;
    Double incrementalErr;
    
    // keep track of residual; if it is not decreasing over several iterations > reset (see below)
    UInt numResets = 0;
    std::deque<Double> resNormHistory(5, 50000);
    
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
    
    // DEBUG -- FORCE delta last ts
    nonLinMethod_ = "HYST_deltaMat_TS";
    
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
    bool lastTS;
    std::string nonLinBackup = nonLinMethod_;
    
    /*
     *  Begin iterations
     */
    do {
      /*
       * Step 0: Increase counter 
       */
      overallIterationCounter++;
      iterationCounter++;
      LOG_TRACE(solvehyst) << "## Iteration: " << iterationCounter ;

      /*
       * Step 0.2: Reset solution UPDATE vector (if necessary)
       * (this vector is the one that gets subtracted from the rhs, i.e. the
       * one the stepping is going to)
       * > NOTE: do not overwrite solVec_ here, as the system still needs to be
       *         assembled using the current solution!
       * 
       */
      if( (steppingTowardsLastTS)&&(!steppingTowardsZero) ){
        // increment has to be towards last ts > solVecForUpdate_ has to hold this
        // value at begin of iteration
        solVecForUpdate_ = solVecLastTS_;
      } else if( (steppingTowardsLastTS)&&(steppingTowardsZero) ){
        // increment has to be towards 0 > solVecForUpdate_ has to be reset to 0
        solVecForUpdate_.Init();
      } else {
        // stepping towards last iteration > solVecForUpdate_ = solVec_
        solVecForUpdate_ = solVecLastIT_;
      }
      
      LOG_DBG(solvehyst) << "Solution vectors at start of iteration: ";
      LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
      LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
      LOG_DBG(solvehyst) << "solVecLastIT_: " << solVecLastIT_.ToString(); 
      LOG_DBG(solvehyst) << "solVecForUpdate_: " << solVecForUpdate_.ToString(); 
      
      LOG_DBG(solvehyst) << "Try following nonlin method: " << nonLinMethod_;
      LOG_DBG(solvehyst) << "Flags: ";
      LOG_DBG(solvehyst) << "steppingTowardsLastTS: " << steppingTowardsLastTS;
      LOG_DBG(solvehyst) << "steppingTowardsZero: " << steppingTowardsZero;

      /*
       * Step 1: setup rhs
       */
      LOG_DBG(solvehyst) << "SetupRHS";
      SetupRHS(solVecForUpdate_);
      
      LOG_DBG(solvehyst) << "RHS to solve for: : " << RhsVec_.ToString(); 
      
      /*
       * Step 2: (re-)assemble matrices according to nonLinMethod
       * Note: During SetupRHS() the matrices already get assembled, but 
       *        most probably with different flags > reassembly needed
       */
      LOG_DBG(solvehyst) << "AssembleSystem";
      AssembleSystem();
      
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
       */
      solVecLastIT_ = solVec_;
            
      lastTS = false;
      // Note: by calling SetPreviousHystVals with flag lastTS = false,
      //        coefFncHyst will evaluate the hysteresis operator with
      //        the current state of solVec_; for this evaluation, the memory of
      //        the hysteresis operator will be locked, i.e. the evaluation will
      //        lead to reversible changes
      PDE_.SetPreviousHystVals(lastTS);
      
      /*
       * Step 4: solve system and retrieve the solution (-update)
       */
      LOG_DBG(solvehyst) << "Solve";
      algsys_->Solve(setIDBC,deltaIDBC);
      /*
       *  RETRIEVE SOLUTION INCREMENT
       *  - nonLinMethod_ = "HYST_debug","HYST_fixPoint_IT","HYST_deltaMat_IT" > increment towards last iteration
       *  - nonLinMethod_ = "HYST_fixPoint_TS","HYST_deltaMat_TS" > increment towards last ts
       *  - steppingTowardsZero > increment towards 0 > full step
       */
      algsys_->GetSolutionVal(solInc, setIDBC, deltaIDBC );
      
      LOG_DBG(solvehyst) << "solInc (computed): " << solInc.ToString();
      
      /*
       * The solution vector will be the increment towards solVecForUpdate_,
       * i.e. the new solVec will be solVecForUpdate + eta*solInc
       */
      solVec_ = solVecForUpdate_;
            
      /*
       * Step 5: check solution by evaluating residual and increment
       */
      Double residualL2Norm = 0.0;
      Double incrementL2Norm = 0.0;
      Double etaLinesearch = 1.0;
      
      LOG_DBG(solvehyst) << "SetupRES/Perform LS";
      if( setIDBC || ( lineSearch_ == "none" ) ){
        /*
         * Linesearch is only allowed if the increment does not hold contributions
         * of the Dirichlet boundary as in that case, we have to step with length 1
         * to get the boundary terms correct!
         * > during first iteration for incremental towards last iteration
         * > during each iteration for incremental towards last ts or 0
         */
        solVec_.Add(etaLinesearch,solInc);
        LOG_DBG(solvehyst) << "Updated solution: " << solVec_.ToString();
        // different to setupRHS, the residual always uses solVec_
        SetupRES(solVec_);
      } else {
        // TODO: implement linsearch
        solVec_.Add(etaLinesearch,solInc);
        LOG_DBG(solvehyst) << "Updated solution: " << solVec_.ToString();
        SetupRES(solVec_);
      }
      
      // update actSol to currently computed value
      actSol = solVec_;
      
      residualL2Norm = ResVec_.NormL2();
      if ( RhsLinL2Norm > 1.0 ){
        residualErr = residualL2Norm / RhsLinL2Norm;
      } else {
        residualErr = residualL2Norm;
      }
      resNormHistory.push_back(residualErr);
      resNormHistory.pop_front();
      
      /*
       * For incremental error checking, compute the actual update between
       *  the current and the previous iteration
       */
      solInc = solVec_;
      solInc.Add(-1.0,solVecLastIT_);
      
      LOG_DBG(solvehyst) << "solInc (towards lastIT): " << solInc.ToString();
      
      incrementL2Norm = solInc.NormL2();
      solL2Norm  = solVec_.NormL2();
      
      if ( solL2Norm > 1.0 ){
        incrementalErr = incrementL2Norm / solL2Norm;
      } else {
        incrementalErr = incrementL2Norm;
      }
      
      performOneMoreStep =
              (incrementalErr > incStopCrit_) || (residualErr > residualStopCrit_);
      
      Double minDecrease = 1e-5;
      // compare residual at front and back of deque; 
      // residual should have decreased (or not enough); otherwise try reset
      if( (performOneMoreStep) && (resNormHistory.front()-minDecrease <= resNormHistory.back()) ){
        solVec_.Init();
        if(numResets == 0){
          LOG_TRACE(solvehyst) << "Residual did not decrease using solution method " << nonLinMethod_;
          LOG_TRACE(solvehyst) << "Reset solVec to 0 and try again."; 
        } else if(numResets == 1){
          LOG_TRACE(solvehyst) << "Residual did not decrease using solution method " << nonLinMethod_;
          LOG_TRACE(solvehyst) << "Reset solVec to 0 and try solution method HYST_deltaMat_IT"; 
          nonLinMethod_ = "HYST_deltaMat_IT";
          steppingTowardsLastTS = false;
        } else {
          EXCEPTION("Error could not be decreases by resetting");
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
        bool forceMemoryLock = true;
        lastTS = steppingTowardsLastTS;
        PDE_.SetPreviousHystVals(lastTS,forceMemoryLock);
        
        LOG_DBG(solvehyst) << "Iteration " << iterationCounter << " of timestep " << timestepCnt;
        
        // reset deque
        for(UInt i = 0; i < resNormHistory.size(); i++){
          resNormHistory[i] = 50000;
        }
        steppingTowardsZero = true;
        numResets++;
      }
            
      // output of norms and data to info.xml
      if ( nonLinLogging_ == true ) {
        // get current step
        UInt actStep = PDE_.GetSolveStep()->GetActStep();
        
        if (PDE_.IsIterCoupled()) {
          WriteNonLinIterToInfoXML(pdename_, couplingIter_, actStep,iterationCounter, residualErr, incrementalErr, etaLinesearch);
        } else {
          //WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErrABS, etaLineSearch);
          WriteNonLinIterToInfoXML(pdename_, actStep,iterationCounter, residualErr, incrementalErr, etaLinesearch);
        }
        // write norm to file
        logFile_ <<  iterationCounter << "\t"
                << residualErr << "\t"
                << incrementalErr << "\t"
                << etaLinesearch << std::endl;
      }
      
      LOG_DBG(solvehyst) << "Solution vectors at end of iteration: ";
      LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
      LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
      LOG_DBG(solvehyst) << "solVecForUpdate_: " << solVecLastTS_.ToString(); 
      
    } while(performOneMoreStep && (iterationCounter < nonLinMaxIter_));
    
    abortOnMaxIter_ = true;
    if (performOneMoreStep && (iterationCounter >= nonLinMaxIter_) && abortOnMaxIter_) {
      EXCEPTION("NON CONVERGENCE error in PDE '" << pdename_
              << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
              << "' at iteration '" << iterationCounter
              << "'.\n ==> incremental error: " << incrementalErr
              << "\n ==> residual error: " << residualErr);
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
    
    // reconver solution strategy
    nonLinMethod_ = nonLinBackup;
    
    /*
     * STORE CURRENT VALUES FOR NEXT Time step
     */
    lastTS = true;
    // Note: by calling SetPreviousHystVals with flag lastTS = true,
    //        coefFncHyst will evaluate the hysteresis operator with
    //        the final state of solVec_; for this evaluation, the memory of
    //        the hysteresis operator will be unlocked, i.e. the evaluation will
    //        lead to irreversible changes
    PDE_.SetPreviousHystVals(lastTS);
    
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
    
    LOG_DBG(solvehyst) << "Solution vectors at end of time step: ";
    LOG_DBG(solvehyst) << "solVec_: " << solVec_.ToString(); 
    LOG_DBG(solvehyst) << "solVecLastTS_: " << solVecLastTS_.ToString(); 
    
  }
 
} // end of namespace
