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
#include "MatVec/SBM_Matrix.hh"

//#include "MatVec/SingleVector.hh"

namespace CoupledField {

  /*
   * Logging system
   *
   * LOG_TRACE = Calls to / steps inside main functions / most important results
   * LOG_TRACE2 = Steps inside sub functions
   * LOG_DBG = Global debugging info like convergence results
   * LOG_DBG2 = Specialized debugging stuff like state of vectors
   * LOG_DBG3 = Very specialized debugging stuff e.g. for checking input parameter
   */
  DEFINE_LOG(solvestephyst_main, "solvestephyst_main")

  // important subfunctions get own logger
  DEFINE_LOG(solvestephyst_bgm, "solvestephyst_bgm")
  DEFINE_LOG(solvestephyst_linesearch, "solvestephyst_linesearch")
  DEFINE_LOG(solvestephyst_residual, "solvestephyst_residual")
  DEFINE_LOG(solvestephyst_assemble_and_build, "solvestephyst_assemble_and_build")
  DEFINE_LOG(solvestephyst_helperfuncs, "solvestephyst_helperfuncs")

  /*
   * I. Constructor and deconstructor
   */
  SolveStepHyst::SolveStepHyst(StdPDE & apde)
  :StdSolveStep(apde)
  {
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== SolveStepHyst                 ==";
    LOG_TRACE(solvestephyst_main) << "===================================";

		trans_ = true;
    materialTensorsHystDepenedent_ = false;
    nonHystereticNonLinearitiesPresent_ = false;
    timestepCnt_ = 0;
    totalIterationCnt_ = 0; // sum over all iterations from all time steps
    maxIterationPerTSCnt_ = 0;
    initStep1Done_ = false;

    /*
     * Time levelflags
     * (-10 is no valid option for later usage and will trigger an initial setup
     *  in each case)
     */
//    timeLevelMaterialCurrent_ = -10;
//    timeLevelDeltaMatPolCurrent_ = -10;
//    timeLevelDeltaMatStrainCurrent_ = -10;
//    timeLevelDeltaMatCouplingCurrent_ = -10;

//    timeLevelRHSPolCurrent_ = -10;
//    timeLevelRHSStrainCurrent_ = -10;
//    timeLevelRHSCouplingCurrent_ = -10;
	}

  //! Destructor
  SolveStepHyst::~SolveStepHyst() {
    // TODO: free memory of arrays
    delete[] resIncrements_history_;
    delete[] solIncrements_history_;
  }

  /*
   * II. Init
   * > these functions are called only ONCE per sequence step
   */
  void SolveStepHyst::InitTimeStepping(){
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== InitTimeStepping              ==";
    LOG_TRACE(solvestephyst_main) << "===================================";

    // call base class first
		StdSolveStep::InitTimeStepping();

    /*
     * Read in defining parameter for solution procedure
     * from input file; needs to be done only once
     */
    if(initStep1Done_ == false){
      Init_Step1_ReadNonLinDataHyst();
    }

    materialTensorsHystDepenedent_ = PDE_.MaterialTensorsHystDependent();
		nonHystereticNonLinearitiesPresent_ = PDE_.IsNonLinAndNonHyst();
    if(nonHystereticNonLinearitiesPresent_){
      WARN("SolveStepHyst::SolveStepHyst --- The combination of hysteretic and non-hysteretic nonlinearities has not been tested yet.")
    }

    // required for each region > has to be performed multiple times
    Init_Step2_InitGlobalArrays();
  }

  LinesearchParameter SolveStepHyst::Init_Step1b_ReadSingleLinesearch(PtrParamNode& singleLinesearchNode,std::string nameTag){

    LinesearchParameter currentLinesearchParameter = LinesearchParameter();
    currentLinesearchParameter.nameTag_ = nameTag;
    bool allowNegEta = false;

    /*
     * Read in specific parameter
     */
    if(nameTag == "Backtracking_SmallestResidual"){
      singleLinesearchNode->GetValue("etaStart",currentLinesearchParameter.etaStart_,ParamNode::PASS);
      singleLinesearchNode->GetValue("etaMin",currentLinesearchParameter.minEta_,ParamNode::PASS);
      singleLinesearchNode->GetValue("stoppingTol",currentLinesearchParameter.stoppingTol_,ParamNode::PASS);
      singleLinesearchNode->GetValue("decreaseFactor",currentLinesearchParameter.decreaseFactor_,ParamNode::PASS);
      singleLinesearchNode->GetValue("residualDecreaseForSuccess",currentLinesearchParameter.residualDecreaseForSuccess_,ParamNode::PASS);

      if(singleLinesearchNode->Has("testOverrelaxation")){
        singleLinesearchNode->GetValue("testOverrelaxation",currentLinesearchParameter.testOverrelaxation_,ParamNode::PASS);
      } else {
        currentLinesearchParameter.testOverrelaxation_ = false;
      }
    }

    if( (nameTag == "Backtracking_Armijo") ||
        (nameTag == "Backtracking_NonMonotonicArmijo") ||
        (nameTag == "Backtracking_GradientFreeArmijo") ){

      singleLinesearchNode->GetValue("etaStart",currentLinesearchParameter.etaStart_,ParamNode::PASS);
      singleLinesearchNode->GetValue("etaMin",currentLinesearchParameter.minEta_,ParamNode::PASS);
      singleLinesearchNode->GetValue("decreaseFactor",currentLinesearchParameter.decreaseFactor_,ParamNode::PASS);
      singleLinesearchNode->GetValue("rho",currentLinesearchParameter.rho_,ParamNode::PASS);

      // only for non-monotonic Armijo
      if(singleLinesearchNode->Has("historyLength")){
        singleLinesearchNode->GetValue("historyLength",currentLinesearchParameter.historyLengthArmijo_,ParamNode::PASS);
      } else {
        currentLinesearchParameter.historyLengthArmijo_ = 1;
      }
    }

    if( (nameTag == "Exact_GoldenSection") ||
        (nameTag == "Exact_QuadraticPolynomial") ){

      singleLinesearchNode->GetValue("stoppingTol",currentLinesearchParameter.stoppingTol_,ParamNode::PASS);

      // exact linesearches may require negative eta, i.e., set default to true
      allowNegEta = true;
    }

    if(nameTag == "SufficientDecrease_ThreePointPolynomial"){

      singleLinesearchNode->GetValue("alphaCheck",currentLinesearchParameter.alphaCheck_,ParamNode::PASS);
      singleLinesearchNode->GetValue("sigma0",currentLinesearchParameter.sigma0_,ParamNode::PASS);
      singleLinesearchNode->GetValue("sigma1",currentLinesearchParameter.sigma1_,ParamNode::PASS);
    }

    if(nameTag == "SufficientDecrease_VectorPolynomial"){

      singleLinesearchNode->GetValue("epsilon",currentLinesearchParameter.epsilon_,ParamNode::PASS);

      // negative stepping explicitly intended for this linesearch
      allowNegEta = true;
    }

    /*
     * Read in common parameter that appear in every linesearch
     */
    // mandatory
    singleLinesearchNode->GetValue("maxIterLS",currentLinesearchParameter.maxIter_,ParamNode::PASS);

    // optional
    if(singleLinesearchNode->Has("measurePerformance")){
      singleLinesearchNode->GetValue("measurePerformance",currentLinesearchParameter.measurePerformance_,ParamNode::PASS);
    } else {
      currentLinesearchParameter.measurePerformance_ = false;
    }

    if( (nameTag == "Exact_GoldenSection") ||
            (nameTag == "Exact_QuadraticPolynomial") ){
      // the exact linesearches require the initial determination of a searching interval; if we measure the performance
      // of these linesearches, we should measure the performance of this initial search, too
      DEBUG_measureLSPerformance_ = currentLinesearchParameter.measurePerformance_;
    }

    if(singleLinesearchNode->Has("allowNegativeSteps")){
      singleLinesearchNode->GetValue("allowNegativeSteps",currentLinesearchParameter.allowNegativeEta_,ParamNode::PASS);
    } else {
      currentLinesearchParameter.allowNegativeEta_ = allowNegEta;
    }

    return currentLinesearchParameter;
  }

  void SolveStepHyst::Init_Step1a_ReadSolutionMethod(PtrParamNode& methodNode){
    /*
     * read in SINGLE solution method for non-linear, hysteretic system
     * > only single method is allowed (currently; maybe one can add some kind of switch later)
     * > store parameter directly in class
     */

    /*
     * flag indicating selected solution approach
     */
    solutionApproach_ = SOLVE_NOTSET;
    deltaMatCoupling_ = false;
    deltaMatStrain_ = false;
    nonLinMethod_ = "None";
    initialNumberFPSteps_ = 0;

    /*
     * for fixpoint methods
     * global and local contraction factors
     */
    // fixpoint contraction factors
    FPLocalC_ = 2.0; // should be larger 1; for most problems a value between 1 and 2 is ok; higher value leads to
    // larger convergence intervals but to slower convergence rates
    FPGlobalC_ = 1.0;

    /*
     * for quasi-Newton methods
     */
    // 1: compute Jacobian every iteration
    // 2: compute Jacobian every second iteration
    // n: compute Jacobian every nth iteration
    numIterationsTillUpdateOfJacobian_ = 1;

    PtrParamNode currentMethod = NULL;
    bool evalJacAtMidpointOnly;
    /*
     * Read in fix-point
     */
    if( methodNode->Has("Fixpoint_Global_B")){
      currentMethod = methodNode->Get("Fixpoint_Global_B");
      nonLinMethod_ = "HYST_FP_GLOBAL";
      solutionApproach_ = SOLVE_GLOBAL_FIXPOINT_B;
      FPApproach_ = FP_GLOBAL_B;
      currentFPApproach_ = FP_GLOBAL_B;
      currentMethod->GetValue( "ContractionFactor", FPGlobalC_, ParamNode::PASS );
    }

    if( methodNode->Has("Fixpoint_Localized_B")){
      currentMethod = methodNode->Get("Fixpoint_Localized_B");
      nonLinMethod_ = "HYST_FP_LOCAL_v2";
      solutionApproach_ = SOLVE_LOCAL_FIXPOINT_B_v2;
      FPApproach_ = FP_LOCAL_B_DIAG;
      currentFPApproach_ = FP_LOCAL_B_DIAG; //FP_LOCAL_B_FULLJAC
      currentMethod->GetValue( "ContractionFactor", FPLocalC_, ParamNode::PASS );
      if(currentMethod->Has("initialNumberGlobalFPSteps")){
        currentMethod->GetValue( "initialNumberGlobalFPSteps", initialNumberFPSteps_, ParamNode::PASS );
      } else {
        initialNumberFPSteps_ = 0;
      }
      if(currentMethod->Has("estimateFPSlopeAtMidpointOnly")){
        currentMethod->GetValue( "estimateFPSlopeAtMidpointOnly", evalJacAtMidpointOnly, ParamNode::PASS );
      } else {
        evalJacAtMidpointOnly = false;
      }
    }

    if( methodNode->Has("Fixpoint_Global_H")){
      currentMethod = methodNode->Get("Fixpoint_Global_H");
      nonLinMethod_ = "HYST_FP_GLOBAL_H";
      solutionApproach_ = SOLVE_GLOBAL_FIXPOINT_H;
      FPApproach_ = FP_GLOBAL_H;
      currentFPApproach_ = FP_GLOBAL_H;
      currentMethod->GetValue( "ContractionFactor", FPGlobalC_, ParamNode::PASS );
    }

    if( methodNode->Has("Fixpoint_Localized_H")){
      currentMethod = methodNode->Get("Fixpoint_Localized_H");
      nonLinMethod_ = "HYST_FP_LOCAL_H";
      solutionApproach_ = SOLVE_LOCAL_FIXPOINT_H;
      FPApproach_ = FP_LOCAL_H_DIAG;
      currentFPApproach_ = FP_LOCAL_H_DIAG; //FP_LOCAL_H_FULLJAC
      currentMethod->GetValue( "ContractionFactor", FPLocalC_, ParamNode::PASS );
      if(currentMethod->Has("initialNumberGlobalFPSteps")){
        currentMethod->GetValue( "initialNumberGlobalFPSteps", initialNumberFPSteps_, ParamNode::PASS );
      } else {
        initialNumberFPSteps_ = 0;
      }
      if(currentMethod->Has("estimateFPSlopeAtMidpointOnly")){
        currentMethod->GetValue( "estimateFPSlopeAtMidpointOnly", evalJacAtMidpointOnly, ParamNode::PASS );
      } else {
        evalJacAtMidpointOnly = false;
      }
    }

    /*
     * Read in quasi-Newton
     */
    if( methodNode->Has("QuasiNewton_ChordMethod")){
      nonLinMethod_ = "HYST_FD_CHORD";
      solutionApproach_ = SOLVE_CHORD;
      FPApproach_ = NOFP;
      currentFPApproach_ = NOFP;
      currentMethod = methodNode->Get("QuasiNewton_ChordMethod");
      if(currentMethod->Has("initialNumberFPSteps")){
        currentMethod->GetValue( "initialNumberFPSteps", initialNumberFPSteps_, ParamNode::PASS );
      } else {
        initialNumberFPSteps_ = 1;
      }
      
      if(currentMethod->Has("calculateFDJacobianAtMidpointOnly")){
        currentMethod->GetValue( "calculateFDJacobianAtMidpointOnly", evalJacAtMidpointOnly, ParamNode::PASS );
      } else {
        evalJacAtMidpointOnly = false;
      }
    }
    if( methodNode->Has("QuasiNewton_FiniteDifferenceJacobian")){
      currentMethod = methodNode->Get("QuasiNewton_FiniteDifferenceJacobian");
      nonLinMethod_ = "HYST_FD_JACOBIAN";
      solutionApproach_ = SOLVE_JACOBIAN;
      FPApproach_ = NOFP;
      currentFPApproach_ = NOFP;
      currentMethod->GetValue( "IterationsTillUpdate", numIterationsTillUpdateOfJacobian_, ParamNode::PASS );
      if(currentMethod->Has("initialNumberFPSteps")){
        currentMethod->GetValue( "initialNumberFPSteps", initialNumberFPSteps_, ParamNode::PASS );
      } else {
        initialNumberFPSteps_ = 1;
      }
      
      if(currentMethod->Has("calculateFDJacobianAtMidpointOnly")){
        currentMethod->GetValue( "calculateFDJacobianAtMidpointOnly", evalJacAtMidpointOnly, ParamNode::PASS );
      } else {
        evalJacAtMidpointOnly = false;
      }

    }
    if( methodNode->Has("QuasiNewton_SecantMethod")){
      currentMethod = methodNode->Get("QuasiNewton_SecantMethod");
      currentMethod->GetValue( "towardsPreviousTimestep", towardsPreviousTimestep_, ParamNode::PASS );

      if(currentMethod->Has("initialNumberFPSteps")){
        currentMethod->GetValue( "initialNumberFPSteps", initialNumberFPSteps_, ParamNode::PASS );
      } else {
        initialNumberFPSteps_ = 1;
      }

      if(currentMethod->Has("includeDeltaStrain")){
        currentMethod->GetValue( "includeDeltaStrain", deltaMatStrain_, ParamNode::PASS );
      } else {
        deltaMatStrain_ = false;
      }

      if(currentMethod->Has("useDeltaInCouplingTensor")){
        currentMethod->GetValue( "useDeltaInCouplingTensor", deltaMatCoupling_, ParamNode::PASS );
      } else {
        deltaMatCoupling_ = false;
      }

      if(currentMethod->Has("calculateDeltaMatAtMidpointOnly")){
        currentMethod->GetValue( "calculateDeltaMatAtMidpointOnly", evalJacAtMidpointOnly, ParamNode::PASS );
      } else {
        evalJacAtMidpointOnly = false;
      }
      
      if(towardsPreviousTimestep_){
        nonLinMethod_ = "HYST_DeltaMat_TS_MOD";
        solutionApproach_ = SOLVE_DELTAMAT_TS_TOWARDS_IT;
        FPApproach_ = NOFP;
        currentFPApproach_ = NOFP;
      } else {
        nonLinMethod_ = "HYST_DeltaMat_IT";
        solutionApproach_ = SOLVE_DELTAMAT_IT;
        FPApproach_ = NOFP;
        currentFPApproach_ = NOFP;
      }

//      // currently we do not distinguish between -1 and +1 when checking
//      // for this flag in CoefFncHyst
//      // we check however, if flag != 0
//      // if 0 > no deltaMat will be computed!
//      // not needed at all in the momement
//      Integer deltaForm;
//      if(nonLinMethod_ == "HYST_DeltaMat_TS"){
//        deltaForm = -1;
//      } else if(nonLinMethod_ == "HYST_DeltaMat_IT"){
//        deltaForm = 1;
//      } else {
//        deltaForm = 0;
//      }
//      // deltaForm not used anymore
//      PDE_.SetFlagInCoefFncHyst("deltaForm",deltaForm);

    }

    if(evalJacAtMidpointOnly){
      PDE_.SetFlagInCoefFncHyst("evalJacAtMidpointOnly",1);
    } else {
      PDE_.SetFlagInCoefFncHyst("evalJacAtMidpointOnly",0);
    }
    
    std::cout << "++ Selected solution method for non-linear system: "  << solveFlagToString(solutionApproach_) << std::endl;

    if(initialNumberFPSteps_ != 0){
      std::cout << "++ " << initialNumberFPSteps_ << " initial fixpoint iterations (Global B, safety-factor 2.0) will be applied " << std::endl;
      FPGlobalC_ = 2.0;
      currentFPApproach_ = FP_GLOBAL_B;
      currentNonLinMethod_ = "HYST_FP_GLOBAL";
      currentSolutionApproach_ = SOLVE_GLOBAL_FIXPOINT_B;
    } else {
      currentNonLinMethod_ = nonLinMethod_;
      currentSolutionApproach_ = solutionApproach_;
      currentFPApproach_ = FPApproach_;
    }

    PDE_.SetDoubleFlagInCoefFncHyst("FPGlobalC",FPGlobalC_);
    PDE_.SetDoubleFlagInCoefFncHyst("FPLocalC",FPLocalC_);
    PDE_.SetFlagInCoefFncHyst("FPApproach",currentFPApproach_);

  }

  void SolveStepHyst::Init_Step1c_ReadStoppingCriteria(PtrParamNode& stoppingNode, bool isFailback){

    UInt numberOfCriteria = 0;

    if(stoppingNode == NULL){
      if(isFailback){
        // no failback criteria defined (ok)
        return;
      } else {
        EXCEPTION("At least one stopping criterion has to be defined.");
      }
    }

    PtrParamNode currentStoppingCriterion;

    if(stoppingNode->Has("residual_relative")){
      currentStoppingCriterion = stoppingNode->Get("residual_relative");

      // Read in parameter from current node
      StoppingCriterion currentCriterionParameter = Init_Step1c_ReadSingleStoppingCriterion(currentStoppingCriterion,"residual_relative");

      // Add counter
      currentCriterionParameter.checkingOrder_ = ++numberOfCriteria;

      // Insert
      if(isFailback){
        failbackCriteria_.insert(currentCriterionParameter);
      } else {
        stoppingCriteria_.insert(currentCriterionParameter);
      }
    }

    if(stoppingNode->Has("residual_absolute")){
      currentStoppingCriterion = stoppingNode->Get("residual_absolute");

      // Read in parameter from current node
      StoppingCriterion currentCriterionParameter = Init_Step1c_ReadSingleStoppingCriterion(currentStoppingCriterion,"residual_absolute");

      // Add counter
      currentCriterionParameter.checkingOrder_ = ++numberOfCriteria;

      // Insert
      if(isFailback){
        failbackCriteria_.insert(currentCriterionParameter);
      } else {
        stoppingCriteria_.insert(currentCriterionParameter);
      }
    }

    if(stoppingNode->Has("increment_relative")){
      currentStoppingCriterion = stoppingNode->Get("increment_relative");

      // Read in parameter from current node
      StoppingCriterion currentCriterionParameter = Init_Step1c_ReadSingleStoppingCriterion(currentStoppingCriterion,"increment_relative");

      // Add counter
      currentCriterionParameter.checkingOrder_ = ++numberOfCriteria;

      // Insert
      if(isFailback){
        failbackCriteria_.insert(currentCriterionParameter);
      } else {
        stoppingCriteria_.insert(currentCriterionParameter);
      }
    }

    if(stoppingNode->Has("increment_absolute")){
      currentStoppingCriterion = stoppingNode->Get("increment_absolute");

      // Read in parameter from current node
      StoppingCriterion currentCriterionParameter = Init_Step1c_ReadSingleStoppingCriterion(currentStoppingCriterion,"increment_absolute");

      // Add counter
      currentCriterionParameter.checkingOrder_ = ++numberOfCriteria;

      // Insert
      if(isFailback){
        failbackCriteria_.insert(currentCriterionParameter);
      } else {
        stoppingCriteria_.insert(currentCriterionParameter);
      }
    }

    if(!isFailback){
      if(stoppingCriteria_.size() == 1){
        std::cout << "++ Selected stopping criterion: "  << stoppingCriteria_.begin()->nameTag_ << " <= " << stoppingCriteria_.begin()->stoppingTolerance_ << std::endl;
      } else {
        std::set< StoppingCriterion >::iterator critIterator;
        std::cout << "++ List of selected stopping criteria: " << std::endl;

        UInt cnt = 1;
        for(critIterator = stoppingCriteria_.begin(); critIterator != stoppingCriteria_.end(); critIterator++){
          std::cout << "   (" << cnt++ << "): " << critIterator->nameTag_ << " <= " << critIterator->stoppingTolerance_ << std::endl;
        }
      }
    } else {
      if(failbackCriteria_.size() == 1){
        std::cout << "++ Selected failback criterion: "  << failbackCriteria_.begin()->nameTag_ << " <= " << failbackCriteria_.begin()->stoppingTolerance_ << " (first check in iteration " << failbackCriteria_.begin()->firstCheck_ << ")" << std::endl;
      } else if(failbackCriteria_.size() > 1){
        std::set< StoppingCriterion >::iterator critIterator;
        std::cout << "++ List of selected failback criteria: " << std::endl;

        UInt cnt = 1;
        for(critIterator = failbackCriteria_.begin(); critIterator != failbackCriteria_.end(); critIterator++){
          std::cout << "   (" << cnt++ << "): " << critIterator->nameTag_ << " <= " << critIterator->stoppingTolerance_ << " (first check in iteration " << critIterator->firstCheck_ << ")" << std::endl;
        }
      }
    }
  }

  StoppingCriterion SolveStepHyst::Init_Step1c_ReadSingleStoppingCriterion(PtrParamNode& singleStoppingNode,std::string nameTag){

    StoppingCriterion currentCriterion = StoppingCriterion();
    currentCriterion.nameTag_ = nameTag;

    // mandatory
    singleStoppingNode->GetValue("value",currentCriterion.stoppingTolerance_,ParamNode::PASS);

    // optional
    if(singleStoppingNode->Has("scaleToSystemSize")){
      singleStoppingNode->GetValue("scaleToSystemSize",currentCriterion.scaleToSystemSize_,ParamNode::PASS);
    } else {
      currentCriterion.scaleToSystemSize_ = false;
    }

    // currently not used
    currentCriterion.scalingFactor_ = 1.0;

    // firstCheck_: first iteration at which criterion shall be checked (for failback criteria which shall be checked
    //  after e.g., the 10th iteration for the first time
    if(singleStoppingNode->Has("firstCheck")){
      singleStoppingNode->GetValue("firstCheck",currentCriterion.firstCheck_,ParamNode::PASS);
    } else {
      currentCriterion.firstCheck_ = 1;
    }

    // checkingFrequency_: number of iterations that has to be passed before criterion is checked again;
    // e.g. 1 = check every iteration
    //      2 = check every second iteration
    if(singleStoppingNode->Has("checkingFrequency")){
      singleStoppingNode->GetValue("checkingFrequency",currentCriterion.checkingFrequency_,ParamNode::PASS);
    } else {
      currentCriterion.checkingFrequency_ = 1;
    }

    // name based parameter
    if(nameTag == "increment_absolute"){
      currentCriterion.isResidual_ = false;
      currentCriterion.isRelative_ = false;
    }
    if(nameTag == "increment_relative"){
      currentCriterion.isResidual_ = false;
      currentCriterion.isRelative_ = true;
    }
    if(nameTag == "residual_absolute"){
      currentCriterion.isResidual_ = true;
      currentCriterion.isRelative_ = false;
    }
    if(nameTag == "residual_relative"){
      currentCriterion.isResidual_ = true;
      currentCriterion.isRelative_ = true;
    }

//    /*
//     * to be set outside
//     */
//    // checkingOrder_: as name suggests, order in which multiple criteria shall be checked
//    UInt checkingOrder_;
//    Double scalingFactor_;

    return currentCriterion;
  }

  void SolveStepHyst::Init_Step1b_ReadLinesearches(PtrParamNode& linesearchNode){

    UInt numberOfLinesearches = 0;
    nonmonotonicArmijo_numEntries_ = 1;
    checkDirection_ = 0;

    if(linesearchNode == NULL){
      // no linesearch at all; overwrites all otherwise defined linesearches
//      std::cout << "linesearchNode = NULL" << std::endl;
      return;
    }

    /*
     * Iterate over all defined linesearches and check if it is defined in linesearchNode
     * If yes: add to set of all linsearches
     */
    if(linesearchNode->Has("None")){
      // no linesearch at all; overwrites all otherwise defined linesearches
//      std::cout << "linesearchNode has None" << std::endl;
      return;
    }

    PtrParamNode currentLinesearch;

    if(linesearchNode->Has("Backtracking_SmallestResidual")){
      currentLinesearch = linesearchNode->Get("Backtracking_SmallestResidual");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_TRIAL_AND_ERROR;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Backtracking_SmallestResidual");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("Backtracking_Armijo")){
      currentLinesearch = linesearchNode->Get("Backtracking_Armijo");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_BACKTRACKING_ARMIJO;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Backtracking_Armijo");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("Backtracking_NonMonotonicArmijo")){
      currentLinesearch = linesearchNode->Get("Backtracking_NonMonotonicArmijo");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_BACKTRACKING_ARMIJO_NONMONOTONIC;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Backtracking_NonMonotonicArmijo");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );

      // set history length
      nonmonotonicArmijo_numEntries_ = currentLinesearchParameter.historyLengthArmijo_;
    }

    if(linesearchNode->Has("Backtracking_GradientFreeArmijo")){
      currentLinesearch = linesearchNode->Get("Backtracking_GradientFreeArmijo");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_BACKTRACKING_GRADFREE;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Backtracking_GradientFreeArmijo");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("Exact_GoldenSection")){
      currentLinesearch = linesearchNode->Get("Exact_GoldenSection");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_GOLDENSECTION;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Exact_GoldenSection");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("Exact_QuadraticPolynomial")){
      currentLinesearch = linesearchNode->Get("Exact_QuadraticPolynomial");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_QUADINTERPOL;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"Exact_QuadraticPolynomial");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("SufficientDecrease_ThreePointPolynomial")){
      currentLinesearch = linesearchNode->Get("SufficientDecrease_ThreePointPolynomial");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_THREEPOINTPOLY;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"SufficientDecrease_ThreePointPolynomial");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    if(linesearchNode->Has("SufficientDecrease_VectorPolynomial")){
      currentLinesearch = linesearchNode->Get("SufficientDecrease_VectorPolynomial");

      // Assign right flag
      linesearchFlag currentLinesearchFlag = LS_VECBASEDPOLY;

      // Read in parameter from current node
      LinesearchParameter currentLinesearchParameter = Init_Step1b_ReadSingleLinesearch(currentLinesearch,"SufficientDecrease_VectorPolynomial");

      // Add counter
      currentLinesearchParameter.checkingOrder_ = ++numberOfLinesearches;

      // Insert
      specifiedLinesearches_.insert(currentLinesearchFlag);
      specifiedLinesearchesMap_.insert( std::pair<linesearchFlag, LinesearchParameter>(currentLinesearchFlag,currentLinesearchParameter) );
    }

    /*
     * finally read in global attributes
     */
    if(linesearchNode->Has("checkSolutionUpdateWithFDJacobian")){
      linesearchNode->GetValue("checkSolutionUpdateWithFDJacobian",checkDirection_,ParamNode::PASS);
    } else {
      checkDirection_ = 0;
    }

//    If multiple linesearches are applied, we have to specify which one to take. <br/>
//            Value = 1 (Default): take eta that leads to minimal residual  <br/>
//            Value = 2: take largest found eta (positive value) <br/>
//            Value = 3: take largest found eta (absolute value)  <br/>
    if(linesearchNode->Has("selectionCriterionForMultipleLS")){
      linesearchNode->GetValue("selectionCriterionForMultipleLS",LSSelectionCriterion_,ParamNode::PASS);
    } else {
      LSSelectionCriterion_ = 1;
    }

    if(specifiedLinesearchesMap_.size() == 1){
      std::cout << "++ Selected line-search: "  << specifiedLinesearchesMap_.begin()->second.nameTag_ << std::endl;
    } else {
      std::map< linesearchFlag, LinesearchParameter >::iterator LSIterator;
      std::cout << "++ List of selected line-searches: " << std::endl;

      UInt cnt = 1;
      for(LSIterator = specifiedLinesearchesMap_.begin(); LSIterator != specifiedLinesearchesMap_.end(); LSIterator++){
        std::cout << "   (" << cnt++ << "): " << LSIterator->second.nameTag_ << std::endl;
      }
    }
  }

  // TODO: change parameter in input; read in new parameter; remove unused input parameter
	void SolveStepHyst::Init_Step1_ReadNonLinDataHyst(){
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== Init_Step1_ReadNonLinDataHyst ==";
    LOG_TRACE(solvestephyst_main) << "===================================";

		PtrParamNode hystNode = solStrat_->GetHystNode();

    /*
     * Debug related flags
     */
    testInversion_ = 0;
    DEBUG_testJacobianApproximations_ = false;
    DEBUG_measureLSPerformance_ = false;
    DEBUG_compareLinesearches_ = 0;

    if( !hystNode ) {
      EXCEPTION("SolveStepHyst::ReadNonLinDataHyst --- No hyst node found in input file");
    } else {
      /*
       * Read in simple elements and attributes
       */
      hystNode->GetValue( "maxIter", nonLinMaxIter_, ParamNode::PASS );
      hystNode->GetValue( "loggingToFile", nonLinLogging_, ParamNode::PASS );
      hystNode->GetValue( "loggingToTerminal", minLoggingToTerminal_, ParamNode::PASS );

      iterationsTillReset_ = 0; // per default, we perform no reset
      PtrParamNode resetNode = NULL;
      if(hystNode->Has("resetSolutionVectorToZero")){
        resetNode = hystNode->Get("resetSolutionVectorToZero");
        resetNode->GetValue("iterationsTillReset",iterationsTillReset_,ParamNode::PASS);
      }

      PtrParamNode evalDepthNode = NULL;
      if(hystNode->Has("evaluationDepth")){
        evalDepthNode = hystNode->Get("evaluationDepth");
        if(evalDepthNode->Has("evaluateAtElementCenter")){
          evalDepth_ = 1;
        } else if(evalDepthNode->Has("evaluateAtIntegrationPoints_OneOperatorPerElement")){
          evalDepth_ = 2;
        } else if(evalDepthNode->Has("evaluateAtIntegrationPoints_OneOperatorPerIntegrationPoint")){
          evalDepth_ = 3;
        }
      } else {
        // Default: evalDepth_ = 2
        evalDepth_ = 2;
      }

      /*
       * Read in more complex subnodes with helper functions
       */
      PtrParamNode methodNode = hystNode->Get("solutionMethod");
      Init_Step1a_ReadSolutionMethod(methodNode);

      PtrParamNode linesearchNode = NULL;
      if (hystNode->Has("lineSearch")){
        linesearchNode = hystNode->Get("lineSearch");
      }
      Init_Step1b_ReadLinesearches(linesearchNode);

      PtrParamNode stoppingNode = NULL;
      if (hystNode->Has("stoppingCriteria")){
        stoppingNode = hystNode->Get("stoppingCriteria");
      }
      Init_Step1c_ReadStoppingCriteria(stoppingNode,false);

      stoppingNode = NULL;
      if (hystNode->Has("failbackCriteria")){
        stoppingNode = hystNode->Get("failbackCriteria");
      }
      Init_Step1c_ReadStoppingCriteria(stoppingNode,true);

      // read in optional tests
      PtrParamNode testParamNode = NULL;
      if (hystNode->Has("performTests")){
        testParamNode = hystNode->Get("performTests");

        if(testParamNode->Has("testJacobianApproximations")){
          DEBUG_testJacobianApproximations_ = true;
        } else {
          DEBUG_testJacobianApproximations_ = false;
        }

        testNode_ = NULL;
        if( testParamNode->Has("testHysteresisOperator") ) {
          testInversion_ = 1;
          testNode_ = testParamNode->Get("testHysteresisOperator");
        } else {
          testInversion_ = 0;
        }
      }

      /*
       * Important: evaluationDepth has to be set first as it triggers the
       * initializstion of the storage
       */
      PDE_.SetFlagInCoefFncHyst("evaluationDepth",evalDepth_);
      //PDE_.SetFlagInCoefFncHyst("deltaForm",deltaForm);
      //PDE_.SetFlagInCoefFncHyst("measurePerformance",measurePerformance_);
    }
    initStep1Done_ = true;
	}

  void SolveStepHyst::Init_Step2_InitGlobalArrays(){
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== Init_Step2_InitGlobalArrays   ==";
    LOG_TRACE(solvestephyst_main) << "===================================";

    // Create a dummy vector first which contains the correct size
    // this one can be used throughout the code to setup local SBM_Vectors, too
    // Important note:
    //  in order to obtain a full copy of this dummy, we must not set
    //  SBM_Vector dummy2 = dummy;  < creates just a pointer!
    //  SBM_Vector dummy2; dummy2 = dummy; < creates copy of dummy
    zeroVec_.Resize(feFunctions_.size());
    systemSize_ = feFunctions_.size();

    matrices = PDE_.GetMatrixDerivativeMap();
    UInt rhsSize = 0;
    for(fncIt = feFunctions_.begin(); fncIt != feFunctions_.end();++fncIt ){
      rhsSize = fncIt->second->GetSingleVector()->GetSize();
      FeFctIdType id = fncIt->second->GetFctId();
			zeroVec_.SetSubVector(new Vector<Double>(),id);
      zeroVec_.GetPointer(id)->Resize(rhsSize);
    }
    zeroVec_.Init();

    // Vectors defined in base class
    //    solVec_ > direct access to current solution vector in algsys
    //    rhsVec_ > direct access to current rhs vector in algsys

    // Global SBM_Vectors and SBM_Vector arrays
    // a) solution vectors
    solVec_n_0_ = zeroVec_; // solution at start of time step
    solVec_n_kPrev_ = zeroVec_; // solution from previous iteration
    solVec_n_k_ = zeroVec_; // solution from current iteration
    solIncrement_n_k_ = zeroVec_; // current solution increment = solVec_n_k - solVec_n_kPrev
    bestSol_ = zeroVec_;

    // b) residuals and rhs vectors
    resVec_n_0_ = zeroVec_;
    resVec_n_kPrev_ = zeroVec_;
    resVec_n_k_ = zeroVec_;
    resIncrement_n_k_ = zeroVec_;

    // c) rhs and rhs-parts
    linRhsVec_n_0_ = zeroVec_;
    predictorCorrectorUpdate_n_0_ = zeroVec_;
    nonLinRhsVec_n_k_ = zeroVec_;

    historyLength_ = 0;
    // to be put into input file
    // FOR BGM as described in Kelly 1995, there is no overflow intended
    //  i.e. we store all increments; in L-BFGS we only store the latest n vectors
    hisotryLengthMax_ = nonLinMaxIter_+1;//10;

    resIncrements_history_ = new SBM_Vector[hisotryLengthMax_];
    solIncrements_history_ = new SBM_Vector[hisotryLengthMax_];
    etaLinsearch_history_ = Vector<Double>(hisotryLengthMax_);

    for(UInt i = 0; i < hisotryLengthMax_; i++){
      solIncrements_history_[i] = zeroVec_;
      resIncrements_history_[i] = zeroVec_;
    }

    // Assembly related
    lastUpdateVecForNonLinRHS_ = zeroVec_;
    lastUpdateVecForHystFreeSystem_ = zeroVec_;
    lastUpdateVecForFDJacobianSystem_ = zeroVec_;
    lastUpdateVecForDeltaMatJacobianSystem_ = zeroVec_;
  }

  /*
   * III. Solving
   * > called during each timestep
   * > outer loop
   * > n = time step counter
   */
  void SolveStepHyst::SolveStepTrans() {
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== SolveStepTrans (" << timestepCnt_+1 << ") ============";
    LOG_TRACE(solvestephyst_main) << "===================================";

		timestepCnt_++;
    iterationCntCurrentTS_ = 0;
    iterationOfLatestReset_ = 0;
    UInt resetCnt = 0;
    bool nonLinConverged = false;

    LOG_DBG(solvestephyst_main) << "--- Norm of solution vector at START of time step: " << solVec_.NormL2() << "---";

    Solve_Step1_PrepareSolveStep();

    // currently no reset implemented
    allowedNumResets_ = 0;
    while(!nonLinConverged){
      LOG_TRACE(solvestephyst_main) << "== Start non-linear solution attempt " << resetCnt+1 << " ==";
      nonLinConverged = Solve_Step2_NonLinIteration(resetCnt);

      resetCnt++;
      if(!nonLinConverged){
        if(resetCnt > allowedNumResets_){
          EXCEPTION("Maximal number of solution attempts ("<<allowedNumResets_+1<< ") for time step "<<timestepCnt_<<" exceeded.");
          break;
        }
      }
    }

    Solve_Step3_FinishSolveStep();

    LOG_DBG(solvestephyst_main) << "--- Norm of solution vector at END of time step: " << solVec_.NormL2() << " ---";;
  }

  void SolveStepHyst::Solve_Step1_PrepareSolveStep(){
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== Solve_Step1_PrepareSolveStep ===";
    LOG_TRACE(solvestephyst_main) << "===================================";
    
    /*
     * Test hysteresis operator if specified in input xml
     */
    if(testInversion_ == 1){
      std::cout << "++ Test hysteresis operator" << std::endl;
			PDE_.TestInversionOfHystOperator(testNode_);
			testInversion_ = 0;
		}
    if(DEBUG_testJacobianApproximations_){
      std::cout << "++ Test different approximations for the Jacobian matrix" << std::endl;
      PDE_.SetFlagInCoefFncHyst("TestJacobianApproximations",0);
      DEBUG_testJacobianApproximations_ = false;
    }
    
    /*
     * we only support incremental stepping, i.e. we just compute
     * du = u_new - u_old in each time step
     * > force incremental
     */
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt ){
      if( fncIt->second->GetTimeScheme()->isIncremental() == false){
        //        WARN("SolveStepHyst: Incremental stepping will be enforced");
        fncIt->second->GetTimeScheme()->forceIncremental();
      }
    }

    /*
     * Init time stepping
     */
    //obtain the number of stages of time stepping scheme
    UInt numStages = feFunctions_.begin()->second->GetTimeScheme()->GetNumStages();
    if ( numStages > 1 ){
      EXCEPTION("SolveStepHyst: currently just single-stage time-stepping allowed");
    }
    UInt stage = 0;

    bool updatePredictor = ( PDE_.IsIterCoupled() == false || couplingIter_ == 0 );
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->BeginStep(updatePredictor);
    }

    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      /*
       * we used to initialize the stageSol vector here, but this is no longer
       * necessary here, as we no longer use this vector;
       * instead we update by passing solVec_
       */
      fncIt->second->GetTimeScheme()->InitStage(stage,actTime_,PDE_.GetDomain());
    }

    /*
     * Initialize/Reset solutions, residual, rhs etc
     */
    solVec_n_0_ = solVec_; // solution at start of time step > take solution from last timestep
    solVec_n_kPrev_ = solVec_; // solution from previous iteration > previous iteration = last timestep
    solVec_n_k_ = zeroVec_; // solution from current iteration > not determined yet; set to 0
    solIncrement_n_k_ = zeroVec_; // current solution increment > not determined yet; set to 0
    bestSol_ = zeroVec_; // > not determined yet; set to 0

    // b) residuals and rhs vectors
    resVec_n_0_ = zeroVec_; // residual at start of time step > has to be evaluated based on new rhs
    resVec_n_kPrev_ = zeroVec_;
    resVec_n_k_ = zeroVec_;
    resIncrement_n_k_ = zeroVec_;

    // c) rhs and rhs-parts
    linRhsVec_n_0_ = zeroVec_;
    predictorCorrectorUpdate_n_0_ = zeroVec_;
    // if we reset the nonlin rhs vector, we must force reassembly to make sure to set it up again

    forceRHSevaluation_ = true;
    nonLinRhsVec_n_k_ = zeroVec_;

		solIncrement_n_k_.Init();
    solVec_n_kPrev_ = solVec_;
    solVec_n_0_ = solVec_;
    bestSol_.Init();
    //
    //    // residual and rhs
    //		resIncrement_n_k_.Init();
    //		rhsVec_.Init();
    //    algsys_->InitRHS(rhsVec_);

    /*
     * Set/Reset flags
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

    // reset all build counters such that each possible system gets assembled at least once if needed
    latestBuildCounter_[HYSTFREE] = -1;
    latestBuildCounter_[FDJACOBIAN] = -1;
    latestBuildCounter_[DELTAMATJACOBIAN_LAST_TS] = -1;
    latestBuildCounter_[DELTAMATJACOBIAN_LAST_IT] = -1;
    latestBuildCounter_[REUSED] = -1;
    latestBuildCounter_[HYSTFREE_GLOBALLY_SCALED_B_VERSION] = -1;
    latestBuildCounter_[HYSTFREE_LOCALLY_SCALED_B_VERSION] = -1;
    latestBuildCounter_[HYSTFREE_SCALED_H_VERSION] = -1;
    latestBuildCounter_[HYSTFREE_LOCALLY_SCALED_H_VERSION] = -1;

    LinRHSRequiresAssembly_ = true;
    PredictorCorrectorRequiresAssembly_ = true;

    /*
     * 24.6.2020 - reset evaluation purpose to 1 (assemble);
     * at the end of a timestep this flag is set to 4 (output); thus the hysteresis operator
     * would only be evaluated at the center of each element but in fact the flag is already set to 1
     * in the precomputePolarization method; however it is reset afterwards to the previous state, i.e., to 4
     * in this case; for estimating the local fp slope flag 4 will be active again which is not wanted
     * > set flag here to 1 so that it will be 1 during the fp slope estimatation further below
     */
    // note: this change actually breaks the testcase to piezo-voltage-excitation; it seems like
    // the deltaMat formulation works more reliable if the deltaMat is evaluated only at the midpoint
    // (which was the case before setting the evaluationPurpose to 1)
    // > use same flag as for localized fp method to force midpoint evaluation
    // > seems to work but is slower ...
    int evaluationPurpose = 1;
    PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose);
    
    // before we evaluate the hyst operators, we have to check if we solve the (magnetic) system
    // by fixpoint scheme in its H-Version
    // in this case the forward hysteresis operator can be applied directly without performing an inversion step
    if( (solutionApproach_ == SOLVE_GLOBAL_FIXPOINT_H) || (solutionApproach_ == SOLVE_LOCAL_FIXPOINT_H) ){
      PDE_.SetFlagInCoefFncHyst("SetFPH",1);
    } else {
      PDE_.SetFlagInCoefFncHyst("SetFPH",0);
    }

    /*
     * Prepare Hyst operators
     */
    /*
     * For magnetic PDE, we have the problem, that the matrix that is used
     * to compute H (i.e. mu/nu) might depend on H/M itself
     * > compute mu/nu for inversion only once each iteration
     */
    // allow setting of matrices
    PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
    // EvaluateHystOperators; flag = 3 > set Matrix for Inverion
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",3);

    /*
     * disable output of switching and rotation state as BMP images (even if flag
     * printOut is set to 1 in material file)
     * > too expensive; could even be deprecated by now (needs checking)
     */
    PDE_.SetFlagInCoefFncHyst("allowBMP",true);
    PDE_.SetFlagInCoefFncHyst("resetEstimator",true);
    PDE_.SetFlagInCoefFncHyst("resetGlobalFPTensors",true);

    /*
     * Assemble linear system, linear rhs and initial residual
     * Note: in case of fixPoint iteration, we might have to use scaled material tensors
     * > retrieve correct assembly flags from solveFlag
     */
    currentAssemblyFlag_ = NOTSET;
    currentBuildFlag_ = NOTSET;
    linearSystemBackuped_ = false;
    forceReassembly_ = true;
    currentRESandRHSFlag_ = NOTSET_VEC;

    matrixFlag flag_for_systemMatrix;
    matrixFlag flag_for_resEvalMatrix;
    vectorFlag flag_for_rhsVector;
    vectorFlag flag_for_resVector;
    bool estimatorUseful;

    // TODO: add check for failback case
    // > use solutionApproachFAILBACK_ in that case
    // general question in that regard: what steps have to be redone after a reset/failback?
    Helper_GetParamsFromSolveFlag(currentSolutionApproach_, flag_for_systemMatrix, flag_for_resEvalMatrix, flag_for_rhsVector, flag_for_resVector, estimatorUseful);


    // for fixpoint methods, we have to scale the material tensor to ensure convergence
    // (see e.g. "Locally Convergent Fixed-Point Method for Solving Time-Stepping Nonlinear Field Problems" - E. Dlala et al. 2007)
    // this has to be done only once each time step
    // for Fixpoint-B we have two options:
    //  a) global approach, i.e. nu_FP = (nu_min + nu_max)/2 \approx nu_0/2
    //  b) local approach, i.e. nu_FP = nu(B_current) \approx FD_Jacobian
    int timeLevelMaterial;
    int FPMaterialTensor;
    int FPRHSCorrectionTensor;
    int timeLevelDeltaMatPol;
    int timeLevelDeltaMatStrain;
    int timeLevelDeltaMatCoupling;
    int inversionApproach;
    Helper_GetHystOperatorParamsFromMatrixFlag(flag_for_systemMatrix, timeLevelMaterial, FPMaterialTensor, FPRHSCorrectionTensor,
            timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling, inversionApproach);
//
//    std::cout << "Matrix flag for system matrix. " << matrixFlagToString(flag_for_systemMatrix) << std::endl;
//    std::cout << "FPMaterialTensor. " << FPMaterialTensor << std::endl;
//    std::cout << "SetFPRHSCorrectionTensors. " << FPRHSCorrectionTensor << std::endl;

    // eval fixpoint reluctivity/permittivity > only once each time step! do this during Solve_Step1_PrepareSolveStep
    // BUT if rhs uses standard matrix, we have to switch off the matrix flag during matrix assembly for rhs and switch on during
    // system assembly!
    ReEvalFPMaterialTensors_ = true;
    if((FPMaterialTensor != 0)&&(ReEvalFPMaterialTensors_)){
      PDE_.SetFlagInCoefFncHyst("EvalFPMaterialTensors",FPMaterialTensor); // only once per timestep OR when switching
      // from initial globalFP to localFP!
      ReEvalFPMaterialTensors_ = false;
    }
    PDE_.SetFlagInCoefFncHyst("SetFPMatrixFlag",FPMaterialTensor); // during each assembly!

    // rhs and system might use different matrices (e.g. for Jacobian or even for Fixpoint), rhsCorrection term must fit
    // to rhs update matrix!
    Helper_GetHystOperatorParamsFromMatrixFlag(flag_for_resEvalMatrix, timeLevelMaterial, FPMaterialTensor, FPRHSCorrectionTensor,
            timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling, inversionApproach);

//    std::cout << "Matrix flag for rhs update matrix. " << matrixFlagToString(flag_for_resEvalMatrix) << std::endl;
//    std::cout << "FPMaterialTensor. " << FPMaterialTensor << std::endl;
//    std::cout << "SetFPRHSCorrectionTensors. " << FPRHSCorrectionTensor << std::endl;
//
//    std::cout << "Solve_Step1_PrepareSolveStep" << std::endl;
    PDE_.SetFlagInCoefFncHyst("SetFPCorrectionFlag",FPRHSCorrectionTensor); // during each rhs assembly


    // todo remove flag; it is never set to true and not really beneficial, too
    // AssembleSystem checks if assembly is required based on non-linearities and so on, i.e. it has more
    // advanced checks
    if(linearSystemBackuped_ == false){
      LOG_DBG(solvestephyst_main) << "--- Setup and store LIN system ---";
      // flag_for_resEvalMatrix usually is HYST_FREE except for Fixpoint cases where
      // it is HYSTFREE_GLOBALLY_SCALED_B_VERSION, HYSTFREE_LOCALLY_SCALED_B_VERSION or HYSTFREE_SCALED_H_VERSION, respective
//      std::cout << "Solve_Step1_PrepareSolveStep - Assemble and build HYSTFREE" << std::endl;
      AssembleSystem(HYSTFREE, solVec_,false);
      //      AssembleSystem(solVec_, 0);
      // backup system > flag = false
      BuildSystem(false);
    }

    if( (flag_for_resEvalMatrix == HYSTFREE_GLOBALLY_SCALED_B_VERSION) ||
        (flag_for_resEvalMatrix == HYSTFREE_LOCALLY_SCALED_B_VERSION) ||
        (flag_for_resEvalMatrix == HYSTFREE_LOCALLY_SCALED_H_VERSION) ||
        (flag_for_resEvalMatrix == HYSTFREE_SCALED_H_VERSION) ){
      LOG_DBG(solvestephyst_main) << "--- Setup and store LIN system ---";
      // flag_for_resEvalMatrix usually is HYST_FREE except for Fixpoint cases where
      // it is HYSTFREE_GLOBALLY_SCALED_B_VERSION, HYSTFREE_LOCALLY_SCALED_B_VERSION or HYSTFREE_SCALED_H_VERSION, respective
//      std::cout << "++++ Assemble and build flag_for_resEvalMatrix = " << matrixFlagToString(flag_for_resEvalMatrix) << std::endl;
      AssembleSystem(flag_for_resEvalMatrix, solVec_,false);
      //      AssembleSystem(solVec_, 0);
      BuildSystem(false);

      // matrix for rhs and system solution might be different but both should be constant if no further nonlinearities
      // are present > evaluate during init phase and reload later
//      std::cout << "++++ Assemble and build flag_for_systemMatrix = " << matrixFlagToString(flag_for_systemMatrix) << std::endl;
      AssembleSystem(flag_for_systemMatrix, solVec_,false);
      //      AssembleSystem(solVec_, 0);
      BuildSystem(false);
    }

    // setup linear rhs
    if(LinRHSRequiresAssembly_){
      LOG_DBG(solvestephyst_main) << "--- Setup LIN rhs ---";

			linRhsVec_n_0_.Init();
      algsys_->InitRHS(linRhsVec_n_0_);
      assemble_->AssembleLinRHS();
      algsys_->GetRHSVal( linRhsVec_n_0_ );

			LinRHSRequiresAssembly_ = false;

      LOG_DBG(solvestephyst_main) << "--- Norm of LIN rhs: " << linRhsVec_n_0_.NormL2() << " ---";
		}

    // setup rhs contribution due to time step
    if(PredictorCorrectorRequiresAssembly_){
      LOG_DBG(solvestephyst_main) << "--- Setup Predictor-Corrector for time stepping ---";

			predictorCorrectorUpdate_n_0_.Init();
      algsys_->InitRHS(predictorCorrectorUpdate_n_0_);

			if(trans_){
				// remember: (currently) we only support single stage timestepping!
				UInt stage = 0;

				std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
				std::map<FEMatrixType,Integer> matrices = PDE_.GetMatrixDerivativeMap();
				std::map<FEMatrixType,Integer>::iterator matIt;

				// update RHS according to time stepping
				// here we only want to get the update due to predictor corrector, but
				// not the update for incremental stepping
        //  > this is done later in step 4
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
			algsys_->GetRHSVal( predictorCorrectorUpdate_n_0_ );

      PredictorCorrectorRequiresAssembly_ = false;

      LOG_DBG(solvestephyst_main) << "--- Norm of Predictor-Corrector update: " << predictorCorrectorUpdate_n_0_.NormL2() << " ---";
		}

    // setup initial residual (i.e. residual at current timestep
    // with solution from last ts
    // > pass iteration counter = 0 to store directly to resVec_n_0_
    LOG_DBG(solvestephyst_main) << "--- Setup Initial residual ---";

//    std::cout << "+++++ Setup initial residual" << std::endl;
    SetupRESorRHS(flag_for_resEvalMatrix, flag_for_resVector, solVec_, solVec_, solVec_, resVec_n_0_, false);
//    SetupRESVector(resVec_n_0_,solVec_, solVec_, solVec_, 0);
    resVec_n_0_Norm_ = resVec_n_0_.NormL2();

//    std::cout << "Norm of current solVec_: " << solVec_.NormL2() << std::endl;
//    std::cout << "Norm of initial residual: " << resVec_n_0_.NormL2() << std::endl;
//
    LOG_DBG(solvestephyst_main) << "--- Norm of initial residual: " << resVec_n_0_.NormL2() << " ---";

    // currently not used
    performEstimatorStep_ = false;
    if(performEstimatorStep_ && estimatorUseful){
      LOG_DBG(solvestephyst_main) << "--- Estimator step - solve linear system for current res ---";

      // backup current solution to storage for last iteration
      // NOTE: previously, hyst values were restored from last TS without
      // setting previous hyst values first; i.e. the later restore command would load
      // zeros in that case
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
      /*
       * Solve the system
       *  K_eff*du = residual
       * to estimate in which direction du will go in order to improve
       * the approximation of later updates using Jacobian/DeltaMatrix approach
       *
       * for this step, we need to set IDBC and deltaIDBC as we calculate du
       * as the update towards the solution of the previous timestep (solVec_)
       */
      bool setIDBC = true;
      bool deltaIDBC = true;

      algsys_->GetRHSVal( rhsVec_ );

      // IMPORTANT: set the correct rhs before calling solve!
      algsys_->InitRHS( resVec_n_0_ );

//      std::cout << "Perform estimatior -- currentBuildFlag = " << matrixFlagToString(currentBuildFlag_) << std::endl;
      if(currentBuildFlag_ != HYSTFREE){
        if(linearSystemBackuped_ == false){
//          std::cout << "Assemble and build hystFree" << std::endl;
          LOG_DBG2(solvestephyst_main) << "--- SETUP and STORE LIN system for Estimator step ---";
          AssembleSystem(HYSTFREE,solVec_,false);//, 0);
          //          AssembleSystem(solVec_, 0);
          BuildSystem(false);
        } else {
          LOG_DBG2(solvestephyst_main) << "--- LOAD LIN system for Estimator step ---";
//          std::cout << "Load backup" << std::endl;
          Helper_LoadSystemBackup(false,false,false,true,SYSTEM_HYSTFREE);
        }
      }

      algsys_->Solve(setIDBC,deltaIDBC);
      algsys_->GetSolutionVal(solIncrement_n_k_, setIDBC, deltaIDBC );

      // full step in order to satisfy BC
      solVec_.Add(1.0,solIncrement_n_k_);

      LOG_DBG(solvestephyst_main) << "--- Norm of estimated solution (full step): " << solVec_.NormL2() << " ---";

      // evaluate hyst operators without writing the memory and without
      // reevaluationg the matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      // store the previously computed values of the hyst operator to the
      // estimator storage in coefFncHyst
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",2);

      // reset solution and hyst current entries of hyst operator
      // to values from lastTS
      solVec_ = solVec_n_0_;

      // relaod backup state
      PDE_.SetFlagInCoefFncHyst("RestorePreviousHystValues",0);
    }

    /*
     * backup solution and linsearch approaches in case that we have to
     * reload them after a reset
     */
    nonLinBackup_ = nonLinMethod_;
    lineSearchBackup_ = lineSearch_;
//    towardsLastIterationBackup_ = towardsLastIteration_;
    deltaMatCouplingBackup_ = deltaMatCoupling_;
    deltaMatStrainBackup_ = deltaMatStrain_;
  }

  bool SolveStepHyst::Solve_Step2_NonLinIteration(UInt resetCnt) {
    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== Solve_Step2_NonLinIteration  ===";
    LOG_TRACE(solvestephyst_main) << "===================================";

    /*
     * Local params
     */
//    Double bestResNorm = 1e15;
    bool converged = false;
    std::ostringstream convergenceLogger;

    /*
     * Setup flags
     */
    //    UInt iterationCounter = 0; // iteration counter for this timestep
    UInt iterationCounterLog = 0; // iteration couter for this timestep; gets reset when system gets reset
    bool performOneMoreStep;
//    bool reset = false;
//    bool failback = false;

//    Double residualErrRel = 0.0;
//    Double residualErrAbs = 0.0;

//    Double incrementalErrRel = 0.0;
//    Double incrementalErrAbs = 0.0;
//    Double solL2Norm = 0.0;

    // setup linear rhs and get norm
    //    Double loadFactor = 1.0;
    //    Double normOfLinRhs = SetLinRHS(loadFactor);
    /*
     * For no-load setups, we might use the norm of the system instead of normOfLinRhs
     */
//    Double normOfSystem; // norm of K_eff(u_k+1)u_k+1

    // keep track of residual; if it is not decreasing over several iterations > reset (see below)
    UInt numResets = 0;
//    std::deque<Double> resNormHistory(minNumberItTillReset_, 50000);

    // initialize with 0; here check for maximum
    resNormHistoryArmijo_ = std::deque<Double>(nonmonotonicArmijo_numEntries_, 0);

    bool steppingTowardsLastTS = false;
    if( (currentNonLinMethod_ == "HYST_FixPoint_TS")||(currentNonLinMethod_ == "HYST_DeltaMat_TS") ){
      steppingTowardsLastTS = true;
    } else if(currentNonLinMethod_ == "HYST_DeltaMat_TS_MOD"){
      // although called stepping towards ts, we step towards last iteration
      // only the deltaMat is computed towards ts
      steppingTowardsLastTS = false;
    }

    // after reset, we step towards zero;
    // if steppingTowardsLastTS = true, we start from 0
    bool steppingTowardsZero = false;

    // if steppingTowardsZero is true, we have to set full idbc
    // otherwise, we just have to set the difference towards the previous timestep
    // as the effect of the previous idbc values are already contained in the solution
    //    bool deltaIDBC;
    bool setIDBC;
//    bool lastTS = false;
    Double etaLinesearch;

    /*
     * Testsection
     * 11.10.2018:
     *  use mixture of deltaMat and fixpoint iteration in the following form:
     *  - compute deltaMat, solve system, check residual and increment
     *    > if increment << incCriterion and residual > resCriterion
     *        compute new deltaMat (as for standard deltaMat approach)
     *    > if increment >= incCriterion and residual > resCriterion
     *        keep deltaMat and do fixpoint iteration towards the current state
     *    > if increment < incCriterion and res < resCriterion
     *        stop
     *
     */
    disableDeltaMat_= false;
    reuseDeltaMat_ = false;

    // reset residual from previous time steps
    resVec_n_k_.Init();

    std::stringstream collectedLoggingMSG;

    residualErrAbsPrev_ = resVec_n_0_Norm_;
//    UInt forceReset = 0;

    bool storeIncrements = false;
    bool setStartingResForBGM = false;
    SBM_Vector startingResForBGM;
    SBM_Vector startingResForBGM_editable;
    if( (solutionApproach_ == SOLVE_BGM)||(solutionApproach_ == SOLVE_BGM_ENHANCED) ){
      storeIncrements = true;
      startingResForBGM = zeroVec_;
      startingResForBGM_editable = zeroVec_;
    }
    bool forceJacobian = false;
    bool forceReuseOfJacobian = false;

    /*
     *  Begin iterations
     */
    do {
//      std::cout << "Begin of iteration " << totalIterationCnt_ + 1 << std::endl;
          // test: check if convergence can be achieved if piezo coupling tensor is updated only once per timestep
//    PDE_.SetFlagInCoefFncHyst("resetReevalCouplingTensor",true);
      iterationCntCurrentTS_++;
      totalIterationCnt_++;
      iterationCounterLog++;

      if((currentNonLinMethod_ == "HYST_DeltaMat_TS_MOD")&&((iterationCntCurrentTS_ - iterationOfLatestReset_) == 1)){
        // although called stepping towards ts, we step towards last iteration
        // only the deltaMat is computed towards ts
        // during initial iteration force Jacobian for this particular method
        forceJacobian = true;
      }

      LOG_TRACE2(solvestephyst_main) << "== Start iteration " << iterationCntCurrentTS_ << " (" << totalIterationCnt_ << ") of time step " << timestepCnt_ << " ==";

      if(steppingTowardsLastTS || steppingTowardsZero || ( iterationCntCurrentTS_ == (iterationOfLatestReset_+1) && couplingIter_ == 0 ) ){
        setIDBC = true;
      } else {
        setIDBC = false;
      }

      /*
       * for BGM, we perform one standard step first to embody IDBC and to
       * get an initial guess; afterwards we compute updates without solving the
       * system
       */
      if( (currentSolutionApproach_ == SOLVE_BGM)||(currentSolutionApproach_ == SOLVE_BGM_ENHANCED) ){
        if(setIDBC){
          // to incorporate boundary conditions and/or to get an initial estimate
          // perform a step using FDJacobian
          forceJacobian = true;
          setStartingResForBGM = true;
        } else {
          // reuse Jacobian
          forceReuseOfJacobian = true;
        }
      } else if(currentSolutionApproach_ == SOLVE_CHORD){
//        std::cout << "Solve via Chord method" << std::endl;
        if((setIDBC) || (iterationCounterLog > initialNumberFPSteps_)){
//          std::cout << "setIDBC > forceJacobian" << std::endl;
          // to incorporate boundary conditions and/or to get an initial estimate
          // perform a step using FDJacobian
          forceJacobian = true;
        } else {
//          std::cout << "NO setIDBC > forceReuseOfJacobian" << std::endl;
          // reuse Jacobian
          forceReuseOfJacobian = true;
        }
      }

      // determine actual solution approach
      if(forceJacobian){
        currentSolutionApproach_ = SOLVE_JACOBIAN;
      } else {
        if(resetCnt != 0){
          currentSolutionApproach_ = solutionApproachFailback_;
        } else {
          currentSolutionApproach_ = solutionApproach_;
        }
      }

//      std::cout << "iterationCounterLog = " << iterationCounterLog << std::endl;
      if((iterationCounterLog <= initialNumberFPSteps_)){
        if((solutionApproach_ == SOLVE_GLOBAL_FIXPOINT_H) || (solutionApproach_ == SOLVE_LOCAL_FIXPOINT_H)){
//          std::cout << "Perform initial fixpoint steps (global H) number << " << iterationCntCurrentTS_ << std::endl;
          currentNonLinMethod_ = "HYST_FP_GLOBAL_H";
          currentSolutionApproach_ = SOLVE_GLOBAL_FIXPOINT_H;
          currentFPApproach_ = FP_GLOBAL_H;
          PDE_.SetFlagInCoefFncHyst("FPApproach",currentFPApproach_);
        } else {
//          std::cout << "Perform initial fixpoint steps (global B) number << " << iterationCntCurrentTS_ << std::endl;
          currentNonLinMethod_ = "HYST_FP_GLOBAL";
          currentSolutionApproach_ = SOLVE_GLOBAL_FIXPOINT_B;
          currentFPApproach_ = FP_GLOBAL_B;
          PDE_.SetFlagInCoefFncHyst("FPApproach",currentFPApproach_);
        }
        forceJacobian = false;
      } else if(iterationCounterLog == initialNumberFPSteps_ + 1) {
        // reset method and fp flags
        currentNonLinMethod_ = nonLinMethod_;
        currentSolutionApproach_ = solutionApproach_;
        currentFPApproach_ = FPApproach_;
        PDE_.SetFlagInCoefFncHyst("FPApproach",currentFPApproach_);
        ReEvalFPMaterialTensors_ = true;
      }
//      else {
//        std::cout << "iterationCntCurrentTS_: " << iterationCntCurrentTS_ << std::endl;
//      }

      if(currentSolutionApproach_ == SOLVE_JACOBIAN){
//        std::cout << "SOLVE JACOBIAN" << std::endl;
//        std::cout << "iterationCntCurrentTS_: " << iterationCntCurrentTS_ << std::endl;
//        std::cout << "iterationOfLatestReset_: " << iterationOfLatestReset_ << std::endl;
//        std::cout << "numIterationsTillUpdateOfJacobian_: " << numIterationsTillUpdateOfJacobian_ << std::endl;
        if((iterationCntCurrentTS_ - iterationOfLatestReset_ - 1 - initialNumberFPSteps_)%numIterationsTillUpdateOfJacobian_ == 0){
//          std::cout << "Iteration " << iterationCntCurrentTS_ << " - Compute new Jacobian approximation" << std::endl;
          forceReuseOfJacobian = false;
        } else {
//          std::cout << "Iteration " << iterationCntCurrentTS_ << " - Reuse previous Jacobian approximation" << std::endl;
          forceReuseOfJacobian = true;
        }
      }

      if(forceJacobian){
        // after reset for example, we really want to recompute the Jacobian even if iteration counter does not fit;
        // forceJacobian > forceReuseOfJacobian
        // TODO: count from this point anew, so that current Jacobian lasts for multiple iterations even if the next
        // one would not trigger forceReuseOfJacobian above
        forceReuseOfJacobian = false;
      }

      /*
       * Step 0.2: Reset solution UPDATE vector (if necessary)
       * (this vector is the one that gets subtracted from the rhs, i.e. the
       * one the stepping is going to)
       * > NOTE: do not overwrite solVec_ here, as the system still needs to be
       *         assembled using the current solution!
       *
       */
      if(steppingTowardsZero){
        solVecForUpdate_.Init();
      } else if(steppingTowardsLastTS){
        solVecForUpdate_ = solVec_n_0_;
      } else {
        // DOES NOT WORK IF COMMENTED OUT, reasons not clear yet
        //if((iterationCounterLog-1)%numIterationsTillUpdateOfLastITValues == 0){
        solVecForUpdate_ = solVec_; // solVec_n_kPrev_; < solVecLastIT is actually the solution from two iterations ago
        // i.e. during iteration k+1, solVec = sol_k and solVecLastIT = sol_k-1
        //}
      }

      /*
       * NOTE regarding BGM
       * the implemented form of the algorithm assumes that the Jacobian of the residual (i.e. the system to be solved)
       * is in the same order of magnitude as the identity; in that case, the solution update and the residual of the
       * system are in the same order of magnetude, too;
       * for our case, however, this condition is not met as the residual scales with K_eff*solution, i.e. K_eff deter-
       * mines the order of magnitude; for magnetics, with nu_0 in the range of 1e6, BGM does a bad job;
       *
       * BGM allows the usage of different start matrices B_0 which are better approximations to the Jacobian than the
       * identity I; to include B_0, however, we have to either know B_0_inv or to solve B_0*d = -residual which more
       * or less means, that we have to solve the system using algsys as for the standard solution process;
       * in that sense, BGM is more a solution improvement than an actual algsys-free solution method
       * > it needs to be checked, if BGM can improve the solution such, that the additional BGM computation steps
       * lead to a significant reduction of outer iterations
       */



      // compute solution update according to selected solution approach
      if(currentSolutionApproach_ == SOLVE_NOTSET){
        EXCEPTION("No valid solve flag selected");
      } else if(currentSolutionApproach_ == SOLVE_BGM_ENHANCED){
//        std::cout << "Use BGM enhancer" << std::endl;
        /*
         * new approach: use BGM only as solution improvmenent, i.e. solve system normally using FD-Jacobian (or Fixpoint?!), i.e.
         * get startingResForBGM_editable = z_0 = -B_0_inv*res
         */
        // not quite sure but in the original algorithm, we use the SAME B_0 for each iteration
        // if we solve via Jacobian, would we not get a new B_0 each iteration? maybe this is a problem
        // > try SOLVE_CHORD (i.e. uses the Jacobian from the start of the timestep (until a reeavl is forced > chord method)
//        std::cout << "Step 1: compute z_0 = -B_0_inv*res (standard solution via algsys)" << std::endl;
        startingResForBGM_editable = ComputeSolutionUpdate_ALGSYS(SOLVE_CHORD,steppingTowardsLastTS,steppingTowardsZero,forceJacobian,forceReuseOfJacobian);
        startingResForBGM_editable.ScalarMult(1.0);
//        std::cout << "z_0.Norm() = " << startingResForBGM_editable.NormL2() << std::endl;
        /*
         * pass z_0 as starting vector to BGM
         */
//        std::cout << "Step 2: get improved search direction/solution increment from BGM" << std::endl;
        solIncrement_n_k_ = ComputeSolutionUpdate_BGM(startingResForBGM_editable);
//        std::cout << "retrieved solution increment (NormL2) = " << solIncrement_n_k_.NormL2() << std::endl;
      } else if(currentSolutionApproach_ == SOLVE_BGM){
//        std::cout << "Use pure BGM" << std::endl;
        /*
         * after fixing some major bugs, try pure BGM again
         * > still not working BUT maybe with scaling to eps/mu?
         * > in mangetics (in hyst material) k_eff is in the order of nu_0
         * > scale increment with mu_0
         */
        Double scaling = 1.0;
        if(pdename_ == "magnetic"){
          scaling = -4*M_PI*1e-7;
        }
        // does not work either
        // > test with fixpoint matrix next time
//        std::cout << "scaling: " << scaling << std::endl;
        startingResForBGM_editable = resVec_n_k_;
        startingResForBGM_editable.ScalarMult(scaling);
//        std::cout << "z_0.Norm() = " << startingResForBGM_editable.NormL2() << std::endl;
        /*
         * pass z_0 as starting vector to BGM
         */
//        std::cout << "get search direction/solution increment from BGM" << std::endl;
        solIncrement_n_k_ = ComputeSolutionUpdate_BGM(startingResForBGM_editable);
//        std::cout << "retrieved solution increment (NormL2) = " << solIncrement_n_k_.NormL2() << std::endl;
      } else {
//        std::cout << "preSolve: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;
        solIncrement_n_k_ = ComputeSolutionUpdate_ALGSYS(currentSolutionApproach_,steppingTowardsLastTS,steppingTowardsZero,forceJacobian,forceReuseOfJacobian);
//        std::cout << "postSolve: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;
      }

//      std::cout << "Specieifed non-linear solution method (deprecated): " << nonLinMethod_ << std::endl;
//      std::cout << "Solution approach as specified in input.xml: " << solveFlagToString(solutionApproach_) << std::endl;
//      std::cout << "Currently used solution approach: " << solveFlagToString(currentSolutionApproach_) << std::endl;
//      std::cout << "--- Norm of computed update: " << solIncrement_n_k_.NormL2() << " ---" << std::endl;
//
      LOG_DBG(solvestephyst_main) << "--- Norm of computed update: " << solIncrement_n_k_.NormL2() << " ---";

      if( (steppingTowardsZero)&&(!steppingTowardsLastTS) ){
        // if stepping towards 0 was active due to a reset, check if
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
       * > can be done AFTER solution of system
       */
      solVec_n_kPrev_ = solVec_;

      /*
       * The solution vector will be the increment towards solVecForUpdate_,
       * i.e. the new solVec will be solVecForUpdate + eta*solIncrement_n_k_
       */
      solVec_ = solVecForUpdate_;

      /*
       * Step 5: check solution by evaluating residual and increment
       */
			// unlock hyst operator (if locked)
			PDE_.SetFlagInCoefFncHyst("lockHystOperator",0);

      // before we evaluate the hyst operator with the actual solution,
      // store its old value for next iteration
      // > do this BEFORE linesearch as operators get evaluated during LS, too
      // Note: by calling SetPreviousHystVals with 0 (=> lastTS = false),
      //        coefFncHyst will evaluate the hysteresis operator with
      //        the current state of solVec_; for this evaluation, the memory of
      //        the hysteresis operator will be locked, i.e. the evaluation will
      //        lead to reversible changes
      //      > NEW VERSION: take the preevaluated values and copy them to the corresponding
      //              arrays; this has to be done before evaluating the final result for the iteration
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);

      matrixFlag flag_for_systemMatrix;
      matrixFlag flag_for_resEvalMatrix;
      vectorFlag flag_for_rhsVector;
      vectorFlag flag_for_resVector;
      bool estimatorUseful;

//      std::cout << "preLinesearch: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;
      Helper_GetParamsFromSolveFlag(currentSolutionApproach_, flag_for_systemMatrix, flag_for_resEvalMatrix, flag_for_rhsVector, flag_for_resVector, estimatorUseful);

      // skip setting of E_H_ and P_H_ in coefFunctionHyst during the linesearch process
      // reason: for eastimating the input to the hysteresis operator (for evaluation purposes), the H-version
      // of the FP method access these two arrays; in consequence, a repeated evaluation of the hyst operator will
      // more and more alter the estimates which is NOT correct
      // note that this flag only has a function if useFPH_ is set to true in CoefFunctionHyst
      PDE_.SetFlagInCoefFncHyst("SkipStorage",1);
      
      //LOG_DBG(solvestephyst_main) << "SetupRESVector/Perform LS";
//      std::cout << "setIDBC? " << setIDBC << std::endl;
      if( setIDBC || ( lineSearch_ == "none" ) || specifiedLinesearches_.empty() ){
        LOG_DBG(solvestephyst_main) << "--- Perform full step ---";
        /*
         * Linesearch is only allowed if the increment does not hold contributions
         * of the Dirichlet boundary as in that case, we have to step with length 1
         * to get the boundary terms correct!
         * > during first iteration for incremental towards last iteration
         * > during each iteration for incremental towards last ts or 0
         */
        etaLinesearch = 1.0;
        currentLinesearchApproach_ = NO_LINESEARCH;
      } else {
        LOG_DBG(solvestephyst_main) << "--- Perform linesearch ---";
        etaLinesearch = LineSearchHyst(flag_for_resEvalMatrix, flag_for_resVector, solVec_, solIncrement_n_k_,iterationCounterLog);
        // for testing: overwrite etaLinesearch by 1.0 i.e. full step; now the any LS version should
        // result in the same steps as the version without LS; if this is not the case,
        // linesearchHyst has some unwanted effects on the system! > seems to work
//        etaLinesearch = 1.0;
      }
//      std::cout << "postLinesearch: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;

      if((etaLinesearch == 0)&&(forceJacobian == false)){
        // linesearch was not successful; solve system again but this time
        // with different approach > full jacobian
        forceJacobian = true;
        LOG_DBG(solvestephyst_main) << "-!!!- linesearch was not successful; discard update and force Jacobian for next iteration -!!!-";
        collectedLoggingMSG << "-!!!- DISCARD UPDATE -!!!-" << std::endl;
        // decrease counter; we track only full iterations
        iterationCounterLog--;
        iterationCntCurrentTS_--;
//        std::cout << "Discard update" << std::endl;
        continue;
      } else {
        forceJacobian = false;
      }

      /*
       * NEW SOLUTION
       */
      solIncrement_n_k_.ScalarMult(etaLinesearch);
      solVec_.Add(1.0,solIncrement_n_k_);
      // reset for next iteration
      linSysMatrix_n_0_TIMES_solVec_n_kPrev_precomputed_ = false;
      linSysMatrix_n_0_TIMES_solIncrement_n_k_precomputed_ = false;
      // allow storage of E_H_ and P_H_ again
      PDE_.SetFlagInCoefFncHyst("SkipStorage",0);
      
      /*
       * For incremental error checking, compute the actual update between
       *  the current and the previous iteration
       */
//      solIncrement_n_k_ = solVec_;
//      solIncrement_n_k_.Add(-1.0,solVec_n_kPrev_);

//      incrementalErrAbs = solIncrement_n_k_.NormL2();
//      solL2Norm  = solVec_.NormL2();

//      if ( solL2Norm > 0.0 ){
//        incrementalErrRel = incrementalErrAbs/solL2Norm;
//      } else {
//        incrementalErrRel = incrementalErrAbs;
//      }

      /*
       * New solution has been found
       * > evaluate hyst operator
       */
      // allow setting of matrices
      // NOTE: delta matrices will be set with values from last iteration so make
      // sure that SetPreviousHystValues was called before
      PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
      // EvaluateHystOperators; flag = 1 > set Matrix for Inverion, too
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",1);

      //LOG_DBG(solvestephyst_main) << "Updated solution: " << solVec_.ToString();
      // although LineSearchHyst will setup the residual, we have to
      // call this function once more with the actual best eta to ensure
      // that the system is assemble correctly!
      //debugOutput_ = true;
      LOG_DBG(solvestephyst_main) << "--- Evaluate current residual ---";

//      normOfSystem =
//      std::cout << "preRes: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;
      SetupRESorRHS(flag_for_resEvalMatrix, flag_for_resVector, solVec_, solVec_, solVec_, resVec_n_k_, false);
//      normOfSystem = SetupRESVector(resVec_n_k_,solVec_,solVec_,solVec_,iterationCounterLog);

      /*
       * NEW checks
       * a) absolute norm of actual residual
       * b) relative norm of actual residual wrt res0
       * c) relative norm of actual residual wrt previous residual
       */
//      //a)
//      residualErrAbs = resVec_n_k_.NormL2();
//      //b)
//      residualErrRel = residualErrAbs/resVec_n_0_.NormL2();
//
//      // a)
//      residualErrAbs = resVec_n_k_.NormL2();
//      residualErrRel = residualErrAbs/resVec_n_0_.NormL2();
//
      // b)
      resIncrement_n_k_ = resVec_n_k_;
      resIncrement_n_k_.Add(-1.0,resVec_n_kPrev_);
//      Double resDiffNorm = resIncrement_n_k_.NormL2();

      // c)
//      Double resQuotNorm = residualErrAbs/resVec_n_kPrev_.NormL2();

//      Double residualErrRelAlt = residualErrAbs;

//      std::string tag = "";
//      std::string tag2 = "";
//      if ( normOfSystem > 0.0 ){
//        residualErrRel = residualErrAbs/normOfSystem;
//        resDiffNorm = resDiffNorm/normOfSystem;
//        tag = "rel";
//      } else {
//        residualErrRel = residualErrAbs;
//      }
//
//      if(resVec_n_0_.NormL2() != 0){
//        tag2 = "rel";
//        residualErrRelAlt = residualErrAbs/resVec_n_0_.NormL2();
//      }
//
//      residualErrRel = residualErrRelAlt; // take residual error relative to initial residual
//
      if(setStartingResForBGM){
        // take current residual as base for BGM update
        // i.e.
        // store current residual and reset incremenets!
        setStartingResForBGM = false;
        startingResForBGM = resVec_n_k_;
        historyLength_ = 0;
      }

      if(storeIncrements){
        UInt currentPosition = historyLength_%hisotryLengthMax_;
//        std::cout << "Current position: " << currentPosition << std::endl;
        resIncrements_history_[currentPosition] = resIncrement_n_k_;
        solIncrements_history_[currentPosition] = solIncrement_n_k_;
        etaLinsearch_history_[currentPosition] = etaLinesearch;
        historyLength_++;

        // testing
//        Vector<UInt> indices = Helper_GetStorageIndices();
        //        std::cout << "Indices: " << indices.ToString() << std::endl;
      }

//      residualErrAbsPrev_ = residualErrAbs;

//      std::stringstream minLoggingMSG;
//      if(iterationCounterLog >= 1){
//        minLoggingMSG << "Iteration " << iterationCounterLog << " (timestep "<<timestepCnt_<<"): "<< solveFlagToString(currentSolutionApproach_)<<std::endl;
//        minLoggingMSG << "etaUsed = " << etaLinesearch << "(" << linesearchFlagToString(currentLinesearchApproach_) << ")" << std::endl;
//        minLoggingMSG << " -> " << "rel norm of increment |solInc_" << iterationCounterLog << "|:" << incrementalErrRel << std::endl;
//        minLoggingMSG << " -> " << "abs norm of residual |r_full_" << iterationCounterLog << "|:" << residualErrAbs << std::endl;
//        //        if(tag == "rel"){
//        //          minLoggingMSG << " -> " << tag << " norm of residual |r_full_" << iterationCounterLog << "| wrt. system Norm:" << residualErrRel << std::endl;
//        //        }
//        if(tag2 == "rel"){
//          minLoggingMSG << " -> " << tag << " norm of residual |r_full_" << iterationCounterLog << "| wrt. initialResidual:" << residualErrRelAlt << std::endl;
//        }
//      }

//      if(!steppingTowardsLastTS){
//        if(((resDiffNorm < minResDecreaseTillReset_))&&(iterationCounterLog >= minNumberItTillReset_)){
//          forceReset++;
//        }
//      }

      resVec_n_kPrev_ = resVec_n_k_;

//      Double incrementalErrCheck;
//      Double residualErrCheck;

//      /*
//       * TODO: adapt residual criteraion for the case of constant excitation!
//       * in this the initial residual will already be very small such that the relative residual will be larger than
//       * tolerance even! > switch to absolute residual or ...
//       */
//
//      if(useRelativeNormForRes_){
//        residualErrCheck = residualErrRel;
//      } else {
//        residualErrCheck = residualErrAbs;
//      }

//      if(residualErrCheck < bestResNorm){
//        bestResNorm = residualErrCheck;
//        bestSol_ = solVec_;
//      }

//      resNormHistory.push_back(residualErrCheck);
//      resNormHistory.pop_front();
//
//      if(useRelativeNormForInc_){
//        incrementalErrCheck = incrementalErrRel;
//      } else {
//        incrementalErrCheck = incrementalErrAbs;
//      }
//
      // reset logger
      convergenceLogger.str("");
      convergenceLogger << std::endl;
      int iterationSuccess = checkStoppingCriteria(iterationCounterLog, convergenceLogger, 1);

//      std::cout << "logger: currentSolutionApproach_ = " << solveFlagToString(currentSolutionApproach_) << std::endl;
      convergenceLogger << "- Used method: " << solveFlagToString(currentSolutionApproach_) << std::endl;
      convergenceLogger << "- Used linesearch: " << linesearchFlagToString(currentLinesearchApproach_) << std::endl;
      convergenceLogger << "- Used step-length: " << etaLinesearch << std::endl;

      if(minLoggingToTerminal_){
        std::cout << convergenceLogger.str() << std::endl;
      }

//      performOneMoreStep =
//              (incrementalErrCheck > incStopCrit_) || (residualErrCheck > residualStopCrit_);
//
      if(iterationSuccess == -1){
        performOneMoreStep = true;
      } else if(iterationSuccess >= 0){
        performOneMoreStep = false;
      }

      bool resetSolVecToZero = false;
      if(((iterationCounterLog >= iterationsTillReset_)&&(numResets == 0))&&(iterationsTillReset_ != 0)){
        resetSolVecToZero = true;
      }
      if(resetSolVecToZero){
        std::cout << "++ " << iterationsTillReset_ << " iterations exceeded; reset solution vector to 0" << std::endl;
//      // sometimes we cannot converge from the current solution; in that case we might
        // try a reset of the current solution vector to 0 first and step from that point

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

        if(steppingTowardsLastTS){
          PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);
        } else {
          PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
        }

        steppingTowardsZero = true;
        numResets++;
        iterationCounterLog = 0;
        iterationOfLatestReset_ = iterationCntCurrentTS_;
        historyLength_ = 0;

        convergenceLogger << "- RESET Solution vector to 0 -" << std::endl;

//        if(minLoggingToTerminal_){
//          std::cout << "-!!!- RESET -!!!-" << std::endl;
//        }
      }

//
//      // compare residual at front and back of deque;
//      // residual should have decreased (at least a bit); otherwise try reset
//      if(( (performOneMoreStep) && (resNormHistory.front()-minResDecreaseTillReset_ <= resNormHistory.back()) ) || (forceReset==1)
//              || (performOneMoreStep && (iterationCounterLog >= nonLinMaxIter_)) ){
//
//        LOG_DBG2(solvestephyst_main) << "-!!!- Decrease in residual not sufficient / check failback -!!!-";
//
//        forceReset++;
//        // check for failback
////        if(bestResNorm < failBackCrit_){
////          if(minLoggingToTerminal_){
////            minLoggingMSG << "-!!!- FAILBACK - bestResNorm = " << bestResNorm << " < failBackCrit_ = " << failBackCrit_ << " - !!!-" << std::endl;
////          }
////
////          solVec_ = bestSol_;
////          residualErrAbs = bestResNorm;
////
////          performOneMoreStep = false;
////        }
//
//        // failback not met, check for possible reset
//        if((numResets == 0)&&allowReset_){
//
//          LOG_DBG2(solvestephyst_main) << "-!!!- Decrease in residual not sufficient / perform reset -!!!-";
//
//          nonLinMethod_ = nonLinMethodAfterReset_;
//          towardsLastIteration_ = towardsLastIterationReset_;
//          deltaMatCoupling_ = deltaMatCouplingReset_;
//          deltaMatStrain_ = deltaMatStrainReset_;
//          lineSearch_ = lineSearchAfterReset_;
//
//          LOG_DBG2(solvestephyst_main) << "---- Switch to  " << nonLinMethod_ << " with LS " << lineSearch_;
//
//          if((nonLinMethod_ == "HYST_DeltaMat_IT")||(nonLinMethod_ == "HYST_FixPoint_IT")){
//            steppingTowardsLastTS = false;
//          } else {
//            steppingTowardsLastTS = true;
//          }
//
//          // important: we have to evaluate the hyst operator with the reset solution value
//          // to be able to compute correct deltaMatrices
//          // for this evaluation, we have to take care of:
//          // a) correct storage
//          //      for stepping towards lastIT, we have to set the storage for lastIT
//          //      for stepping towards lastTS, we have to set the storage for lastTS
//          //  > this is adjusted by the flag lastTS (as used prior)
//          // b) memory of hyst operator must not be overwritten (we do not want up with
//          //    some minor loops due to iterations)
//          //  > force memory lock (required if lastTS = true)
//          solVec_.Init();
//
//          // sol vec was reset > evaluate hyst operator
//          // + reset matrix for inversion
//          //      > here with flag 2 > set solution to 0 first, then set matrix, then eval hyst operator
//          // allow setting of matrices
//          PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
//          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",2);
//
//          lastTS = steppingTowardsLastTS;
//
//          if(lastTS){
//            PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);
//          } else {
//            PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",0);
//          }
//          //bool forceMemoryLock = true;
//          //PDE_.SetPreviousHystVals(lastTS,forceMemoryLock);
//
//          //PDE_.EstimateCurrentSlopeForHysteresis(steppingLength, scaling);
//          //LOG_DBG(solvestephyst_main) << "Iteration " << iterationCounter << " of timestep " << timestepCnt_;
//
//          // reset deque
//          for(UInt i = 0; i < resNormHistory.size(); i++){
//            resNormHistory[i] = 50000;
//          }
//          steppingTowardsZero = true;
//          numResets++;
//          iterationCounterLog = 0;
//          iterationOfLatestReset_ = iterationCntCurrentTS_;
//          historyLength_ = 0;
//
//          if(minLoggingToTerminal_){
//            minLoggingMSG << "-!!!- RESET -!!!-" << std::endl;
//          }
//        }
//      }

      /*
       * Logging
       */
//      if(minLoggingToTerminal_){
//        collectedLoggingMSG << minLoggingMSG.str();
//      }
      LOG_DBG(solvestephyst_main) << "--- Convergence information for iteration " << iterationCntCurrentTS_ << " (" << totalIterationCnt_ << ") of time step " << timestepCnt_ << " ---";
      LOG_DBG(solvestephyst_main) << convergenceLogger.str();

      // output of norms and data to info.xml
      if ( nonLinLogging_ == true ) {
        // get current step
        UInt actStep = PDE_.GetSolveStep()->GetActStep();

        if (PDE_.IsIterCoupled() == false) {
          couplingIter_ = 0;
        }

        Output_WriteHystIterToInfoXML(pdename_,
                solveFlagToString(currentSolutionApproach_),
                linesearchFlagToString(currentLinesearchApproach_),
                couplingIter_,
                actStep,
                iterationCntCurrentTS_,
                convergenceLogger.str());
//                iterationCounterLog,
//                iterationCntCurrentTS_,
//                residualErrAbs, residualErrRel, useRelativeNormForRes_,
//                incrementalErrAbs, incrementalErrRel, useRelativeNormForInc_,
//                etaLinesearch, reset, failback);
      }

//      if(minLoggingToTerminal_){
//        std::cout << minLoggingMSG.str() << std::endl;
//      }

    } while(performOneMoreStep && (iterationCounterLog < nonLinMaxIter_));

//    LOG_TRACE(solvestephyst_main) << "== Collected convergence information for time step " << timestepCnt_ << " ==";
//    LOG_TRACE(solvestephyst_main) << collectedLoggingMSG.str();

//    if(minLoggingToTerminal_ == 2){
//      std::cout << "== Collected convergence information for time step " << timestepCnt_ << " ==" << std::endl;
//      std::cout << collectedLoggingMSG.str() << std::endl;
//    }

    if(iterationCntCurrentTS_ > maxIterationPerTSCnt_){
      maxIterationPerTSCnt_ = iterationCntCurrentTS_;
    }

    if(minLoggingToTerminal_ == 2){
      std::cout << "== Iteration information for time step " << timestepCnt_ << " ==" << std::endl;
      std::cout << "- number of non-linear iterations: " << iterationCntCurrentTS_ << std::endl;
      std::cout << "- max number of non-linear iterations for a single TS (sofar): " << maxIterationPerTSCnt_ << std::endl;
      std::cout << "- total number of non-linear iterations for all TS (sofar): " << totalIterationCnt_ << std::endl;
    }

    if (performOneMoreStep && (iterationCntCurrentTS_ >= nonLinMaxIter_)){
      converged = false;
    } else {
      converged = true;
    }

    if(!converged){
      WARN("NON CONVERGENCE error in PDE '" << pdename_
              << "' in step no '" << PDE_.GetSolveStep()->GetActStep()
              << "' at iteration '" << iterationCntCurrentTS_ << "'"
              << convergenceLogger.str());
//              << "'.\n ==> incremental error (abs): " << incrementalErrAbs
//              << "'.\n ==> incremental error (rel): " << incrementalErrRel
//              << "\n ==> residual error (abs): " << residualErrAbs
//              << "\n ==> residual error (rel): " << residualErrRel);
    }

    return converged;
  }

  int SolveStepHyst::checkStoppingCriteria(UInt iterationCounter, std::ostringstream& logger, UInt loggingLevel){
    logger << "Convergence results for iteration " << iterationCounter << std::endl;

    /*
     * Check standard stopping criteria first;
     * if ALL are passed > return 0 (SUCCESS!)
     */
    std::set< StoppingCriterion >::iterator criterionIterator;
    bool allPassed = true;
    Double valueToCheck = 0.0;
    Double precision = 1e-12;

    logger << "- Checked stopping criteria: " << std::endl;
    for(criterionIterator = stoppingCriteria_.begin(); criterionIterator != stoppingCriteria_.end(); criterionIterator++){
      // shall criterion be checked already?
      if(criterionIterator->firstCheck_ > iterationCounter){
        continue;
      }

      // is it time for a recheck?
      int needsRecheck = (iterationCounter - criterionIterator->firstCheck_)%criterionIterator->checkingFrequency_;
//      std::cout << needsRecheck << std::endl;
      if(needsRecheck != 0){
        continue;
      }

      // do we check residual or increment?
      if(criterionIterator->isResidual_){
        valueToCheck = resVec_n_k_.NormL2();

        // relative w.r.t. initial residual?
        if(criterionIterator->isRelative_){
          if(resVec_n_0_.NormL2() > precision){
            valueToCheck /= resVec_n_0_.NormL2();
          } else {
            logger << "Info: Norm of initial residual (" << resVec_n_0_.NormL2() << ") < " << precision << "; taking absolute value of residual" << std::endl;
          }
        }

      } else {
        valueToCheck = solIncrement_n_k_.NormL2();

        // relative w.r.t. L2-norm of actual solution vector?
        if(criterionIterator->isRelative_){
          if(solVec_.NormL2() > precision){
            valueToCheck /= solVec_.NormL2();
            if(minLoggingToTerminal_ > 1){
              logger << "Info: Norm of current solution vector = " << solVec_.NormL2() << std::endl;
            }
          } else {
            logger << "Info: Norm of current solution vector (" << solVec_.NormL2() << ") < " << precision << "; taking absolute value of increment" << std::endl;
          }
        }
      }

      // additional scaling (e.g., to number of unknowns)?
      if( (criterionIterator->scaleToSystemSize_) && (systemSize_ != 0) ){
        valueToCheck /= systemSize_;
      }

      if(valueToCheck > criterionIterator->stoppingTolerance_){
        allPassed = false;
        logger << " (" << criterionIterator->checkingOrder_ << "): " << criterionIterator->nameTag_ << " \t FAILED";
        if(loggingLevel > 0){
          logger << " ( " << valueToCheck << " > " << criterionIterator->stoppingTolerance_ << " ) " << std::endl;
        }
      } else {
        logger << " (" << criterionIterator->checkingOrder_ << "): " << criterionIterator->nameTag_ << " \t PASSED";
        if(loggingLevel > 0){
          logger << " ( " << valueToCheck << " <= " << criterionIterator->stoppingTolerance_ << " ) " << std::endl;
        }
      }

      // store for later
      criterionIterator->lastCheckedValue_ = valueToCheck;
    }

    if(allPassed){
      logger << " -> ALL convergence criteria satisfied -> STOP " << std::endl;
      return 0;
    }

    bool anyPassed = false;
    bool firstFailbackLog = true;

    if(failbackCriteria_.size() == 0){
      return -1;
    }

//    std::cout << "Checking " << failbackCriteria_.size() << " failback criteria " << std::endl;
    for(criterionIterator = failbackCriteria_.begin(); criterionIterator != failbackCriteria_.end(); criterionIterator++){
      // shall criterion be checked already?
      if(criterionIterator->firstCheck_ > iterationCounter){
//        std::cout << "no check for " << criterionIterator->nameTag_ << std::endl;
//        std::cout << "Reason: firstCheck = " << criterionIterator->firstCheck_ << " but current it = " << iterationCounter << std::endl;
        continue;
      }

      // is it time for a recheck?
      int needsRecheck = (iterationCounter - criterionIterator->firstCheck_)%criterionIterator->checkingFrequency_;
      if(needsRecheck != 0){
//        std::cout << "no recheck for " << criterionIterator->nameTag_ << std::endl;
//        std::cout << "Reason: iterationCounter - criterionIterator->firstCheck_ = " << iterationCounter - criterionIterator->firstCheck_ << "; criterionIterator->checkingFrequency_ = " << criterionIterator->checkingFrequency_ << std::endl;
//        std::cout << "(iterationCounter - criterionIterator->firstCheck_)%criterionIterator->checkingFrequency_ = " << (iterationCounter - criterionIterator->firstCheck_)%criterionIterator->checkingFrequency_ << std::endl;
        continue;
      }

      if(firstFailbackLog){
        logger << "- Checked failback criteria: " << std::endl;
        firstFailbackLog = false;
      }

      // do we check residual or increment?
      if(criterionIterator->isResidual_){
        valueToCheck = resVec_n_k_.NormL2();

        // relative w.r.t. initial residual?
        if(criterionIterator->isRelative_){
          if(resVec_n_0_.NormL2() > precision){
            valueToCheck /= resVec_n_0_.NormL2();
          } else {
            logger << "Info: Norm of initial residual (" << resVec_n_0_.NormL2() << ") < " << precision << "; taking absolute value of residual" << std::endl;
          }
        }

      } else {
        valueToCheck = solIncrement_n_k_.NormL2();

        // relative w.r.t. L2-norm of actual solution vector?
        if(criterionIterator->isRelative_){
          if(solVec_.NormL2() > precision){
            valueToCheck /= solVec_.NormL2();
          } else {
            logger << "Info: Norm of current solution vector (" << solVec_.NormL2() << ") < " << precision << "; taking absolute value of increment" << std::endl;
          }
        }
      }

      // additional scaling (e.g., to number of unknowns)?
      if( (criterionIterator->scaleToSystemSize_) && (systemSize_ != 0) ){
        valueToCheck /= systemSize_;
      }

      if(valueToCheck > criterionIterator->stoppingTolerance_){
        logger << " (" << criterionIterator->checkingOrder_ << "): " << criterionIterator->nameTag_ <<  " \t FAILED";
        if(loggingLevel > 0){
          logger << " ( " << valueToCheck << " > " << criterionIterator->stoppingTolerance_ << " ) " << std::endl;
        }
      } else {
        anyPassed = true;
        logger << " (" << criterionIterator->checkingOrder_ << "): " << criterionIterator->nameTag_ <<  " \t PASSED";
        if(loggingLevel > 0){
          logger << " ( " << valueToCheck << " <= " << criterionIterator->stoppingTolerance_ << " ) " << std::endl;
        }
      }

      // store for later
      criterionIterator->lastCheckedValue_ = valueToCheck;

      // failback satisfied; stop
      if(anyPassed){
        logger << " -> ANY failback criterion satisfied -> STOP " << std::endl;
        return 1;
      }
    }

    return -1;

  }

  void SolveStepHyst::Solve_Step3_FinishSolveStep(){
//    std::cout << "Finish Solve Step" << std::endl;

    LOG_TRACE(solvestephyst_main) << "===================================";
    LOG_TRACE(solvestephyst_main) << "== Solve_Step3_FinishSolveStep  ===";
    LOG_TRACE(solvestephyst_main) << "===================================";

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
    std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();++fncIt){
      fncIt->second->GetTimeScheme()->FinishStep();
    }

    // recover solution strategy
    nonLinMethod_ = nonLinBackup_;
//    towardsLastIteration_ = towardsLastIterationBackup_;
    deltaMatCoupling_ = deltaMatCouplingBackup_;
    deltaMatStrain_ = deltaMatStrainBackup_;
    lineSearch_ = lineSearchBackup_;

    /*
     * solution has been found > set values AND overwrite memory
     * (-1 > same as +1 but additionally overwrite hyst memory
     */
//    std::cout << "Finish Solve Step -- allowSettingOfMatForLocalInversion" << std::endl;
    PDE_.SetFlagInCoefFncHyst("allowSettingOfMatForLocalInversion",true);
//    std::cout << "Finish Solve Step -- EvaluateHystOperators" << std::endl;
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",-1);

    /*
     * STORE CURRENT VALUES FOR NEXT Time step
     */
    // Note: by calling SetPreviousHystVals with flag lastTS = true,
    //        coefFncHyst will evaluate the hysteresis operator with
    //        the final state of solVec_; for this evaluation, the memory of
    //        the hysteresis operator will be unlocked, i.e. the evaluation will
    //        lead to irreversible changes
    //PDE_.SetPreviousHystVals(lastTS);
    // do this AFTER the final value was computed
//    std::cout << "Finish Solve Step -- SetPreviousHystValues" << std::endl;
    PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",1);

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
  }

  /*
   * IV. Methods for computing solution updates
   * > called during each iteration
   * > inner loop
   * > k = iteration counter
   */

  // TODO: test
  SBM_Vector SolveStepHyst::ComputeSolutionUpdate_BGM(SBM_Vector& z0){
    /*
     * Compute solution update via Globalized Broydens (Good) Method
     * source:
     *  "Iterative Methods for Linear and NOnlinear Equations" - C.T.Kelly - Siam 1995
     *  Algorithm 8.3.1.; p. 145; chapter 8
     *
     * some notes:
     *  Broydens Method can be used to solve nonlinear, nonsymmetric systems of equationts;
     *  Convergence is only assured locally and only in case of Lipschitz continuous Jacobian.
     *  > we need global convergence and doe not have Lipschitz continuity as Preisach hyst
     *    model has no continous derivative
     *  However, it might still work (we just cannot ensure it). For globalization, we
     *  combine Broydens Method with a linesearch. In case of a recursive, matrix free
     *  implementation, like the one which is used here, the recursion/iteration formula
     *  is affected by the current and past linesearch steps, i.e. do not simply
     *  use the standard recursion formula (as in chapter 7 of aboves source).
     *  The matrix form should be the same, however.
     *
     * Notes:
     * - use matrix-free, recursive formula
     * - adapt formula to allow for linesearch
     *
     * Additional notes/adaptions to Algorithm 8.3.1 - Step 2e) and 2f)
     * - 2e) source implicitely assumes B_0 = Identity, i.e. it solves
     *        B_0*z = -F(x) = residual
     *    > however, Identity might not be so good for our purpose, as actual Jacobian
     *        is quite different
     *    > allow two modifications:
     *      i) Provide z = -B_0^-1*F(x) as input
     *      or
     *      ii) scale Identity with nu_0 / eps_0, such that
     *        z = -mu_0*F(x) / -1/eps_0*F(x)
     *    > i) should be prefered as this allows for the usage of a FD Jacobian
     *          approximation of B_0 AND it allows to incorporate IDBC
     *        -> during first iteration (or after reset), solve B_0*z0 = -F(x) using algsys,
     *            during later iterations, use BGM
     *
     * NOTE 12.03.2019
     * BGM does not work at all; after reading Kelley again it seems to be clear why.
     * According to Algorithm 8.3.1 z = -F(x)   (as already stated above)
     * In case of B_0 != I, they use B_0 as preconditioner directly inside of F(x), i.e.
     * their F(x) is actually B_0^-1 F(x)
     * In other words, each time F(x) appears in their algorithm, it has to be replaced by B_0^-1 F(x) for our
     * purpose > it is not sufficient to compute B_0*z_0 = -F(x) during the FIRST iteration but we have to do this
     * in EACH iteration.
     * > needs testing
     * Another possibly better algorithm which can be found in Deuflhard, Freund, Walter "Solution of Large NonSysmmetric Linear Systems"
     * that was used by Kelley for their version. Note that Deulhard et al. have the algorithm for linear system, but it should
     * be extendable to non-linear systems as well.
     *
     */

    // we assume that 2e)i) is already computed, either from solving via algsys or
    // by setting z0 = scaling*residual
    // compute new d using recursive formula, i.e. 2e)ii)
    SBM_Vector newUpdate;
    newUpdate = zeroVec_;

    Double a,b,c;
    Double sjTz;
    Double sjTsj;

    Vector<UInt> indices = Helper_GetStorageIndices();
    UInt numIndices = indices.GetSize();

//    std::cout << "Pre iteation: " << std::endl;
//    std::cout << "z0.NormL2() = " << z0.NormL2() << std::endl;
//    std::cout << "Num indices: " << numIndices << std::endl;

    UInt j = 0;
    UInt jpp = 0;
    for(UInt i = 0; i < numIndices-1; i++){
//      std::cout << "Iteration i = " << i << std::endl;
      j = indices[i];
      jpp = indices[i+1];
//      std::cout << "etaLinsearch_history_[j] = " << etaLinsearch_history_[j] << std::endl;
//      std::cout << "etaLinsearch_history_[jpp] = " << etaLinsearch_history_[jpp] << std::endl;
      a = -etaLinsearch_history_[j]/etaLinsearch_history_[jpp];
      b = 1.0 - etaLinsearch_history_[j];

//      std::cout << "a: " << a << std::endl;
//      std::cout << "b: " << b << std::endl;

      sjTsj = solIncrements_history_[j].NormL2();
      sjTsj *= sjTsj;
      solIncrements_history_[j].Inner(z0,sjTz);

//      std::cout << "sjTsj = " << sjTsj << std::endl;
//      std::cout << "sjTz = " << sjTz << std::endl;

      c = sjTz/sjTsj;
//      std::cout << "c = sjTz/sjTsj: " << c << std::endl;
//      std::cout << "a*c = " << a*c << std::endl;
//      std::cout << "b*c = " << b*c << std::endl;
      // NOTE: do not use Add(c1,vec1,c2,vec2)!
      //      this will not add c1*vec1 and c2*vec2 to z0 but will replace z0 with the sum!
//      z0.Add(a*c,solIncrements_history_[jpp],b*c,solIncrements_history_[j]);
      z0.Add(a*c,solIncrements_history_[jpp]);
      z0.Add(b*c,solIncrements_history_[j]);
//      std::cout << "solIncrements_history_[jpp].NormL2() = " << solIncrements_history_[jpp].NormL2() << std::endl;
//      std::cout << "solIncrements_history_[j].NormL2() = " << solIncrements_history_[j].NormL2() << std::endl;
//      std::cout << "z0.NormL2() = " << z0.NormL2() << std::endl;
    }
    // 2f)
    j = indices[numIndices-1];
    sjTsj = solIncrements_history_[j].NormL2();
    sjTsj *= sjTsj;

    solIncrements_history_[j].Inner(z0,sjTz);

//    std::cout << "After iteation: " << std::endl;
//    std::cout << "z0.NormL2() = " << z0.NormL2() << std::endl;
//    std::cout << "n = indices[numIndices-1]: " << j << std::endl;
//    std::cout << "etaLinsearch_history_[n]: " << etaLinsearch_history_[j] << std::endl;
//    std::cout << "snTsn = " << sjTsj << std::endl;
//    std::cout << "(1-etaLinsearch_history_[n]): " << (1-etaLinsearch_history_[j]) << std::endl;
//    std::cout << "(1-etaLinsearch_history_[n])*snTz: " << (1-etaLinsearch_history_[j])*sjTz << std::endl;
//    std::cout << "snTsn - etaLinsearch_history_[n]*snTz: " << sjTsj - etaLinsearch_history_[j]*sjTz << std::endl;
    // here we can use add with two vectors, as newUpdate is zero up to this point
    newUpdate.Add(sjTsj,z0,(1-etaLinsearch_history_[j])*sjTz,solIncrements_history_[j]);
    newUpdate.ScalarDiv(sjTsj - etaLinsearch_history_[j]*sjTz);

//    currentSolutionApproach_ = SOLVE_BGM;
    return newUpdate;
  }

  SBM_Vector SolveStepHyst::ComputeSolutionUpdate_ALGSYS(solveFlag flagForSolving,
          bool forceStepTowardsPrevTS, bool forceFullStep, bool forceJacobian, bool forceReuseOfJacobian){
    /*
     * Compute solution update via ALGSYS (standard CFS solution approach)
     * i.e. solve
     *
     *  K_eff(u_k) du = res(u_k)
     *
     * with
     *  K_eff(u_k,u_PREV1)    ... effective stiffness matrix evaluated at current solution u_k and
     *                                          possibly different previous solutions u_PREV1
     *  du = u_k+1 - u_PREV2  ... solution increment towards u_PREV2
     *
     *  res(u_k) = linRHS + nonLinRHS(u_k) + hystereticRHS(u_PREV3) - K_eff_rhs(u_k)*u_PREV2
     *                        ... current residual vector consisting of current parts (u_k) and parts
     *                              corresponding to u_PREV3; K_eff_rhs = K_eff but without any approximations to the Jacobian
     *
     * > K_eff, K_eff_rhs, u_PREV1-3 depend on the solveFlag
     *
     * SOLVE_GLOBAL_FIXPOINT_B. SOLVE_LOCAL_FIXPOINT_B, SOLVE_GLOBAL_FIXPOINT_H:
     *  u_PREV1-3 = u_k
     *  K_eff = (non-)linear system matrix using SCALED material tensors evaluated at current iterate u_k
     *  K_eff_rhs = K_eff
     * > note: due to nu_FP != nu_0, we have to subtract (nu_0-nu_FP)*B from output of hyst operator > see CoefFunctionHyst
     *
     * SOLVE_JACOBIAN:
     *  u_PREV1-3 = u_k
     *  K_eff = Approximation to Jacobian of current residual, calculated via Finite Differences on element level
     *  K_eff_rhs = (non-)linear system matrix using standard material tensors evaluated at current iterate u_k
     *
     * SOLVE_DELTAMAT_IT
     *  u_PREV1 = u_k-1
     *  u_PREV2 = u_k
     *  u_PREV3 = u_k (should be u_k-1 from theory but this does not converge properly)
     *  K_eff = Approximation to Jacobian of current residual, calculated via deltaMat approach on element level using
     *              delta between u_k and uPrev1
     *  K_eff_rhs = (non-)linear system matrix using standard material tensors evaluated at current iterate u_k
     *
     * SOLVE_DELTAMAT_TS
     *  u_PREV1 = u_0
     *  u_PREV2 = u_0
     *  u_PREV3 = u_0
     *  K_eff = Approximation to Jacobian of current residual, calculated via deltaMat approach on element level using
     *              delta between u_k and uPrev1
     *  K_eff_rhs = (non-)linear system matrix using standard material tensors evaluated at current iterate u_k
     *
     * SOLVE_DELTAMAT_TS_TOWARDS_IT
     *  u_PREV1 = u_0
     *  u_PREV2 = u_k
     *  u_PREV3 = u_k (should be u_0 from theory but this does not converge properly)
     *  K_eff = Approximation to Jacobian of current residual, calculated via deltaMat approach on element level using
     *              delta between u_k and uPrev1
     *  K_eff_rhs = (non-)linear system matrix using standard material tensors evaluated at current iterate u_k
     */

    /*
     * Initial phase:
     *  setup variables
     *  determine kind of solution update
     *  get flags for rhs and system assembly
     */
    SBM_Vector newUpdate;
    newUpdate = zeroVec_;

    matrixFlag flag_for_systemMatrix;
    matrixFlag flag_for_resEvalMatrix;
    vectorFlag flag_for_rhsVector;
    vectorFlag flag_for_resVector;
    bool estimatorUseful;
    Helper_GetParamsFromSolveFlag(flagForSolving, flag_for_systemMatrix, flag_for_resEvalMatrix, flag_for_rhsVector, flag_for_resVector, estimatorUseful);

    /*
     * Determine kind of solution update
     * default: solution update = incremenet towards previous iterate
     * steppingTowardsLastTS = true: soluiton update = increment towards previous time step
     * steppingTowardsZero = true: solution update = full solution
     */
    bool steppingTowardsLastTS = false;
    if((flag_for_rhsVector == PREVIOUS_TS_RES)||forceStepTowardsPrevTS){
      steppingTowardsLastTS = true;
    }

    bool steppingTowardsZero = false;
    if(forceFullStep){
      steppingTowardsZero = true;
    }

    /*
     * SteppingTowardsZero > full step; solution increment = full solution
     * SteppingTowardsLastTS > incremental step; solution increment = difference to solution at beginning of time step
     * Else > incremenetal step; solution incremenet = update to current solution state
     */
    if(steppingTowardsZero){
      solVecForUpdate_ = zeroVec_;
    } else if(steppingTowardsLastTS){
      solVecForUpdate_ = solVec_n_0_;
    } else {
      solVecForUpdate_ = solVec_;
    }

    /*
     * Step 1: setup rhs
     */
    LOG_DBG2(solvestephyst_main) << "--- Substep a - setupRHS ---";
//    std::cout << "ComputeSolutionUpdate_ALGSYS -----> Setup RHS" << std::endl;
    SetupRESorRHS(flag_for_resEvalMatrix, flag_for_rhsVector, solVecForUpdate_, solVec_, solVec_, rhsVec_, false);

    // do not forget to pass rhsVec_ to algsys
    algsys_->InitRHS(rhsVec_);
//    std::cout << "ComputeSolutionUpdate_ALGSYS > Norm of current rhs: " << rhsVec_.NormL2() << std::endl;


    LOG_DBG2(solvestephyst_main) << "--- Norm of rhs: " << rhsVec_.NormL2() << " ---";

    /*
     * Step 2: (re-)assemble matrices according to nonLinMethod
     * Note: During SetupRHSVector() the matrices already get assembled, but
     *        most probably with different flags > reassembly needed
     */
    LOG_DBG2(solvestephyst_main) << "--- Substep b - assembleSystem ---";

//    std::cout << "ComputeSolutionUpdate_ALGSYS -----> Assemble System" << std::endl;
    if(forceJacobian){
      AssembleSystem(FDJACOBIAN,solVec_,forceReuseOfJacobian);
      currentSolutionApproach_ = SOLVE_JACOBIAN;
    } else {
      AssembleSystem(flag_for_systemMatrix, solVec_,forceReuseOfJacobian);
      currentSolutionApproach_ = flagForSolving;
    }

    /*
     * Step 3: setup system, including boundary conditions, solver, etc.
     */
    LOG_DBG2(solvestephyst_main) << "--- Substep c - buildSystem ---";

//    std::cout << "----- Build System FOR SOLVING" << std::endl;
    BuildSystem(false);

    /*
     * In the following cases, we must set DBCs
     * a) Incremental stepping towards last iteration
     *     a1) for iterationCntCurrentTS_ == 1
     *        > for later iterations, the contribution from IDBC is already con-
     *          tained in solVec
     *     a2) after reset (i.e. setting the solution vector back to 0 and start from there
     * b) Incremental stepping towards last timestep
     *     b1) for EACH iteration, as we reset the solution to solVec_n_0_ again
     *     b2) after reset
     *
     */
    bool setIDBC;
    if(steppingTowardsLastTS || steppingTowardsZero || ( iterationCntCurrentTS_ == (iterationOfLatestReset_+1) && couplingIter_ == 0 ) ){
      setIDBC = true;
    } else {
      setIDBC = false;
    }

    bool deltaIDBC;
    if(steppingTowardsZero){
      // if we step towards zero, we do not subtract the old idbc values
      deltaIDBC = false;
    } else {
      deltaIDBC = true;
    }

    /*
     * Step 4: solve system and retrieve the solution (-update)
     */
    //LOG_DBG(solvestephyst_main) << "Solve";
    LOG_DBG2(solvestephyst_main) << "--- Substep d - solveSystem ";
    LOG_DBG3(solvestephyst_main) << "---- setIDBC = " << setIDBC << " / deltaIDBC = " << deltaIDBC;

    algsys_->Solve(setIDBC,deltaIDBC);

    LOG_DBG2(solvestephyst_main) << "--- Substep e - getSolution ";

    // retrieve solution from system;
    // note: solVec_ is not updated automatically (although it holds the reference
    //        to the storage array)
    algsys_->GetSolutionVal(newUpdate, setIDBC, deltaIDBC );

    return newUpdate;
  }

  /*
   * V. System and rhs setup, line search
   * > called during each iteration
   * > inner loop
   * > k = iteration counter
   */
  bool SolveStepHyst::AssembleSystem(matrixFlag requestedAssemblyFlag, SBM_Vector& solutionForMatrixAssembly, bool forceReuseOfJacobian){ //, int iterationCounter){
    LOG_TRACE(solvestephyst_assemble_and_build) << "====================================";
    LOG_TRACE(solvestephyst_assemble_and_build) << "== SolveStepHyst - AssembleSystem ==";
    LOG_TRACE(solvestephyst_assemble_and_build) << "====================================";

    SBM_Vector diff;
    bool reassemblyPerformed = false;
    systemAssemblyRequired_ = false;
    systemBuildRequired_ = false;
    systemLoadRequired_ = false;

    // check if solution has changed since last assembly
    // note: it is not enough to check the buildCounter as the system might require multiple builds during one iteation
    //      (during linesearch!)
    diff = latestUpdateVectors_[requestedAssemblyFlag];
    diff.Add(-1.0,solutionForMatrixAssembly);
    bool solutionChanged = (diff.NormL2() > 1e-16);

    LOG_DBG2(solvestephyst_assemble_and_build) << "---- Difference between previously used assembly-solution and requested assembly-solution: " << diff.NormL2() << " ---";

    LOG_TRACE(solvestephyst_assemble_and_build) << "=== Check if assembly is required ===";
    if(forceReassembly_){
      //    reassembly forced from outside e.g. due to a change in the otherwise constant
      //    permeability (e.g. by switching from B- to H-iteration (FixPoint)
      systemAssemblyRequired_ = true;

//      if(DEBUG_resetForceReassembly_){
//        LOG_DBG(solvestephyst_assemble_and_build) << "--- Force assembly (only this time)  ---";
//        forceReassembly_ = false;
//      } else {
//        LOG_DBG(solvestephyst_assemble_and_build) << "--- Force assembly in general ---";
//      }
    }
    // check if requested system was build at least once during this time step
    if(currentAssemblyFlag_ == NOTSET){
      //    matrix has not been setup at all
      LOG_DBG(solvestephyst_assemble_and_build) << "--- Assemble system for the general first time ---";
      systemAssemblyRequired_ = true;
    } else if(latestBuildCounter_[requestedAssemblyFlag] == -1){
      LOG_TRACE(solvestephyst_assemble_and_build) <<  matrixFlagToString(requestedAssemblyFlag) << " for the first time at iteration " << iterationCntCurrentTS_;
      systemAssemblyRequired_ = true;
    }

    // check individual flags
    if( (requestedAssemblyFlag == HYSTFREE) ||
            (requestedAssemblyFlag == HYSTFREE_GLOBALLY_SCALED_B_VERSION) ||
            (requestedAssemblyFlag == HYSTFREE_LOCALLY_SCALED_B_VERSION) ||
            (requestedAssemblyFlag == HYSTFREE_LOCALLY_SCALED_H_VERSION) ||
            (requestedAssemblyFlag == HYSTFREE_SCALED_H_VERSION) ){
      //  system is free of hysteresis
      //  > assembly required if
      //    a) there are other nonlinearities
      //    b) linear system has not been setup yet (already checked above)
      if(nonHystereticNonLinearitiesPresent_ || materialTensorsHystDepenedent_){
        //        std::cout << "SolveStepHyst::AssembleSystem -- Additonal nonlinearities present "<< std::endl;
        // > reassemble if solution changed
        if(solutionChanged){
          LOG_DBG(solvestephyst_assemble_and_build) << "--- Assemble NONLINEAR-HYSTFREE system due to change in solution ---";
          systemAssemblyRequired_ = true;
        }
      }
    } else if(requestedAssemblyFlag == REUSED){
      //          hyst operator shall return the last computed approximation to the Jacobian; can be
      //            deltaMat or FD approximation (therefore it should be replaced by something more clear in the future)
      //            reassembly only equired if there are additional nonlinearities which were already ruled out here
      //         > no reassembly
      //          std::cout << "SolveStepHyst::AssembleSystem -- REUSED "<< std::endl;
      //          std::cout << " --> no Assemble "<< std::endl;
      LOG_DBG(solvestephyst_assemble_and_build) << "--- Reuse previously assembled system as Jacobian shall be REUSED ---";
      systemAssemblyRequired_ = false;
    } else {
      // FDJACOBIAN, DELTAMATJACOBIAN_LAST_TS, DELTAMATJACOBIAN_LAST_IT
      // > system matrix approximates Jacobian of residual and thus needs to be reassembled if solution changed
      if(solutionChanged){
        LOG_DBG(solvestephyst_assemble_and_build) << "--- Assemble "<<matrixFlagToString(currentAssemblyFlag_,0)<<" system due to change in SOLUTION ---";
        systemAssemblyRequired_ = true;
      }
    }

    if(forceReuseOfJacobian){
      systemAssemblyRequired_ = false;
      systemLoadRequired_ = true;
    }


    // at this point, we have checked, if system needs to be reassembled
    // if yes, we also have to build the system afterwards
    // if no, we might have to load an older version of the system matrix
    if(systemAssemblyRequired_ == false){
      LOG_TRACE(solvestephyst_assemble_and_build) << "=== NO reassembly required ===";
      // we do not reassemble, but we want to have a different matrix (which should already be
      // computed as otherwise systemAssemblyRequired should be true!)
      if(requestedAssemblyFlag != currentBuildFlag_){
        LOG_DBG2(solvestephyst_assemble_and_build) << "---- Current system ("<<matrixFlagToString(currentBuildFlag_,0)<<" does not fit > LOAD BACKUP ---";
        systemLoadRequired_ = true;
      } else {
        LOG_DBG2(solvestephyst_assemble_and_build) << "---- Current system ("<<matrixFlagToString(currentBuildFlag_,0)<<" already fitting > DONE  ---";
      }
    } else {
      LOG_TRACE(solvestephyst_assemble_and_build) << "=== Perform reassembly ===";

      /*
       *  Assemble system according to requestedAssemblyFlag
       *  > retrieve timeLevels for coefFunctionHyst that determine which hysteretic values (if any) shall be
       *      used during the assembly process
       */
      int timeLevelMaterial;
      int FPMaterialTensor;
      int FPRHSCorrectionTensor;
      int timeLevelDeltaMatPol;
      int timeLevelDeltaMatStrain;
      int timeLevelDeltaMatCoupling;
      int inversionApproach;
      Helper_GetHystOperatorParamsFromMatrixFlag(requestedAssemblyFlag, timeLevelMaterial, FPMaterialTensor, FPRHSCorrectionTensor,
              timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling, inversionApproach);

      // TODO: Add option to pass inversionApproach and scalingOfMaterialTensor
      //      > needs to be implemented in CoefFunctionHyst
      //      > also implement usage of these values!
      // evalPurpose 1 > assemble
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelMaterial",timeLevelMaterial);
      //			PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMat",timeLevelDeltaMat);

      if(forceReuseOfJacobian){
        timeLevelDeltaMatPol = 2; // force return of previously stored state
      }

      // if one of the following flags is set to a value != -2, the deltaMatrix will be activated (even if it just returns 0)!
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatPol",timeLevelDeltaMatPol);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatStrain",timeLevelDeltaMatStrain);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelDeltaMatCoupling",timeLevelDeltaMatCoupling);

      if(timeLevelDeltaMatPol == -2){
        // switch of deltaMat computation completely
        PDE_.SetFlagInCoefFncHyst("DeactivateDeltaMat",timeLevelDeltaMatPol);
      }

      // switch special FP tensor on or off
      if((FPMaterialTensor != 0)&&(ReEvalFPMaterialTensors_)){
        PDE_.SetFlagInCoefFncHyst("EvalFPMaterialTensors",FPMaterialTensor); // only once per timestep OR when switching
        // from initial globalFP to localFP!
        ReEvalFPMaterialTensors_ = false;
      }
      PDE_.SetFlagInCoefFncHyst("SetFPMatrixFlag",FPMaterialTensor); // during each assembly!

      // important: allow reevaluation of hyst operator > might have changed
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);

      SBM_Vector solBackup;
      solBackup = solVec_;

      // AssembleMatrices accesses solVec_ for assembly > set solVec_ accordingly
      solVec_ = solutionForMatrixAssembly;
      assemble_->AssembleMatrices(false);

      currentAssemblyFlag_ = requestedAssemblyFlag;
      latestUpdateVectors_[requestedAssemblyFlag] = solutionForMatrixAssembly;

      // if system became reassembled, we also should built it!
      systemAssemblyRequired_ = false;
      systemBuildRequired_ = true;
      reassemblyPerformed = true;

      LOG_TRACE2(solvestephyst_assemble_and_build) << "=== Reassemble finished ===";

      // reset solVec_
      solVec_ = solBackup;
//      timeLevelMaterialCurrent_ = timeLevelMaterial;
//      timeLevelDeltaMatPolCurrent_ = timeLevelDeltaMatPol;
//      timeLevelDeltaMatStrainCurrent_ = timeLevelDeltaMatStrain;
//      timeLevelDeltaMatCouplingCurrent_ = timeLevelDeltaMatCoupling;
    }

    requestedBuildFlag_ = requestedAssemblyFlag;

    return reassemblyPerformed;
  }

  void SolveStepHyst::BuildSystem(bool probibitStorage){
    /*
     * Setup system matrix, boundary conditions etc
     * based upon last call to AssembleSystem
     */
    LOG_TRACE(solvestephyst_assemble_and_build) << "=================================";
    LOG_TRACE(solvestephyst_assemble_and_build) << "== SolveStepHyst - BuildSystem ==";
    LOG_TRACE(solvestephyst_assemble_and_build) << "=================================";

    LOG_DBG(solvestephyst_assemble_and_build) << "--- Current assembly flag: " << matrixFlagToString(currentAssemblyFlag_,0)  << " ---";
    LOG_DBG(solvestephyst_assemble_and_build) << "--- Current build flag: " << matrixFlagToString(currentBuildFlag_,0)  << " ---";
    LOG_DBG(solvestephyst_assemble_and_build) << "--- Requestd build flag: " << matrixFlagToString(requestedBuildFlag_,0) << " ---";

//    std::cout << "--- Current assembly flag: " << matrixFlagToString(currentAssemblyFlag_,0)  << " ---" << std::endl;
//    std::cout << "--- Current build flag: " << matrixFlagToString(currentBuildFlag_,0)  << " ---" << std::endl;
//    std::cout << "--- Requestd build flag: " << matrixFlagToString(requestedBuildFlag_,0) << " ---" << std::endl;
//

    if(systemBuildRequired_ == false){
      // check if system might have to be loaded
      if(systemLoadRequired_){
        LOG_TRACE(solvestephyst_assemble_and_build) << "=== LOAD previously built system ===";
//        std::cout << "=== LOAD previously built system ===" << std::endl;
        FEMatrixType systemStorage;
        if(requestedBuildFlag_ == REUSED){
          switch(currentBuildFlag_){
            case NOTSET:
              EXCEPTION("requestedBuildFlag = NOTSET");
              break;
            case HYSTFREE:
              systemStorage = SYSTEM_HYSTFREE;
              break;
            case HYSTFREE_GLOBALLY_SCALED_B_VERSION:
            case HYSTFREE_LOCALLY_SCALED_B_VERSION:
            case HYSTFREE_LOCALLY_SCALED_H_VERSION:
            case HYSTFREE_SCALED_H_VERSION:
              systemStorage = SYSTEM_FIXPOINT;
              break;
            case FDJACOBIAN:
              systemStorage = SYSTEM_FD_JACOBIAN;
              break;
            case DELTAMATJACOBIAN_LAST_TS:
            case DELTAMATJACOBIAN_LAST_IT:
              systemStorage = SYSTEM_DELTAMAT_JACOBIAN;
              break;
            case REUSED:
              EXCEPTION("requestedBuildFlag = REUSED not allowed if currentBuildFlag is already REUSED");
              break;
            default:
              EXCEPTION("requestedBuildFlag unknown");
              break;
          }
        } else {
          switch(requestedBuildFlag_){
            case NOTSET:
              EXCEPTION("requestedBuildFlag = NOTSET");
              break;
            case HYSTFREE:
              systemStorage = SYSTEM_HYSTFREE;
              break;
            case HYSTFREE_GLOBALLY_SCALED_B_VERSION:
            case HYSTFREE_LOCALLY_SCALED_B_VERSION:
            case HYSTFREE_SCALED_H_VERSION:
            case HYSTFREE_LOCALLY_SCALED_H_VERSION:
              systemStorage = SYSTEM_FIXPOINT;
              break;
            case FDJACOBIAN:
              systemStorage = SYSTEM_FD_JACOBIAN;
              break;
            case DELTAMATJACOBIAN_LAST_TS:
            case DELTAMATJACOBIAN_LAST_IT:
              systemStorage = SYSTEM_DELTAMAT_JACOBIAN;
              break;
            case REUSED:
              EXCEPTION("requestedBuildFlag = REUSED should not occur here");
              break;
            default:
              EXCEPTION("requestedBuildFlag unknown");
              break;
          }
          currentBuildFlag_ = requestedBuildFlag_;
        }
//        std::cout << "LoadSystemBackup " << FEMatrixTypeToString(systemStorage) << std::endl;
        Helper_LoadSystemBackup(false,false,false,true,systemStorage);

      } else {
        LOG_TRACE(solvestephyst_assemble_and_build) << "=== NO LOAD required ===";
//        std::cout << "=== NO LOAD required ===" << std::endl;
      }

      return;
    }

//        if(requestedBuildFlag_ == HYSTFREE){
//          // linear, fixpoint or nonlinear system without hysteresis
//          systemStorage = SYSTEM_HYSTFREE;
//        } else if(requestedBuildFlag_ == FDJACOBIAN){
//          // FD-Jacobian
//          systemStorage = SYSTEM_FD_JACOBIAN;
//        } else if(requestedBuildFlag_ == DELTAMATJACOBIAN_LAST_TS){
//          // DeltaMat-Jacobian
//          systemStorage = SYSTEM_DELTAMAT_JACOBIAN;
//        } else {
//          EXCEPTION("SolveStepHyst::BuildSystem - Cannot backup system; unexpected value of currentAssemblyFlag_ ("<<matrixFlag(currentAssemblyFlag_)<<")");
//        }
//        Helper_LoadSystemBackup(false,false,false,true,systemStorage);
//
//        currentBuildFlag_ = requestedBuildFlag_;
//      }
//      else {
//        LOG_TRACE(solvestephyst_assemble_and_build) << "=== NO LOAD required ===";
//      }
//
//      return;
//    }

    if(currentAssemblyFlag_ == NOTSET){
      WARN("SolveStepHyst::BuildSystem - Cannot build system; please call AssembleSystem first");
      return;
    }

    /*
     * Build system based on latest assembly
     */
    LOG_TRACE(solvestephyst_assemble_and_build) << "=== BUILD system based on latest assembly ===";
//    std::cout << "=== BUILD system based on latest assembly ===" << std::endl;
    matrix_factor_.clear();
    algsys_->InitMatrix(SYSTEM);

    for(fncIt = feFunctions_.begin();fncIt != feFunctions_.end();fncIt++){
      FeFctIdType fctId = fncIt->second->GetFctId();
      fncIt->second->GetTimeScheme()->AddMatFactors(0,matrices,matrix_factor_[fctId]);

      LOG_DBG2(solvestephyst_assemble_and_build) << "---- Current fctId = " << fctId << " / content of matrix_factor_[fctId]: ";

      std::map<FEMatrixType,Double>::iterator innermapit;
      UInt cnt = 0;
      for(innermapit = matrix_factor_[fctId].begin(); innermapit != matrix_factor_[fctId].end(); innermapit++){
        LOG_DBG2(solvestephyst_assemble_and_build) << cnt << ": FeMatrixType = " << FEMatrixTypeToString(innermapit->first) << " / factor = " << innermapit->second;
        cnt++;
      }

      algsys_->ConstructEffectiveMatrix(fctId, matrix_factor_[fctId]);
    }

    PDE_.SetBCs();
    algsys_->BuildInDirichlet();
    algsys_->SetupPrecond();
    algsys_->SetupSolver();

    /*
     * Create backups/store built system for later usage
     */
    LOG_TRACE(solvestephyst_assemble_and_build) << "=== BACKUP system for later usage ===";

    FEMatrixType systemStorage;
    switch(requestedBuildFlag_){
      case NOTSET:
        EXCEPTION("requestedBuildFlag = NOTSET");
        break;
      case HYSTFREE:
        systemStorage = SYSTEM_HYSTFREE;
        break;
      case HYSTFREE_GLOBALLY_SCALED_B_VERSION:
      case HYSTFREE_LOCALLY_SCALED_B_VERSION:
      case HYSTFREE_LOCALLY_SCALED_H_VERSION:
      case HYSTFREE_SCALED_H_VERSION:
        systemStorage = SYSTEM_FIXPOINT;
        break;
      case FDJACOBIAN:
        systemStorage = SYSTEM_FD_JACOBIAN;
        break;
      case DELTAMATJACOBIAN_LAST_TS:
      case DELTAMATJACOBIAN_LAST_IT:
        systemStorage = SYSTEM_DELTAMAT_JACOBIAN;
        break;
      case REUSED:
        EXCEPTION("requestedBuildFlag = REUSED should not occur here");
        break;
      default:
        EXCEPTION("requestedBuildFlag unknown");
        break;
    }
//    std::cout << "BackupSystem to " << FEMatrixTypeToString(systemStorage) << std::endl;
    if(!probibitStorage){
//      std::cout << "BackupSystem to " << FEMatrixTypeToString(systemStorage) << std::endl;
      algsys_->BackupCurrentSystemMatrix(systemStorage);
    } else {
//      std::cout << "Prohibit Backup of System to " << FEMatrixTypeToString(systemStorage) << std::endl;
    }
//
//    if(requestedBuildFlag_ == HYSTFREE){
//      // linear, fixpoint or nonlinear system without hysteresis
//      algsys_->BackupCurrentSystemMatrix(SYSTEM_HYSTFREE);
//      hystFreeSysMatrix_n_k_lastBuildCounter_ = iterationCntCurrentTS_;
//    } else if(requestedBuildFlag_ == FDJACOBIAN){
//      // FD-Jacobian
//      algsys_->BackupCurrentSystemMatrix(SYSTEM_FD_JACOBIAN);
//      FDJacobianSysMatrix_n_k_lastBuildCounter_ = iterationCntCurrentTS_;
//    } else if(requestedBuildFlag_ == DELTAMATJACOBIAN_LAST_TS){
//      // DeltaMat-Jacobian
//      algsys_->BackupCurrentSystemMatrix(SYSTEM_DELTAMAT_JACOBIAN);
//      DeltaMatJacobianSysMatrix_n_k_lastBuildCounter_ = iterationCntCurrentTS_;
//    } else {
//      EXCEPTION("SolveStepHyst::BuildSystem - Cannot backup system; unexpected value of currentAssemblyFlag_ ("<<matrixFlag(currentAssemblyFlag_)<<")");
//    }
    latestBuildCounter_[requestedBuildFlag_] = iterationCntCurrentTS_;
    currentBuildFlag_ = requestedBuildFlag_;
  }

  bool SolveStepHyst::LineSearchHyst_GetInitialSearchInterval(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
          SBM_Vector& currentSol, SBM_Vector& searchDirection,
          Double alpha, Double h, Double t, UInt kMax, Double& a, Double& b, Double& m,
          UInt iterationCounter){
    /*
     * Helper function to determine an intial search interval for linesearch
     * via Forward-Backward method
     *
     * Source: Optimization Theory and Methods - Nonlinear Programming -- Sun, Yuan -- Springer 2006
     *
     * Parameter:
     *  currentSol > solution from previous iteration; starting point for line search
     *  searchDirection > current search/update direction that shall be used to
     *                    reduce the objective function (e.g. the residual norm)
     *                    searchDirection = solution update
     *  alpha > initial stepping distance along searchDirection, that shall be tested
     *           value in [0,infty)
     *  h > initial increment in stepping distance
     *       value > 0
     *  t > scaling factor for h0
     *      value > 1
     *  kMax > maximal number of search steps
     *      value > 1
     *  a,b > lower and upper value of 1d search interval
     *  m > point between a,b ([a,b] is always composed of 2 subintervals [a,m] and [m,b])
     *
     *  iterationCounter > iteration number of nonlinear iteration; used for evaluating the residual
     */
    LOG_TRACE(solvestephyst_linesearch) << "=============================================";
    LOG_TRACE(solvestephyst_linesearch) << "== LineSearchHyst_GetInitialSearchInterval ==";
    LOG_TRACE(solvestephyst_linesearch) << "=============================================";

    Timer* SearchIntervalTimer = new Timer();
    UInt SearchIntervalHystEvalCounter = 0;
    Double startTime = 0.0;
    Double endTime;
    if(DEBUG_measureLSPerformance_){
      startTime = SearchIntervalTimer->GetCPUTime();
      SearchIntervalTimer->Start();
    }

    Double resNorm_k, resNorm_kPlus1;
    Double alpha_k, alpha_kPlus1;
    Double h_k;

    SBM_Vector solBak;
    solBak = solVec_;
    alpha_k = alpha;
    h_k = h;

    SBM_Vector tmpSol;
    tmpSol = currentSol;
    tmpSol.Add(alpha_k,searchDirection);

    solVec_ = tmpSol;
    // set solution vector to tmp solution, then evaluate hyst operators
    // do not reevaluate matrices for inversion > flag = 0
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
    SearchIntervalHystEvalCounter++;

    solVec_ = solBak;

    // SetupRes will write the residual directly to resVec_n_k_
    // It will also reassemble the system! > we have to call this function
    // again after the linesearch to get the curret system set up!
    // TODO:
    // set switch depending on nonlinearity of system; if linear > we can reuse alot
    //    bool fullEval = true;
    //
    //    if(fullEval){
    //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
    //    } else {
    SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, searchDirection, alpha_k);
    //    }
    resNorm_k = resVec_n_k_.NormL2();

//    std::cout << "Initial alpha_k = " << alpha_k << std::endl;
//    std::cout << "Residual for alpha_k = " << resNorm_k << std::endl;
//    
    bool success = false;

    for(UInt k = 0; k < kMax; k++){
      alpha_kPlus1 = alpha_k + h_k;

      SBM_Vector tmpSol;
      tmpSol = currentSol;
      tmpSol.Add(alpha_kPlus1,searchDirection);

      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      SearchIntervalHystEvalCounter++;

      solVec_ = solBak;

      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, searchDirection, alpha_kPlus1);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      resNorm_kPlus1 = resVec_n_k_.NormL2();

//      std::cout << "Initial alpha_k+1 = " << alpha_kPlus1 << std::endl;
//      std::cout << "Residual for alpha_k+1 = " << resNorm_kPlus1 << std::endl;
//      
      if(resNorm_kPlus1 < resNorm_k){
        // Increase stepping distance
        h_k *= t;
        alpha = alpha_k;
        alpha_k = alpha_kPlus1;
        resNorm_k = resNorm_kPlus1;
      } else {
        if(k == 0){
          // Revert stepping direction
          h_k *= -1.0;
          // according to the source, alpha_k should be set to alpha_kPlus1
          // but during the next itartion, we would get alpha_k back and as
          // resNorm_k is not updated (neither is k) according to source, we
          // would obtain the same value again and thus end up in the else again
          // which then would lead to a and b = alpha
          // > instead, alpha_k should be set back to alpha_0 (which the algorithm
          //    actually does by not increasing k) and alpha should be set to
          //    alpha_k+1
          alpha_k = alpha;
          alpha = alpha_kPlus1;
        } else {
          a = std::min(alpha,alpha_kPlus1);
          m = alpha_k;
          b = std::max(alpha,alpha_kPlus1);
          success = true;
          break;
        }
      }
    }

    if(success){
      LOG_DBG(solvestephyst_linesearch) << "--- Initial interval for linesearch [a,m,b] = " << a << "," << m << "," << b;
    } else {
      LOG_DBG(solvestephyst_linesearch) << "--- NO suitable interval for linesearch found after " << kMax << " iterations";
    }

    if(DEBUG_measureLSPerformance_){
      SearchIntervalTimer->Stop();
      endTime = SearchIntervalTimer->GetCPUTime();

      //      std::cout << "=== Statistics for LineSearchHyst_GetInitialSearchInterval ===" << std::endl;
      //      std::cout << "Number of hyst evaluations: " << SearchIntervalHystEvalCounter << std::endl;
      //      std::cout << "Required runtime: " << endTime-startTime << std::endl;

      LOG_DBG(solvestephyst_linesearch) << "=== Statistics for LineSearchHyst_GetInitialSearchInterval ===";
      LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << SearchIntervalHystEvalCounter;
      LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
    }

    delete SearchIntervalTimer;
    return success;

  }

  //TODO:
  //      if FD-Jacobian gets assembled, it might be required to backup and reload current system state (unlikely as we
  //      just need the linear system for residual computation which is actually reload on demand)
  //      change from linesearch_ string to linesearchFlag; allow list of multiple linesearches that shall be usable
  Double SolveStepHyst::LineSearchHyst(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
      SBM_Vector& currentSol, SBM_Vector& solIncrement_n_k_local,UInt iterationCounter){
    /*
     * NEW Linesearch methods
     *
     * EXACT linesearches, i.e. iterate until a certain tolerance is met
     * 1) Golden Section
     * Source: "Optimization Theory and Methods - Nonlinear Programming" - Sun, Yuan -- Springer 2006
     * 2) QuadraticInterpolation
     * Source: "Optimization Theory and Methods - Nonlinear Programming" - Sun, Yuan -- Springer 2006
     *
     *
     * INEXACT linesearches, i.e. iterate until some sufficient decrease in residual was achieved
     * 3) VectorBasedPolynomial
     * Source: "The 'Global' Convergence of Broyden-like methods with a suitable line search" - A. Griewank -- J. Austral. Mathj. 1986
     *
     * 4) ThreePointPolynomial
     * Source: "Iterative Methods for Linear and Nonlinear Equations" - C.T. Kelly -- SIAM 1995
     *
     * 5) Backtracking
     * Sourcess:
     *  "Optimization Theory and Methods - Nonlinear Programming" - Sun, Yuan -- Springer 2006
     *  "Limited memory BFGS method with backtracking for symmetric nonlinear equations" - Yuan, Wei, Lu -- 2011
     *
     * 5.1) Armijo Rule
     * 5.2) Non-monotonic Armijo Rule
     * 5.3) Gradient Free Criterion
     *
     * 6) Trial-And-Error
     *  Source: "Numerical Simulation of Mechatronic Sensors and Actuators" - M. Kaltenbacher -- Springer
     */
    LOG_TRACE(solvestephyst_linesearch) << "========================";
    LOG_TRACE(solvestephyst_linesearch) << "== PERFORM LINESEARCH ==";
    LOG_TRACE(solvestephyst_linesearch) << "========================";

    SBM_Vector solBak;
    SBM_Vector tmpSol;
    SBM_Vector resInitial;

    solBak = solVec_;

    Timer* LinesearchTimer = new Timer();
    Double startTime = 0.0;
    Double endTime;
    Double etaReturn = 1.0;
//    Double tol = 1e-3;
    UInt k = 0;

    bool validLSperformed = false;
    bool validEtaFound = false;

    // store name of ls as key; found alpha, corresponding residual, num evals, evaltime inside the vector
    std::stack< LSResults > testedLSMethods;

    /*
     * Evaluate initial residual, i.e. for solIncrement_n_k_local = 0
     * > should be ths same as the currently stored value of resVec_n_k_
     */
    LOG_DBG(solvestephyst_linesearch) << "--- Determine initial residual ---";

    tmpSol = currentSol;
    //      tmpSol.Add(0,solIncrement_n_k_local);
    solVec_ = tmpSol;
    // set solution vector to tmp solution, then evaluate hyst operators
    // do not reevaluate matrices for inversion > flag = 0
    PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
    solVec_ = solBak;

    // SetupRes will write the residual directly to resVec_n_k_
    // It will also reassemble the system! > we have to call this function
    // again after the linesearch to get the current system set up!
    resInitial = zeroVec_;
//    Double normOfSystemFix =
    SetupRESorRHS(flag_for_matrix_assembly, flag_for_vector_assembly,
      tmpSol, tmpSol, tmpSol, resInitial, false);

    Double referenceValue = resVec_n_0_.NormL2();

    LOG_DBG2(solvestephyst_linesearch) << "---- norm of initial residual = " << resVec_n_k_.NormL2() << " ----";
//    std::cout << "Start of LS - Initial residual (everything reassembled): " << resVec_n_k_.NormL2() << std::endl;
//    std::cout << "Start of LS - reference residual resVec_n_0_ (from start of timestep): " << resVec_n_0_.NormL2() << std::endl;
//    
    if(checkDirection_ >= 1){
      /*
       * According to Kelly 1995 (p.138), check if search direction satisfies the
       * inexact Newton property:
       * ||Jac*solIncrement_n_k_local + res|| < eta*||res||
       *
       * > what a stupid idea in my case:
       *      we actually solved Jac*solIncrement_n_k_local = -res
       *      i.e. the condition has to be satisfied
       *
       * exceptions:
       *      SystemMatrix != Jacobian!
       *      > if Broyden (BGM) or deltaMat was used
       *      > but in this case, how to check then? Jacobian would have to be
       *        aseembled! > would be possible, but if we consider to assemble
       *        the Jacobian, why don't we just solve the system with it?
       *      > furthermore: if we assemble using the Jacobian, then we only have
       *        a FiniteDifference approximation, so the check might be in vain
       *
       */
      LOG_TRACE(solvestephyst_linesearch) << "=== Check for inexact-Newton update direction ===";

      if(DEBUG_measureLSPerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      Double eta = 0.9;
      SBM_Vector JacD;
      JacD = resVec_n_k_;
      JacD.Init();

      solVec_ = currentSol;

      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);

      if(currentBuildFlag_ != FDJACOBIAN){
        LOG_DBG(solvestephyst_linesearch) << "--- FD-Jacobian needed for check; assembly of system required ---";

        // build up system using FD-Jacobian
        AssembleSystem(FDJACOBIAN,currentSol,false);//, -1);
        // do not store the system to backup storage
        // this is important for e.g. chord method; here we do not want to overwrite the FD_JACOBIAN storage
        BuildSystem(true);
      }
      algsys_->ComputeSysMatTimesVector(SYSTEM,solIncrement_n_k_local,JacD,false);

      solVec_ = solBak;

      //remember: we have switched sign between residual and jacobian -> subtract
      JacD.Add(-1.0,resInitial);

      if(DEBUG_measureLSPerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for checking inexact-Newton update direction ===" << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for checking inexact-Newton update direction ===";
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      }

      if(JacD.NormL2() < eta*resInitial.NormL2()){
        LOG_TRACE(solvestephyst_linesearch) << "=== Solution increment satisfies condition for inexact-Newton step; PROCEED ===";
      } else {
//        std::cout << "=== Solution increment does not satisfy condition for inexact-Newton step; ";
        LOG_TRACE(solvestephyst_linesearch) << "=== Solution increment does not satisfy condition for inexact-Newton step; ";
        if(checkDirection_ >= 2){
//          std::cout << "DISMISS UPDATE ===" << std::endl;
          LOG_TRACE(solvestephyst_linesearch) << "DISMISS UPDATE ===";
          return 0.0;
        } else {
//          std::cout << "PROCEED ANYWAY ===" << std::endl;
          LOG_TRACE(solvestephyst_linesearch) << "PROCEED ANYWAY ===";
        }
      }
    }

    /*
     * Preparation step for GoldenSection and QuadraticInterpolation
     */
    Double aBAK = 0.0;
    Double bBAK = 0.0;
    Double mBAK = 0.0;

     bool initialIntervalFound = false;
    if( (specifiedLinesearches_.find(LS_GOLDENSECTION) != specifiedLinesearches_.end()) ||
        (specifiedLinesearches_.find(LS_QUADINTERPOL) != specifiedLinesearches_.end()) ){

      LOG_TRACE(solvestephyst_linesearch) << "=== Determine initial search interval for EXACT_GoldenSection and/or EXACT_QuadraticInterpolation ===";
      // LineSearchHyst_GetInitialSearchInterval has built-in timer > no timer needed here
      UInt kMaxInitialInterval = 20;

      initialIntervalFound = LineSearchHyst_GetInitialSearchInterval(flag_for_matrix_assembly, flag_for_vector_assembly,
              currentSol, solIncrement_n_k_local, 0.0, 0.1, 1.5, kMaxInitialInterval, aBAK, bBAK, mBAK, iterationCounter);
    }

    if( (specifiedLinesearches_.find(LS_GOLDENSECTION) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "EXACT_GoldenSection") ){
      /*
       * Source: Optimization Theory and Methods - Nonlinear Programming -- Sun, Yuan -- Springer 2006
       * (Algorithm 2.3.1 on page 86)
       *
       * Base facts:
       *  - linear convergence (rate of 0.618)
       *  - objective function = residual norm
       *  - requires no derivatives
       *
       * Base idea:
       *  - find optimal step length by reducing initial search interval by
       *    golden ration (0.618)
       *
       */

      /*
       * Load read-in parameter
       */
      LinesearchParameter GoldenSectionParameter = specifiedLinesearchesMap_[LS_GOLDENSECTION];

//      std::cout << "GoldenSectionLS - initial interval: " << std::endl;
//      std::cout << "a = " << aBAK << std::endl;
//      std::cout << "m = " << mBAK << std::endl;
//      std::cout << "b = " << bBAK << std::endl;
//      
      /*
       * Check initial interval
       */
      bool initialResidualFail = false;
      if(!initialIntervalFound){
        LOG_TRACE(solvestephyst_linesearch) << "=== No initial interval could be found ===";
        initialResidualFail = true;
      }

      if( (aBAK < 0 && bBAK < 0)&&(!GoldenSectionParameter.allowNegativeEta_)){
        LOG_TRACE(solvestephyst_linesearch) << "=== Initial interval is completely negative ===";
        initialResidualFail = true;
      }

      if( (initialResidualFail)&&(specifiedLinesearches_.size() <= 1)){
        LOG_TRACE(solvestephyst_linesearch) << "=== No valid starting interval found and no backup linesearch defined; retturn full step ===";
        return 1.0;
      }

      // Golden section has no minEta!
//      // upper boundary of interval is smaller than minimal allowed eta > use minimal eta directly
//      if(bBAK<GoldenSectionParameter.minEta_){
//        LOG_TRACE(solvestephyst_linesearch) << "=== Upper bound of starting interval is already below etaMin " << GoldenSectionParameter.minEta_ << " ===";
//        bBAK = GoldenSectionParameter.minEta_;
//        aBAK = GoldenSectionParameter.minEta_;
//        // search interval restricted to length 0 > exact linesearches will stop after first iteration
//      }

//      std::cout << "Initial Interval [a,b] = " << aBAK << "," << bBAK << std::endl;

      LOG_TRACE(solvestephyst_linesearch) << "=== Perform GoldenSection linesarch ===";
      if(GoldenSectionParameter.measurePerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      UInt EXACT_GoldenSectionEvalCounter = 0;

      Double alphaEXACT_GoldenSection = 1.0;
      Double resOptEXACT_GoldenSection = 1.0e16;
      Double lambda,mu;
      Double resLambda,resMu;
      bool successEXACT_GoldenSection = false;
      Double a,b;
      a = aBAK;
      b = bBAK;

      /*
       * Compute two initial observer at lambda and mu
       */
      lambda = a + 0.382*(b-a);
      mu = a + 0.618*(b-a);

      /*
       * for debugging and comparison: compute residual for stepping 1.0 (full step)
       * 
       */
      bool checkFullStep = !true;
      // Test changes result! why?! > due to EXACT_GoldenSectionEvalCounter++! it was just one iteration short!
      // Thats not all! After checking the trial and error values below resiudal increases and increases
      // > are solution values not reset correctly?
      if(checkFullStep){
        std::cout << "Test full step first" << std::endl;
//        std::cout << "Norms of tmpSol, currentSol, solVec_, solBak at start: " << std::endl;
//        std::cout << "tmpSol: " << tmpSol.NormL2() << std::endl;
//        std::cout << "currentSol: " << currentSol.NormL2() << std::endl;
//        std::cout << "solVec_: " << solVec_.NormL2() << std::endl;
//        std::cout << "solBak: " << solBak.NormL2() << std::endl;
        
        tmpSol = currentSol;
        tmpSol.Add(1.0,solIncrement_n_k_local);
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
//        EXACT_GoldenSectionEvalCounter++;
        solVec_ = solBak;

        // SetupRes will write the residual directly to resVec_n_k_
        // It will also reassemble the system! > we have to call this function
        // again after the linesearch to get the curret system set up!
        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, 1.0);
        //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        Double resFullstep = resVec_n_k_.NormL2();
        std::cout << "Residual for full-step: " << resFullstep << std::endl;
//        std::cout << "Norms of tmpSol, currentSol, solVec_, solBak at end: " << std::endl;
//        std::cout << "tmpSol: " << tmpSol.NormL2() << std::endl;
//        std::cout << "currentSol: " << currentSol.NormL2() << std::endl;
//        std::cout << "solVec_: " << solVec_.NormL2() << std::endl;
//        std::cout << "solBak: " << solBak.NormL2() << std::endl;
      }
      bool checkSomeTrialSteps = !true;
      // Test changes result! why?! > due to EXACT_GoldenSectionEvalCounter++! it was just one iteration short!
      // issue detected for FP-H-Version: 
      // although system, solution etc were reset correctly, the storages E_H_ and P_H_ in coefFunction hyst were
      // overwritten each iteration; this is usually wanted to check if evaluation is needed at all but in case of
      // FP-H-Version it is not! here E_H_ and P_H_ are used for computing the new state; if E_H_ and P_H_ are set
      // to temporal values from within the Linesearch, we muddle with the next iterates; this can lead to searious
      // increases in the computed residuals!
      if(checkSomeTrialSteps){
        Vector<Double> testAlphas = Vector<Double>(10);
        testAlphas[0] = 1.0;
        testAlphas[1] = 0.5;
        testAlphas[2] = 0.1;
        testAlphas[3] = 1.0;
        testAlphas[4] = 0.5;
        testAlphas[5] = 0.1;
        testAlphas[6] = 1.0;
        testAlphas[7] = 0.5;
        testAlphas[8] = 0.1;
        testAlphas[9] = 0.000;
        for(UInt i = 0; i<10; i++){
          Double testAlpha = testAlphas[i];
          std::cout << "Current ALPHA = " << testAlpha << std::endl;
//          std::cout << "Norms of tmpSol, currentSol, solVec_, solBak at start: " << std::endl;
//          std::cout << "tmpSol: " << tmpSol.NormL2() << std::endl;
//          std::cout << "currentSol: " << currentSol.NormL2() << std::endl;
//          std::cout << "solVec_: " << solVec_.NormL2() << std::endl;
//          std::cout << "solBak: " << solBak.NormL2() << std::endl;
          tmpSol = currentSol;
          tmpSol.Add(testAlpha,solIncrement_n_k_local);
          solVec_ = tmpSol;
          // set solution vector to tmp solution, then evaluate hyst operators
          // do not reevaluate matrices for inversion > flag = 0
          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
  //        EXACT_GoldenSectionEvalCounter++;
          solVec_ = solBak;

          // SetupRes will write the residual directly to resVec_n_k_
          // It will also reassemble the system! > we have to call this function
          // again after the linesearch to get the curret system set up!
          SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, testAlpha);
          //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
          std::cout << "Residual: " << resVec_n_k_.NormL2() << std::endl;
//          std::cout << "Norms of tmpSol, currentSol, solVec_, solBak at end: " << std::endl;
//          std::cout << "tmpSol: " << tmpSol.NormL2() << std::endl;
//          std::cout << "currentSol: " << currentSol.NormL2() << std::endl;
//          std::cout << "solVec_: " << solVec_.NormL2() << std::endl;
//          std::cout << "solBak: " << solBak.NormL2() << std::endl;
        }
        
      }

      tmpSol = currentSol;
      tmpSol.Add(lambda,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      EXACT_GoldenSectionEvalCounter++;
      solVec_ = solBak;

//      std::cout << "Initial lambda (lambda = a + 0.382*(b-a)): " << lambda << std::endl;
//      resLambda = resVec_n_k_.NormL2();
//      std::cout << "resVec_n_k_.NormL2() pre Call To SetupResVector_Linesearch: " << resLambda << std::endl;
      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, lambda);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      resLambda = resVec_n_k_.NormL2();

//      std::cout << "Residual for lambda: " << resLambda << std::endl;
      
      tmpSol = currentSol;
      tmpSol.Add(mu,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      EXACT_GoldenSectionEvalCounter++;
      solVec_ = solBak;
//      std::cout << "Initial mu (mu = a + 0.618*(b-a)): " << mu << std::endl;
//      resMu = resVec_n_k_.NormL2();
//      std::cout << "resVec_n_k_.NormL2() pre Call To SetupResVector_Linesearch: " << resMu << std::endl;
      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, mu);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      resMu = resVec_n_k_.NormL2();
      
//      std::cout << "Residual for mu: " << resMu << std::endl;
      
      LOG_DBG2(solvestephyst_linesearch) << "Check optimality condition ( ((mu - a)<=tol)||((b - lambda)<=tol) ) for GoldenSection: ";
      for(k = 0; k < GoldenSectionParameter.maxIter_; k++){
        if(resLambda > resMu){
          LOG_DBG2(solvestephyst_linesearch) << "lambda = " << lambda << "; b = " << b << "; GoldenSectionParameter.stoppingTol_ = " << GoldenSectionParameter.stoppingTol_;
          if((b - lambda)<=GoldenSectionParameter.stoppingTol_){
            alphaEXACT_GoldenSection = mu;
            resOptEXACT_GoldenSection = resMu;
            successEXACT_GoldenSection = true;
            break;
          } else {
            a = lambda;
            // b = b
            lambda = mu;
            resLambda = resMu;
            mu = a + 0.618*(b - a);

            tmpSol = currentSol;
            tmpSol.Add(mu,solIncrement_n_k_local);
            solVec_ = tmpSol;
            // set solution vector to tmp solution, then evaluate hyst operators
            // do not reevaluate matrices for inversion > flag = 0
            PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
            EXACT_GoldenSectionEvalCounter++;
            solVec_ = solBak;

            // SetupRes will write the residual directly to resVec_n_k_
            // It will also reassemble the system! > we have to call this function
            // again after the linesearch to get the curret system set up!
//            std::cout << "k = " << k << ": mu (mu = a + 0.618*(b-a)): " << mu << std::endl;
//            resMu = resVec_n_k_.NormL2();
//            std::cout << "resVec_n_k_.NormL2() pre Call To SetupResVector_Linesearch: " << resMu << std::endl;
            SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, mu);
            //            SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
            resMu = resVec_n_k_.NormL2();
//            std::cout << "k = " << k << ": mu (mu = a + 0.618*(b-a)): " << mu << std::endl;
//            std::cout << "Residual for mu: " << resMu << std::endl;
          }
        } else {
          LOG_DBG2(solvestephyst_linesearch) << "mu = " << mu << "; a = " << a << "; GoldenSectionParameter.stoppingTol_ = " << GoldenSectionParameter.stoppingTol_;
          if((mu - a)<=GoldenSectionParameter.stoppingTol_){
            alphaEXACT_GoldenSection = lambda;
            resOptEXACT_GoldenSection = resLambda;
            successEXACT_GoldenSection = true;
            break;
          } else {
            // a = a
            b = mu;
            mu = lambda;
            resMu = resLambda;
            lambda = a + 0.382*(b-a);

            tmpSol = currentSol;
            tmpSol.Add(lambda,solIncrement_n_k_local);
            solVec_ = tmpSol;
            // set solution vector to tmp solution, then evaluate hyst operators
            // do not reevaluate matrices for inversion > flag = 0
            PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
            EXACT_GoldenSectionEvalCounter++;
            solVec_ = solBak;
//            std::cout << "k = " << k << ": lambda (lambda = a + 0.382*(b-a)): " << lambda << std::endl;
//            resLambda = resVec_n_k_.NormL2();
//            std::cout << "resVec_n_k_.NormL2() pre Call To SetupResVector_Linesearch: " << resLambda << std::endl;
            // SetupRes will write the residual directly to resVec_n_k_
            // It will also reassemble the system! > we have to call this function
            // again after the linesearch to get the curret system set up!
            SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, lambda);
            //            SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
            resLambda = resVec_n_k_.NormL2();

//            std::cout << "Residual for lambda: " << resLambda << std::endl;
          }
        }
      }

      if(GoldenSectionParameter.measurePerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for GoldenSection-Linesearch ===" << std::endl;
        //        std::cout << "Number of iteraion steps: " << k << std::endl;
        //        std::cout << "Number of hyst evaluations: " << EXACT_GoldenSectionEvalCounter << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for GoldenSection-Linesearch ===";
        LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
        LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << EXACT_GoldenSectionEvalCounter;
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      } else {
        endTime = 0.0;
      }

      // flag indicating if at least one linesearch was performed (this is the case at this point)
      // this does not mean, that linesearch was successful!
      validLSperformed = true;

	  /*
       * check if linesearch was successful, too
       */
      if((alphaEXACT_GoldenSection < 0)&&(!GoldenSectionParameter.allowNegativeEta_)){
//        std::cout << "Negative eta not allowed" << std::endl;
        validEtaFound = false;
        successEXACT_GoldenSection = false;
      }

      if(successEXACT_GoldenSection){
//        std::cout << "GoldenSection: validEtaFound; " << alphaEXACT_GoldenSection << std::endl;
        validEtaFound = true;
      }
//      else {
//        std::cout << "Remaining intervals: " << std::endl;
//        std::cout << "mu-a: " << mu-a << std::endl;
//        std::cout << "b-lambda: " << b-lambda << std::endl;
//        std::cout << "mu: " << mu << std::endl;
//        std::cout << "Remaining residual for mu: " << resMu << std::endl;
//        std::cout << "lambda: " << lambda << std::endl;
//        std::cout << "Remaining residual for lambda: " << resLambda << std::endl;
//        std::cout << "--- stop ---" << std::endl;
//        exit(0);
//      }

      LSResults curMethod;
      curMethod.name = "EXACT_GoldenSection";
      curMethod.lsFlag = LS_GOLDENSECTION;
      curMethod.success = successEXACT_GoldenSection;
      curMethod.alphaFound = alphaEXACT_GoldenSection;
      curMethod.resultingRes = resOptEXACT_GoldenSection/referenceValue;
      curMethod.totalEvalTime = endTime-startTime;
      curMethod.numResidualEvals = EXACT_GoldenSectionEvalCounter;

      testedLSMethods.push(curMethod);

      etaReturn = alphaEXACT_GoldenSection;
    }

    if( (specifiedLinesearches_.find(LS_QUADINTERPOL) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "EXACT_QuadraticInterpolation") ){
      /*
       * Source: Optimization Theory and Methods - Nonlinear Programming -- Sun, Yuan -- Springer 2006
       * (Fig. 2.4.1 on page 97)
       *
       * Base facts:
       *  - linear convergence (rate of 1.312)
       *  - objective function = residual norm
       *  - requires no derivatives but more function evaluations than EXACT_GoldenSection
       *
       * Base idea:
       *  - approximate residual/objective function by quadratic polynomial
       *    (can be done from three points or two points and one derivative)
       *
       * Notes: Fig 2.4.1 checks k = 2 at the beginning; it is not said was k actually denotes
       *        but it looks like it is the number of initial intervals
       *      > from our helper function we get 1 interval
       *      > set alpha2 to a and alpha3 to b
       *
       * Note: let the interval method return the inner point m, too
       *      > we start with case k=2 from sketch
       *
       */
       /*
       * Load read-in parameter
       */
      LinesearchParameter QuadInterpolParameter = specifiedLinesearchesMap_[LS_QUADINTERPOL];

      /*
       * Check initial interval
       */
      bool initialResidualFail = false;
      if(!initialIntervalFound){
        LOG_TRACE(solvestephyst_linesearch) << "=== No initial interval could be found ===";
        initialResidualFail = true;
      }

      if( (aBAK < 0 && bBAK < 0)&&(!QuadInterpolParameter.allowNegativeEta_)){
        LOG_TRACE(solvestephyst_linesearch) << "=== Initial interval is completely negative ===";
        initialResidualFail = true;
      }

      if( (initialResidualFail)&&(specifiedLinesearches_.size() <= 1)){
        LOG_TRACE(solvestephyst_linesearch) << "=== No valid starting interval found and no backup linesearch defined; retturn full step ===";
        return 1.0;
      }

      // QuadInterpol has no minEta
//      // upper boundary of interval is smaller than minimal allowed eta > use minimal eta directly
//      if(bBAK<QuadInterpolParameter.minEta_){
//        LOG_TRACE(solvestephyst_linesearch) << "=== Upper bound of starting interval is already below etaMin " << QuadInterpolParameter.minEta_ << " ===";
//        bBAK = QuadInterpolParameter.minEta_;
//        aBAK = QuadInterpolParameter.minEta_;
//        // search interval restricted to length 0 > exact linesearches will stop after first iteration
//      }

      LOG_TRACE(solvestephyst_linesearch) << "=== Perform QuadraticInterpolation linesarch ===";
      if(QuadInterpolParameter.measurePerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      UInt EXACT_QuadraticInterpolationEvalCounter = 0;

      Double alphaQuadInterpol = 1.0;
      Double resOptQuadInterpol = 1.0e16;
      Double alpha, alpha1, alpha2, alpha3;
      Double resAlpha, resAlpha1, resAlpha2, resAlpha3;
      resAlpha = -1.0;
      alpha = -1e16;

      UInt caseK = 1;
      if(caseK == 1){
        alpha1 = -1e16;
        alpha2 = aBAK;
        alpha3 = bBAK;
        resAlpha1 = 0;
      } else {
        alpha1 = aBAK;
        alpha2 = mBAK;
        alpha3 = bBAK;

        // evaluate with alpha1
        tmpSol = currentSol;
        tmpSol.Add(alpha1,solIncrement_n_k_local);
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        EXACT_QuadraticInterpolationEvalCounter++;
        solVec_ = solBak;

        // SetupRes will write the residual directly to resVec_n_k_
        // It will also reassemble the system! > we have to call this function
        // again after the linesearch to get the curret system set up!
        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alpha1);
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        resAlpha1 = resVec_n_k_.NormL2();
      }

      bool successQuadInterpolation = false;

      tmpSol = currentSol;
      tmpSol.Add(alpha2,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      EXACT_QuadraticInterpolationEvalCounter++;
      solVec_ = solBak;

      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alpha2);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      resAlpha2 = resVec_n_k_.NormL2();

      tmpSol = currentSol;
      tmpSol.Add(alpha3,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      EXACT_QuadraticInterpolationEvalCounter++;
      solVec_ = solBak;

      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alpha3);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      resAlpha3 = resVec_n_k_.NormL2();

      Double tmp1,tmp2, A;

      LOG_DBG2(solvestephyst_linesearch) << "Check optimality condition ((abs(alpha - alpha2) < tol) for QuadraticInterpolation: ";
      for(k = 0; k < QuadInterpolParameter.maxIter_; k++){

        if( (caseK == 2)&&(k == 0)){
          // jump directly to evaluation of A
        } else {
          if( (caseK == 1) ){
            // only one interval known; use midpoint
            alpha = (alpha2+alpha3)/2.0;
            // use midpoint alpha
            tmpSol = currentSol;
            tmpSol.Add(alpha,solIncrement_n_k_local);
            solVec_ = tmpSol;
            // set solution vector to tmp solution, then evaluate hyst operators
            // do not reevaluate matrices for inversion > flag = 0
            PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
            EXACT_QuadraticInterpolationEvalCounter++;
            solVec_ = solBak;

            // SetupRes will write the residual directly to resVec_n_k_
            // It will also reassemble the system! > we have to call this function
            // again after the linesearch to get the curret system set up!
            SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alpha);
            //            SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
            resAlpha = resVec_n_k_.NormL2();
          }

          if(resAlpha <= -1.0){
            EXCEPTION("resAlpha = -1.0, i.e it was not set yet; this must not happen at this point");
          }
          if(alpha <= -1e16){
            EXCEPTION("alpha = -1e16, i.e it was not set yet; this must not happen at this point");
          }

          LOG_DBG2(solvestephyst_linesearch) << "alpha = " << alpha << "; alpha2 = " << alpha2 << "; QuadInterpolParameter.stoppingTol_ = " << QuadInterpolParameter.stoppingTol_;
          if( abs(alpha - alpha2) < QuadInterpolParameter.stoppingTol_){
            resOptQuadInterpol = resAlpha;
            alphaQuadInterpol = alpha;
            successQuadInterpolation = true;
            break;
          }

          tmp1 = alpha - alpha2;
          tmp2 = alpha - alpha3;

          if(tmp1*tmp2 < 0){
            // case1
            // search in interval alpha2 to alpha3
            if(resAlpha < resAlpha2){
              // case1.1
              // new three point interval
              // [alpha2, alpha, alpha3]
              alpha1 = alpha2;
              resAlpha1 = resAlpha2;

              alpha2 = alpha;
              resAlpha2 = resAlpha;
            } else {
              // case1.2
              // new three point interval
              // [alpha1, alpha2, alpha]
              alpha3 = alpha;
              resAlpha3 = resAlpha;
            }
          } else {
            // case2
            // search in interval alpha1 to alpha2
            if(resAlpha < resAlpha2){
              // case2.1
              // new three point interval
              // [alpha1, alpha, alpha2]
              alpha3 = alpha2;
              resAlpha3 = resAlpha2;

              alpha2 = alpha;
              resAlpha2 = resAlpha;
            } else {
              // case2.2
              // new three point interval
              // [alpha, alpha2, alpha3]
              alpha1 = alpha;
              resAlpha1 = resAlpha;
            }
          }
        }

        // NOTE: for computation of A we need alpha1 and resAlpha1
        // during the first iteration however we start with only alpha2 and alpha3
        // due to alpha = (alpha2+alpha3)/2 at the beginning, tmp1*tmp2 above
        // will be negative > case1
        // Case 1.1 is ok as alpha1 will be set here, but case 1.2 is problematic
        // as up to this point neither alpha1 nor res(alpha1) are known
        // > case 1.2 can occur, if residual at midpoint is larger than residual
        //    at alpha2 which is not too unlikely
        // > new attempt: replace alpha3 with alpha (which should at least work better than
        //        one of the corner points and start again with caseK=1; else set caseK to 2
        //  > works
        if(alpha1 <= -1e16){
          // alpha1 <= -1e16 only for caseK = 1 > for this case, alpha and resAlpha are known
          //          WARN("Case 1.2 during first iteration!");
          alpha3 = alpha;
          resAlpha3 = resAlpha;
          caseK = 1;
        } else {
          caseK = 2;
          // compute A = denominator for roots of quadratic approximation
          A = 2*(resAlpha1*(alpha2-alpha3) + resAlpha2*(alpha3-alpha1) + resAlpha3*(alpha1-alpha2));

          if(A == 0){
            //            std::cout << "A == 0" << std::endl;
            // Failback case > during next iteration alpha-alpha2 < tol and thus we stop
            alpha = alpha2;
            resAlpha = resAlpha2; //source uses resAlpha2 here (which is not really needed at all); we set to resAlpha2 so that it fits to alphas
          } else {

            alpha = resAlpha1*(alpha2*alpha2-alpha3*alpha3) +  resAlpha2*(alpha3*alpha3-alpha1*alpha1) + resAlpha3*(alpha1*alpha1-alpha2*alpha2);
            alpha /= A;

            tmp1 = alpha-alpha1;
            tmp2 = alpha-alpha3;

            if(tmp1*tmp2 > 0){
              //              std::cout << "tmp1*tmp2 > 0" << std::endl;
              // alpha would be outside of interval [alpha1,alpha3] > fail
              alpha = alpha2;
              resAlpha = resAlpha2; //source uses resAlpha2 here (which is not really needed at all); we set to resAlpha2 so that it fits to alphas
            } else {
              // alpha inside interval [alpha1,alpha3] > ok!
              tmpSol = currentSol;
              tmpSol.Add(alpha,solIncrement_n_k_local);
              solVec_ = tmpSol;
              // set solution vector to tmp solution, then evaluate hyst operators
              // do not reevaluate matrices for inversion > flag = 0
              PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
              EXACT_QuadraticInterpolationEvalCounter++;
              solVec_ = solBak;

              // SetupRes will write the residual directly to resVec_n_k_
              // It will also reassemble the system! > we have to call this function
              // again after the linesearch to get the curret system set up!
              SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alpha);
              //              SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
              resAlpha = resVec_n_k_.NormL2();
            }
          }
        }
      }

      if(QuadInterpolParameter.measurePerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for QuadraticInterpolation-Linesearch ===" << std::endl;
        //        std::cout << "Number of iteraion steps: " << k << std::endl;
        //        std::cout << "Number of hyst evaluations: " << EXACT_QuadraticInterpolationEvalCounter << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for QuadraticInterpolation-Linesearch ===";
        LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
        LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << EXACT_QuadraticInterpolationEvalCounter;
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      } else {
        endTime = 0.0;
      }

      // flag indicating if at least one linesearch was performed (this is the case at this point)
      // this does not mean, that linesearch was successful!
      validLSperformed = true;

      /*
       * check if linesearch was successful, too
       */
      if((alphaQuadInterpol < 0)&&(!QuadInterpolParameter.allowNegativeEta_)){
        validEtaFound = false;
        successQuadInterpolation = false;
      }

      if(successQuadInterpolation){
        validEtaFound = true;
      }

      LSResults curMethod;
      curMethod.name = "EXACT_QuadraticInterpolation";
      curMethod.lsFlag = LS_QUADINTERPOL;
      curMethod.success = successQuadInterpolation;
      curMethod.alphaFound = alphaQuadInterpol;
      curMethod.resultingRes = resOptQuadInterpol/referenceValue;
      curMethod.totalEvalTime = endTime-startTime;
      curMethod.numResidualEvals = EXACT_QuadraticInterpolationEvalCounter;

      testedLSMethods.push(curMethod);

      etaReturn = alphaQuadInterpol;
    }

    if( (specifiedLinesearches_.find(LS_VECBASEDPOLY) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "INEXACT_VectorBasedPolynomial") ) {
      /*
       * Source: "The 'Global' Convergence of Broyden-like methods with a suitable line search"
       * - A. Griewank - J. Austral. Mathj. 1986
       *
       * INEXACT LINESEARCH, i.e. take largest steplength that satisfies sufficient decrease
       *
       * Three major points:
       * a) allow search direcation to be in the NULL space of the current (Broyden-)approximation
       *    to the Jacobian (unfortunately it is not stated, when to do this)
       * b) linesearch is based on a polynomial approximation of the VECTOR residual (not its norm)
       *    > here we implement this linesearch (see paper p. 80-81)
       * c) resulting step length can be negative (or even zero!)
       *
       * Notation:
       * alphaK = step length
       * gK = residual for alphaK = 0
       * gK_c = residual for some alphaC != 0 (take 1 here)
       * gK_alphaK = gK + (gK_c - gK)*alphaK/alphaC
       * qK_alphaK = gK^T(gK - gK_alphaK)/|gK - gK_alphaK|^2
       *
       * > qK_alphaK acts as criterion for selecting the next alphaK
       *
       * alphaK_new = alphaK*h(qK_alphaK)
       *
       * with h(x) =
       *  1             if x >= 0.5 + epsilon
       *  x             if epsilon <= x < 0.5 + epsilon
       *  epsilon       if 0 <= x < epsilon
       * -epsilon       if -epsilon <= x < 0
       * x              if -0.5-epsilon <= x < -epsilon
       * -0.5-epsilon   if x < -0.5-epsilon
       *
       * Stopping criterion (2.7):
       *  qK_alphaK >= 0.5 + epsilon
       *
       * Notes:
       * > due to usage of linear polynomial, we just need to evaluate two residuals here
       *  (gK and gK_c)
       * > function h(x) not given explicitely as formula but in form of Figure 2
       * > epsilon in (0,1/6)
       *
       */
      /*
       * Load read-in parameter
       */
      LinesearchParameter VecPolyParameter = specifiedLinesearchesMap_[LS_VECBASEDPOLY];

      LOG_TRACE(solvestephyst_linesearch) << "=== Perform vector-based polynomial linesarch ===";
      if(VecPolyParameter.measurePerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      UInt INEXACT_VectorBasedPolynomialEvalCounter = 0;

      Double alphaC = 1.0; // second value required for interpolation // keep fix
      Double alphaK = 1.0; // current stepping length > start with etaStart
      Double epsilon = VecPolyParameter.epsilon_; // slightly smaller than 1/6

      /*
       * step 1: compute residuals
       */
      SBM_Vector gK;
      SBM_Vector gK_alphaK;
      SBM_Vector gK_c;
      Double qK_alphaK;

      /*
       * Eval initial residual (i.e. for no stepping)
       * > already done at beginning of linesearch
       * > reload should be possible
       */
      //      /*
      //       * TODO: BEGIN RELOAD
      //       */
      //      tmpSol = currentSol;
      //      solVec_ = tmpSol;
      //      // set solution vector to tmp solution, then evaluate hyst operators
      //      // do not reevaluate matrices for inversion > flag = 0
      //      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      //      INEXACT_VectorBasedPolynomialEvalCounter++;
      //      // restore solution
      //      solVec_ = solBak;
      //
      //      // SetupRes will write the residual directly to resVec_n_k_
      //      // It will also reassemble the system! > we have to call this function
      //      // again after the linesearch to get the curret system set up!
      //      if(iterationCounter == 0){
      //        WARN("Iteration counter should not be 0 at this point");
      //      }
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter
      //      gK = resVec_n_k_;
      //      /*
      //       * END RELAOD
      //       */
      gK = resInitial;

      tmpSol = currentSol;
      tmpSol.Add(alphaC,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      INEXACT_VectorBasedPolynomialEvalCounter++;
      solVec_ = solBak;

      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alphaC);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      gK_c = resVec_n_k_;

      /*
       * step 2: perform iteration for alphaK
       */
      Double scalingFactor;
      SBM_Vector gK_Minus_gK_alphaK;
      bool INEXACT_VectorBasedPolynomialSuccess = false;

      LOG_DBG2(solvestephyst_linesearch) << "Check sufficient decrease condition (qK_alphaK >= 0.5 + epsilon) for VectorBasedPolynomial: ";
      for(k = 0; k < VecPolyParameter.maxIter_; k++){
        scalingFactor = alphaK/alphaC;
        if(scalingFactor == 1){
          gK_alphaK = gK_c;
        } else {
          gK_alphaK = zeroVec_;
          gK_alphaK.Add(scalingFactor,gK_c,1.0-scalingFactor,gK);
        }

        gK_Minus_gK_alphaK = zeroVec_;
        gK_Minus_gK_alphaK.Add(1.0,gK,-1.0,gK_alphaK);

        if(gK_Minus_gK_alphaK.NormL2() == 0){
          EXCEPTION("gK_Minus_gK_alphaK.NormL2() should not become 0!");
        }

        qK_alphaK = 0;
        gK.Inner(gK_Minus_gK_alphaK,qK_alphaK);
        qK_alphaK /= (gK_Minus_gK_alphaK.NormL2()*gK_Minus_gK_alphaK.NormL2());

        LOG_DBG2(solvestephyst_linesearch) << "Iteration " << k << "; alphaK = " << alphaK << "; qK_alphaK = " << qK_alphaK << "; 0.5 + epsilon = " << 0.5 + epsilon;

        // check stopping criterion
        if(qK_alphaK >= 0.5 + epsilon){
          INEXACT_VectorBasedPolynomialSuccess = true;
          break;
        } else {
          /* determine new alphaK based on h(x)
           *  1             if x >= 0.5 + epsilon => stopping criterion
           *  x             if epsilon <= x < 0.5 + epsilon
           *  epsilon       if 0 <= x < epsilon
           * -epsilon       if -epsilon <= x < 0
           * x              if -0.5-epsilon <= x < -epsilon
           * -0.5-epsilon   if x < -0.5-epsilon
           */
          Double h = 1;
          if(qK_alphaK >= epsilon){
            h = qK_alphaK;
          } else if(qK_alphaK >= 0){
            h = epsilon;
          } else if(qK_alphaK >= -epsilon){
            h = -epsilon;
          } else if(qK_alphaK >= -0.5-epsilon){
            h = qK_alphaK;
          } else {
            h = -0.5-epsilon;
          }
          alphaK *= h;
        }
      }

      if(INEXACT_VectorBasedPolynomialSuccess){
        // as we do not evaluate the residual itself during the iterations, we have
        // to do it afterwards (only for this debugging messsage, though)
        tmpSol = currentSol;
        tmpSol.Add(alphaK,solIncrement_n_k_local);
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        INEXACT_VectorBasedPolynomialEvalCounter++;
        solVec_ = solBak;

        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, alphaK);
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      }

      if(VecPolyParameter.measurePerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for VectorBasedPolynomial-Linesearch ===" << std::endl;
        //        std::cout << "Number of iteraion steps: " << k << std::endl;
        //        std::cout << "Number of hyst evaluations: " << INEXACT_VectorBasedPolynomialEvalCounter << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for VectorBasedPolynomial-Linesearch ===";
        LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
        LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << INEXACT_VectorBasedPolynomialEvalCounter;
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      } else {
        endTime = 0.0;
      }

      // flag indicating if at least one linesearch was performed (this is the case at this point)
      // this does not mean, that linesearch was successful!
      validLSperformed = true;

      /*
       * check if linesearch was successful, too
       */
      if((alphaK < 0)&&(!VecPolyParameter.allowNegativeEta_)){
        validEtaFound = false;
        INEXACT_VectorBasedPolynomialSuccess = false;
      }

      if(INEXACT_VectorBasedPolynomialSuccess){
        validEtaFound = true;
      }

//      // valid eta found but maybe too small!
//      if(validEtaFound){
//        if(abs(alphaK) < VecPolyParameter.minEta_){
//          if(alphaK < 0){
//            alphaK = -VecPolyParameter.minEta_;
//          } else {
//            alphaK = VecPolyParameter.minEta_;
//          }
//        }
//      }

      LSResults curMethod;
      curMethod.name = "INEXACT_VectorBasedPolynomial";
      curMethod.lsFlag = LS_VECBASEDPOLY;
      curMethod.success = INEXACT_VectorBasedPolynomialSuccess;
      curMethod.alphaFound = alphaK;
      curMethod.resultingRes = resVec_n_k_.NormL2()/referenceValue;
      curMethod.totalEvalTime = endTime-startTime;
      curMethod.numResidualEvals = INEXACT_VectorBasedPolynomialEvalCounter;

      testedLSMethods.push(curMethod);
      etaReturn = alphaK;
    }

    if( (specifiedLinesearches_.find(LS_THREEPOINTPOLY) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "INEXACT_ThreePointPolynomial") ){
      /*
       * Source: Iterative Methods for Linear and Nonlinear Equations - C.T. Kelly SIAM 1995
       * (Three-Point-Parabolic-Model page 143)
       *
       * INEXACT LINESEARCH, i.e. take largest steplength that satisfies sufficient decrease
       *
       * Notation:
       * F(lambda) = r(lambda) = |residual(x_k + lambda*dx)|
       * f(lambda) = r(lambda)^2 = = |residual(x_k + lambda*dx)|^2
       *
       *
       */
      /*
       * Load read-in parameter
       */
      LinesearchParameter ThreePointPolyParameter = specifiedLinesearchesMap_[LS_THREEPOINTPOLY];

      LOG_TRACE(solvestephyst_linesearch) << "=== Perform ThreePointPolynomial linesarch ===";
      if(ThreePointPolyParameter.measurePerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      UInt ThreePointEvalCounter = 0;
      k = 0;

      Double sigma0,sigma1;

      /*
       * Factor for sufficient decrease
       * > default: take value as suggested by Kelly (1e-4)
       */
      Double alphaCheck = ThreePointPolyParameter.alphaCheck_;
      bool ThreePointSuccess = false;

      /*
       * Safeguarding parameter for ste increase/decrease
       * > take values as suggested by Kelly
       */
      sigma0 = ThreePointPolyParameter.sigma0_;
      sigma1 = ThreePointPolyParameter.sigma1_;

      // residual norms and squared residual norms
      Double r0,rCur; //,rPrev;
      Double f0,fCur,fPrev;

      // step values
      Double lambdaCur = 1.0;
      Double lambdaPrev = -1.0;

      /*
       * Step 1: a) Compute f(lambda) for full step, i.e. lambda = 1
       *          and for zero step, i.e. previous residual
       *         b) check if full step is acceptable
       */
      //      /*
      //       * zero step should be loadable from beginning of linesearch as residual is precomputed
      //       *
      //       * TODO: BEGIN LOAD
      //       */
      //      tmpSol = currentSol;
      //      solVec_ = tmpSol;
      //      // set solution vector to tmp solution, then evaluate hyst operators
      //      // do not reevaluate matrices for inversion > flag = 0
      //      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      //      ThreePointEvalCounter++;
      //      solVec_ = solBak;
      //
      //      // SetupRes will write the residual directly to resVec_n_k_
      //      // It will also reassemble the system! > we have to call this function
      //      // again after the linesearch to get the curret system set up!
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      //      r0 = resVec_n_k_.NormL2();
      //      /*
      //       * END LOAD
      //       */
      r0 = resInitial.NormL2();
      f0 = r0*r0;

      tmpSol = currentSol;
      tmpSol.Add(1.0,solIncrement_n_k_local);
      solVec_ = tmpSol;
      // set solution vector to tmp solution, then evaluate hyst operators
      // do not reevaluate matrices for inversion > flag = 0
      PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
      ThreePointEvalCounter++;
      solVec_ = solBak;

      // SetupRes will write the residual directly to resVec_n_k_
      // It will also reassemble the system! > we have to call this function
      // again after the linesearch to get the curret system set up!
      SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, 1.0);
      //      SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
      rCur = resVec_n_k_.NormL2();
      fCur = rCur*rCur;

      /*
       * Check if full step is acceptable
       */
      LOG_DBG2(solvestephyst_linesearch) << "Check sufficient decrease condition (rCur < (1-alphaCheck*lambdaCur)*r0) for ThreePointPolynomial: ";
      LOG_DBG2(solvestephyst_linesearch) << "lambdaCur - rCur - (1-alphaCheck*lambdaCur)*r0";
      LOG_DBG2(solvestephyst_linesearch) << lambdaCur << " - " << rCur << " - " << (1-alphaCheck*lambdaCur)*r0;

      if(rCur < (1-alphaCheck)*r0){
        ThreePointSuccess = true;
      } else {
        //        rPrev = rCur;
        fPrev = fCur;
        lambdaPrev = lambdaCur;

        /*
         * Evaluate at sigma1 and check again
         */
        lambdaCur = sigma1;

        tmpSol = currentSol;
        tmpSol.Add(lambdaCur,solIncrement_n_k_local);
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        ThreePointEvalCounter++;
        solVec_ = solBak;

        // SetupRes will write the residual directly to resVec_n_k_
        // It will also reassemble the system! > we have to call this function
        // again after the linesearch to get the curret system set up!
        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, lambdaCur);
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        rCur = resVec_n_k_.NormL2();
        fCur = rCur*rCur;

        LOG_DBG2(solvestephyst_linesearch) << lambdaCur << " - " << rCur << " - " << (1-alphaCheck*lambdaCur)*r0;
        if(rCur < (1-alphaCheck*lambdaCur)*r0){
          ThreePointSuccess = true;
        }
      }

      /*
       * Still no success > begin with polynomial calculation
       */

      // first and second derivative of polynomial
      Double pd1,pd2;
      Double term1,term2,term3;
      Double lambdaOpt;
      while(ThreePointSuccess == false){
        k++;

        /*
         * p(lambda) = f(0) + lambda/[lambdaCur*lambdaPrev*(lambdaCur - lambdaPrev)] *
         *            [ lambdaPrev*(lambda-lambdaPrev)*(fCur-f0) +
         *              lambdaCur*(lambdaCur-lambda)*(fPrev-f0) ]
         *
         * pd1(lambda) = 1/[lambdaCur*lambdaPrev*(lambdaCur - lambdaPrev)] *
         *            [ lambdaPrev*(2*lambda-lambdaPrev)*(fCur-f0) +
         *              lambdaCur*(lambdaCur-2*lambda)*(fPrev-f0) ]
         *
         * pd1(0) = 1/[lambdaCur*lambdaPrev*(lambdaCur - lambdaPrev)] *
         *            [ -lambdaPrev*lambdaPrev*(fCur-f0) +
         *               lambdaCur*lambdaCur*(fPrev-f0) ]
         *
         * pd2(lambda) = 1/[lambdaCur*lambdaPrev*(lambdaCur - lambdaPrev)] *
         *            [ 2*lambdaPrev*(fCur-f0) - 2*lambdaCur*(fPrev-f0) ]
         *
         * lambdaMinimizer = -pd1(0)/pd2(0) if pd2 > 0
         */

        /*
         * Compute pd2(0)
         */
        term1 = 1.0/(lambdaCur*lambdaPrev*(lambdaCur-lambdaPrev));
        term2 = lambdaPrev*(fCur-f0);
        term3 = lambdaCur*(fPrev-f0);

        pd2 = 2*term1*(term2 - term3);

        if(pd2 > 0){
          /*
           * compute pd1 and minimizer lampdaOpt
           */
          pd1 = term1*(-lambdaPrev*term2 + lambdaCur*term3);
          lambdaOpt = -pd1/pd2;

          //          std::cout << "LambdaOpt: " << lambdaOpt << std::endl;

          /*
           * apply safeguarding to lambdaOpt to get new lambdaCur
           */
          if(lambdaOpt < sigma0*lambdaCur){
            lambdaOpt = sigma0*lambdaCur;
          } else if(lambdaOpt > sigma1*lambdaCur){
            lambdaOpt = sigma1*lambdaCur;
          }

        } else {
          /*
           * parabolic model does not lead to minimizer
           * > take sigma1*lambdaCur as suggested by Kelly
           */
          lambdaOpt = sigma1*lambdaCur;
        }

        lambdaPrev = lambdaCur;
        //        rPrev = rCur;
        fPrev = fCur;
        lambdaCur = lambdaOpt;

        /*
         * Check lambdaOpt
         */
        tmpSol = currentSol;
        tmpSol.Add(lambdaCur,solIncrement_n_k_local);
        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        ThreePointEvalCounter++;
        solVec_ = solBak;

        // SetupRes will write the residual directly to resVec_n_k_
        // It will also reassemble the system! > we have to call this function
        // again after the linesearch to get the curret system set up!
        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, lambdaCur);
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        rCur = resVec_n_k_.NormL2();
        fCur = rCur*rCur;

        LOG_DBG2(solvestephyst_linesearch) << lambdaCur << " - " << rCur << " - " << (1-alphaCheck*lambdaCur)*r0;

        if(rCur < (1-alphaCheck*lambdaCur)*r0){
          ThreePointSuccess = true;
        }

        if(k == ThreePointPolyParameter.maxIter_){
          /*
           * maximal number of iterations reached
           */
//          std::cout << "ThreePoint: max number iterations exceeded" << std::endl;
//          std::cout << "rCur = " << rCur << std::endl;
//          std::cout << "r0 = " << r0 << std::endl;
//          std::cout << "alphaCheck = " << alphaCheck << std::endl;
//          std::cout << "lambdaCur = " << lambdaCur << std::endl;
//          std::cout << "(1-alphaCheck*lambdaCur)*r0 = " << (1-alphaCheck*lambdaCur)*r0 << std::endl;
          break;

        }
      }

      if(ThreePointPolyParameter.measurePerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for ThreePointPolynomial-Linesearch ===" << std::endl;
        //        std::cout << "Number of iteraion steps: " << k << std::endl;
        //        std::cout << "Number of hyst evaluations: " << ThreePointEvalCounter << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for ThreePointPolynomial-Linesearch ===";
        LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
        LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << ThreePointEvalCounter;
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      } else {
        endTime = 0.0;
      }

      // flag indicating if at least one linesearch was performed (this is the case at this point)
      // this does not mean, that linesearch was successful!
      validLSperformed = true;

      /*
       * check if linesearch was successful, too
       */
      if((lambdaCur < 0)&&(!ThreePointPolyParameter.allowNegativeEta_)){
        validEtaFound = false;
        ThreePointSuccess = false;
      }

      if(ThreePointSuccess){
        validEtaFound = true;
      }

      LSResults curMethod;
      curMethod.name = "INEXACT_ThreePointPolynomial";
      curMethod.lsFlag = LS_THREEPOINTPOLY;
      curMethod.success = ThreePointSuccess;
      curMethod.alphaFound = lambdaCur;
      curMethod.resultingRes = rCur/referenceValue;
      curMethod.totalEvalTime = endTime-startTime;
      curMethod.numResidualEvals = ThreePointEvalCounter;

      testedLSMethods.push(curMethod);

      etaReturn = lambdaCur;
    }

    if( (specifiedLinesearches_.find(LS_BACKTRACKING_ARMIJO) != specifiedLinesearches_.end()) ||
        (specifiedLinesearches_.find(LS_BACKTRACKING_ARMIJO_NONMONOTONIC) != specifiedLinesearches_.end()) ||
        (specifiedLinesearches_.find(LS_BACKTRACKING_GRADFREE) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "INEXACT_Backtracking") ){
      /*
       * Sources:
       * Armijo:
       *  Optimization Theory and Methods - Nonlinear Programming -- Sun, Yuan -- Springer 2006
       *  Limited memory BFGS method with backtracking for symmetric nonlinear equations -- Yuan, Wei, Lu -- 2011
       *  (equ. (1.6) on page 368
       * Non-monotonic-Armijo:
       *  Limited memory BFGS method with backtracking for symmetric nonlinear equations -- Yuan, Wei, Lu -- 2011
       *  (equ. (1.7) on page 368
       * Gradient-free:
       *  Limited memory BFGS method with backtracking for symmetric nonlinear equations -- Yuan, Wei, Lu -- 2011
       *  (equ. (1.8) on page 368
       *
       * Common idea:
       *  - inexact linesearch, i.e. stop if certain estimations/criterions are satisfied
       *  - backtracking: start at eta=etaStart(usually 1 for Newton) and decresase until
       *      criterion(s) is(are) satisfied
       *  - objective function = residual norm
       *
       * Armijo:
       *  Stop if
       *    f(x + eta*dx) <= f(x) + rho*eta*gradF(x)^T*dx
       *
       *    f(x) = 0.5*residual(x).NormL2*residual(x).NormL2
       *    gradF(x) = gradResidual(x)*residual(x) = Jacobian(x)*residual(x)
       *    > for K_eff*x = f(x) (i.e. nonlinearity only on rhs)
       *      Our Nomenclature
       *        residual_our(x) = f(x) - K_eff*x
       *        Jacobian_our = df/dx - K_eff = -1.0*Systemmatrix in case of deltaMat formulation
       *      Their Nomenclature
       *        residual_th(x) = K_eff*x - f(x)
       *        Jacobian_th = K_eff - df/dx
       *
       *     > we solve:
       *        K_eff*dx = f(x) - K_eff*x = residual_our(x) = -residual_th(x)
       *        or in deltaMat/Jacobian case
       *        (K_eff + df/dx)*dx = f(x) - K_eff*x
       *        >  Jacobian_th = residual_our(x) = -residual_th(x)
       *        > the formulas are based on their nomenclature, so if we have
       *            Jacobian_th*residual_th it is the same as Jacobian_our*residual_our,
       *            but
       *            residual_th*dx = -1.0*residual_our*dx       *
       *
       * Non-monotonic-Armijo:
       *  Stop if
       *    f(x + eta*dx) <= f_max_M(x) + rho*eta*gradF(x)^T*dx
       *
       *    f(x) = 0.5*residual(x).NormL2*residual(x).NormL2
       *    f_max_M(x) = maximal residual(x).NormL2 from the last M outer iterations
       *
       * Gradient free:
       *  Stop if
       *    residual(x+eta*dx).NormL2*residual(x+eta*dx).NormL2 <= residual(x).NormL2*residual(x).NormL2 + rho*eta*eta*residual(x)^T*dx
       *
       *
       */
      std::stack< std::string > modes;
//      if(useAllLinsearches_){
//        modes.push("Armijo");
//        modes.push("non-monotonic-Armijo");
//        modes.push("gradientFree");
//      } else {
      if( specifiedLinesearches_.find(LS_BACKTRACKING_ARMIJO) != specifiedLinesearches_.end()){
        modes.push("Armijo");
      }
      if( specifiedLinesearches_.find(LS_BACKTRACKING_ARMIJO_NONMONOTONIC) != specifiedLinesearches_.end()){
        modes.push("non-monotonic-Armijo");
      }
      if( specifiedLinesearches_.find(LS_BACKTRACKING_GRADFREE) != specifiedLinesearches_.end()){
        modes.push("gradientFree");
      }

      UInt maxIter;
      Double eta,etaMax,etaMin;
      Double rho,omega;
//      omega = 0.5; // decrease factor
//      rho = 0.9; // stopping tolerance

      std::string mode;
      linesearchFlag lsFlagLocal;
      while(!modes.empty()){
        mode = modes.top();

        if(mode == "Armijo"){
          lsFlagLocal = LS_BACKTRACKING_ARMIJO;
        } else if(mode == "non-monotonic-Armijo"){
          lsFlagLocal = LS_BACKTRACKING_ARMIJO_NONMONOTONIC;
        } else if(mode == "gradientFree"){
          lsFlagLocal = LS_BACKTRACKING_GRADFREE;
        } else {
          lsFlagLocal = INVALID_LINESEARCH;
        }

        /*
         * Load read-in parameter
         */
        LinesearchParameter ArmijoParameter = specifiedLinesearchesMap_[lsFlagLocal];
        maxIter = ArmijoParameter.maxIter_;
        rho = ArmijoParameter.rho_;
        omega = ArmijoParameter.decreaseFactor_;

        LOG_TRACE(solvestephyst_linesearch) << "=== Perform Backtracking linesarch based on Armijo-type rules ===";
        LOG_TRACE(solvestephyst_linesearch) << "=== Armijo-Type criterion: "<<mode<<" ===";

        if(ArmijoParameter.measurePerformance_){
          startTime = LinesearchTimer->GetCPUTime();
          LinesearchTimer->Start();
        }

        UInt INEXACT_BacktrackingEvalCounter = 0;

        etaMax = ArmijoParameter.etaStart_;
        etaMin = ArmijoParameter.minEta_;
        eta = etaMax;

        std::stringstream fullName;
        fullName << "INEXACT_Backtracking (mode = " << mode << ")";
        bool successINEXACT_Backtracking = false;

        /*
         * get current residual
         */
        SBM_Vector resK;
        Double resOptINEXACT_Backtracking = 1e15;
        Double alphaINEXACT_Backtracking = 1.0;
        Double resNormOLD;
        Double resNormNEW = 1e15;

        //        /*
        //         * same as initial residual
        //         * > LOAD
        //         *
        //         * TODO: BEGIN LOAD
        //         */
        //        tmpSol = currentSol;
        //        solVec_ = tmpSol;
        //        // set solution vector to tmp solution, then evaluate hyst operators
        //        // do not reevaluate matrices for inversion > flag = 0
        //        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        //        INEXACT_BacktrackingEvalCounter++;
        //        solVec_ = solBak;
        //
        //        // SetupRes will write the residual directly to resVec_n_k_
        //        // It will also reassemble the system! > we have to call this function
        //        // again after the linesearch to get the curret system set up!
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        //        resK = resVec_n_k_;
        //        /*
        //         * END LOAD
        //         */
        resK = resInitial;
        resNormOLD = resK.NormL2();
        resNormOLD *= resNormOLD;

        if(mode == "non-monotonic-Armijo"){
//          std::cout << "Current history for non-monotonic-Armijo" << std::endl;
          LOG_DBG2(solvestephyst_linesearch) << "--- non-monotonic-Armijo > old residual norm = maximum of previous " << resNormHistoryArmijo_.size() << " residual norms";
          resNormHistoryArmijo_.push_back(resNormOLD);
          resNormHistoryArmijo_.pop_front();
          for(UInt ii = 0; ii < resNormHistoryArmijo_.size(); ii++ ){
//            std::cout << "resNormHistoryArmijo_["<<ii<<"] = " << resNormHistoryArmijo_[ii] << std::endl;
            LOG_DBG2(solvestephyst_linesearch) << "--- resNormHistoryArmijo_["<<ii<<"] = " << resNormHistoryArmijo_[ii];
            if(resNormHistoryArmijo_[ii] > resNormOLD){
              resNormOLD = resNormHistoryArmijo_[ii];
            }
          }
//          std::cout << "Maximum of current history = " << resNormOLD << std::endl;
//          std::cout << ">> residual = sqrt(Maximum of current history) = " << std::sqrt(resNormOLD) << std::endl;
        }

        LOG_DBG2(solvestephyst_linesearch) << "--- old residual norm (squared): "<<resNormOLD;

        /*
         * Precompute products with solIncrement_n_k_local and residual
         */
        Double productTerm = 0.0;
        if(mode == "gradientFree"){
          // here we need residual^T*solIncrement_n_k_local
          resK.Inner(solIncrement_n_k_local,productTerm);
          // see note above for Armijo; we have flipped sign of residual
          // if we multiply with Jacobian it is ok (as this one has flipped sign, too)
          // but if me multiply by solIncrement_n_k_local, we have to consider the changed sign
          productTerm *= -1.0;
//          std::cout << "precomputed productTerm = Residual^T*solincrement = " << productTerm << std::endl;
          LOG_DBG2(solvestephyst_linesearch) << "precomputed productTerm = Residual^T*solincrement = " << productTerm;
        } else {
          SBM_Vector JacRes;
          JacRes = resVec_n_k_;

          solVec_ = currentSol;

          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
          INEXACT_BacktrackingEvalCounter++;

          if(currentBuildFlag_ != FDJACOBIAN){
            LOG_DBG(solvestephyst_linesearch) << "--- FD-Jacobian needed for armijo rule; assembly of system required ---";

            // build up system using FD-Jacobian
            AssembleSystem(FDJACOBIAN,currentSol,false);//, -1);
            // do not backup temp system
            BuildSystem(true);
          }
          algsys_->ComputeSysMatTimesVector(SYSTEM,resK,JacRes,false);

          solVec_ = solBak;

          // now compute JacRes^T*solIncrement_n_k_remnt
          JacRes.Inner(solIncrement_n_k_local,productTerm);
          productTerm *= -1.0; // remember: our systemmatrix = -1 jacobian

          LOG_DBG2(solvestephyst_linesearch) << "precomputed productTerm = (Jacobian*Residual)^T*solincrement = " << productTerm;
        }

        LOG_DBG2(solvestephyst_linesearch) << "Check sufficient decrease condition (resNormNEW <= criterion) for Backtracking: ";
        for(k = 0; k < maxIter; k++){
          // compute current residual
          tmpSol = currentSol;
          tmpSol.Add(eta,solIncrement_n_k_local);
          solVec_ = tmpSol;
          // set solution vector to tmp solution, then evaluate hyst operators
          // do not reevaluate matrices for inversion > flag = 0
          PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
          INEXACT_BacktrackingEvalCounter++;
          solVec_ = solBak;
          // SetupRes will write the residual directly to resVec_n_k_
          // It will also reassemble the system! > we have to call this function
          // again after the linesearch to get the curret system set up!
          SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, eta);
          //          SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
          resNormNEW = resVec_n_k_.NormL2();
          resNormNEW *= resNormNEW;

          // define criterion
          Double criterion = resNormOLD;
          if(mode == "gradientFree"){
            criterion += rho*eta*eta*productTerm;
//            std::cout << "eta = "<<eta<<"; resNormNEW = " << resNormNEW << "; criterion = resNormOLD + rho*eta*eta*productTerm = " << criterion << std::endl;
            LOG_DBG2(solvestephyst_linesearch) << "eta = "<<eta<<"; resNormNEW = " << resNormNEW << "; criterion = rho*eta*eta*productTerm = " << criterion;
          } else {
            // factor 2.0 would be needed as we work with res*res instead of 1/2 res*res
            criterion += rho*eta*productTerm;
            LOG_DBG2(solvestephyst_linesearch) << "eta = "<<eta<<"; resNormNEW = " << resNormNEW << "; criterion = rho*eta*productTerm = " << criterion;
          }

          if(resNormNEW <= criterion){
//            std::cout << "SUCCESS for eta = "<<eta<< std::endl;
            successINEXACT_Backtracking = true;
            break;
          } else {
            eta = omega*eta;

            if(eta < etaMin){
//              std::cout << "STOP as eta = "<<eta<< " < etaMin = " << etaMin <<std::endl;
//              std::cout << "BacktrackingLS - eta < etaMin (" << minEta_ << ") > return etaMin" << std::endl;
              break;
            }
          }
        }

        resOptINEXACT_Backtracking = std::sqrt(resNormNEW);
        alphaINEXACT_Backtracking = eta;

        if(ArmijoParameter.measurePerformance_){
          LinesearchTimer->Stop();
          endTime = LinesearchTimer->GetCPUTime();

          //          std::cout << "=== Statistics for "<<fullName.str()<<" ===" << std::endl;
          //          std::cout << "Number of iteraion steps: " << k << std::endl;
          //          std::cout << "Number of hyst evaluations: " << INEXACT_BacktrackingEvalCounter << std::endl;
          //          std::cout << "Required runtime: " << endTime-startTime << std::endl;

          LOG_DBG(solvestephyst_linesearch) << "=== Statistics for "<<fullName.str()<<" ===";
          LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
          LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << INEXACT_BacktrackingEvalCounter;
          LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
        } else {
          endTime = 0.0;
        }

        // flag indicating if at least one linesearch was performed (this is the case at this point)
        // this does not mean, that linesearch was successful!
        validLSperformed = true;

        /*
         * check if linesearch was successful, too
         */
        if((alphaINEXACT_Backtracking < 0)&&(!ArmijoParameter.allowNegativeEta_)){
          validEtaFound = false;
          successINEXACT_Backtracking = false;
        }

        if(successINEXACT_Backtracking){
          validEtaFound = true;
        }

        // valid eta found but maybe too small!
        if(validEtaFound){
          if(abs(alphaINEXACT_Backtracking) < ArmijoParameter.minEta_){
            if(alphaINEXACT_Backtracking < 0){
              alphaINEXACT_Backtracking = -ArmijoParameter.minEta_;
            } else {
              alphaINEXACT_Backtracking = ArmijoParameter.minEta_;
            }
          }
        }

        LSResults curMethod;
        curMethod.name = fullName.str();
        curMethod.lsFlag = lsFlagLocal;
        curMethod.success = successINEXACT_Backtracking;
        curMethod.alphaFound = alphaINEXACT_Backtracking;
        curMethod.resultingRes = resOptINEXACT_Backtracking/referenceValue;
        curMethod.totalEvalTime = endTime-startTime;
        curMethod.numResidualEvals = INEXACT_BacktrackingEvalCounter;

        testedLSMethods.push(curMethod);

        etaReturn = alphaINEXACT_Backtracking;

        modes.pop();
      }

    }

    if( (specifiedLinesearches_.find(LS_TRIAL_AND_ERROR) != specifiedLinesearches_.end()) ||
        (lineSearch_ == "INEXACT_TrialAndError") ){
      /*
       * trial and error backtracking, i.e. decrease eta and take the one that minimizes the residual
       */
      LOG_TRACE(solvestephyst_linesearch) << "=== Perform TrialAndError linesarch ===";

            /*
       * Load read-in parameter
       */
      LinesearchParameter TrialAndErrorParameter = specifiedLinesearchesMap_[LS_TRIAL_AND_ERROR];

      bool successClassicApproach = false;

      if(TrialAndErrorParameter.measurePerformance_){
        startTime = LinesearchTimer->GetCPUTime();
        LinesearchTimer->Start();
      }

      UInt INEXACT_TrialAndErrorEvalCounter = 0;

      /*
       * start with eta = 1;
       * test residual
       *  > success if res-norm < c1*resNormOld
       *  > success if res-norm < stoppingTol
       *  > stop if res-norm increase twice in a row
       *  > otherwise: eta = eta/c2
       */
      Double c1 = TrialAndErrorParameter.residualDecreaseForSuccess_; //  0.1; //residualDecreaseForSuccess
      Double c2 = TrialAndErrorParameter.decreaseFactor_; // 0.5; //decreaseFactor
      Double eta = TrialAndErrorParameter.etaStart_; // 1.0;
      Double bestNorm = 1e16;
      Double bestEta = 1.0;
      Double resNorm = 0.0;
      Double resNormPrev = 1e16;
      //Double minEta = 1/256.0;
      Double minEta = TrialAndErrorParameter.minEta_; //  std::max(1/64.0,minEta_);
      Double resOptClassic = 1e15;
      Double alphaClassic = 1.0;
      Double stoppingTol = TrialAndErrorParameter.stoppingTol_; // 1e-4;
      residualErrAbsPrev_ = resInitial.NormL2();
      bool tryOverrelaxation = TrialAndErrorParameter.testOverrelaxation_;
      bool revertedDirection = false;

      // new: define one fixed reference point for relative error calculation
      //      Double normOfSystemFix = SetupRESVector(currentSol,currentSol,currentSol,iterationCounter);

      UInt cnt = 0;
      k = 0;
      SBM_Vector tmpSol;

      LOG_DBG2(solvestephyst_linesearch) << "Check decrease condition (resNorm < c1*residualErrAbsPrev_) for TrialAndError: ";
//      std::cout << "Initial residual: "<<residualErrAbsPrev_<< std::endl;
      while(eta > minEta){
//        std::cout << "Trial and error: current eta = " << eta << std::endl;
        tmpSol = currentSol;
        tmpSol.Add(eta,solIncrement_n_k_local);

        solVec_ = tmpSol;
        // set solution vector to tmp solution, then evaluate hyst operators
        // do not reevaluate matrices for inversion > flag = 0
        PDE_.SetFlagInCoefFncHyst("EvaluateHystOperators",0);
        solVec_ = solBak;

        // SetupRes will write the residual directly to resVec_n_k_
        // It will also reassemble the system! > we have to call this function
        // again after the linesearch to get the curret system set up!
        // NOTE: normOfSystem changes with eta!
        SetupRESVector_Linesearch(flag_for_matrix_assembly, flag_for_vector_assembly,resVec_n_k_,tmpSol, tmpSol, currentSol, solIncrement_n_k_local, eta);
        //        SetupRESVector(tmpSol,tmpSol,tmpSol,iterationCounter);
        INEXACT_TrialAndErrorEvalCounter++;
        resNorm = resVec_n_k_.NormL2();

//        std::cout << "resNorm: "<<resNorm<< std::endl;
        if(resNorm < bestNorm){
          // abs values
          bestNorm = resNorm;
          bestEta = eta;
        }

//        if((useRelativeNormForRes_)&&(normOfSystemFix != 0)){
//          resNorm = resNorm/normOfSystemFix;
//        }

        LOG_DBG2(solvestephyst_linesearch) << "eta = "<<eta<<"; resNorm = "<<resNorm<<"; c1*residualErrAbsPrev_ "<<c1*residualErrAbsPrev_;

        if(resNorm < c1*residualErrAbsPrev_){

//          std::cout << "resNorm ("<<resNorm<<") < c1 ("<<c1<<")*residualErrAbsPrev_ ("<< residualErrAbsPrev_ <<")" << std::endl;

          resOptClassic = resNorm;
          alphaClassic = eta;
          break;
          //return eta;
        }
        if((resNorm < stoppingTol)&&(resNorm < 0.9*residualErrAbsPrev_ )){
//          std::cout << "resNorm ("<<resNorm<<") < stoppingTol ("<<stoppingTol<<") && 0.9*residualErrAbsPrev_" << std::endl;
          resOptClassic = resNorm;
          alphaClassic = eta;
          break;
        }


        UInt Ntimes = 2;
        if ((resNorm >= resNormPrev)){
          cnt++;
          // if residual is larger than previous one at least N-times in a row,
          // consider that we step away from optimum
          // case a) if we already hat some decrease in residual, take the current value; it seems to be a local min
          // case b) optimal stepping length so far is etaStart, maybe overrelaxtion can help > change c2 to 1/c2
          //
          if(cnt >= Ntimes){
//            LOG_DBG2(solvestephyst_linesearch) << "resNorm increased twice in a row; take bestEta = "<< bestEta;
            resOptClassic = bestNorm;
            alphaClassic = bestEta;

//            std::cout << "Residual increased twice in a row" << std::endl;
            if( (tryOverrelaxation == true)&&(bestEta == TrialAndErrorParameter.etaStart_) ){
//              std::cout << "Currently best eta ("<<bestEta<<") is still etaStart; maybe overrelaxation helps" << std::endl;
              eta = bestEta;
              c2 = 1.0/c2;
//              tryOverrelaxation = true;
              if(revertedDirection == true){
                // flip direction only once! otherwise we could do this over and over again until kmax is exceeded
                break;
              } else {
                revertedDirection = true;
              }
            } else {
//              std::cout << "We already had an improvement in residual; take best eta" << std::endl;
              break;
            }
          }
        } else {
          cnt = 0;
        }

        resNormPrev = resNorm;
        eta = eta*c2;
        k++;

        if(k >= TrialAndErrorParameter.maxIter_){

//          std::cout << "k ("<<k<<") >= TrialAndErrorParameter.maxIter_ ("<<TrialAndErrorParameter.maxIter_<<")" << std::endl;

          break;
        }
      }

//      std::cout << "Taking best eta ("<<bestEta<<") leading to smallest absolute residual ("<<bestNorm<<")" << std::endl;

      bool failAtMinEta = false;
      if(bestEta > minEta){
        successClassicApproach = true;
      } else if(!failAtMinEta) {
        successClassicApproach = true;
      } else {
        successClassicApproach = false;
      }

      resOptClassic = bestNorm;
      alphaClassic = bestEta;

      if(TrialAndErrorParameter.measurePerformance_){
        LinesearchTimer->Stop();
        endTime = LinesearchTimer->GetCPUTime();

        //        std::cout << "=== Statistics for TrialAndErrorEval ===" << std::endl;
        //        std::cout << "Number of iteraion steps: " << k << std::endl;
        //        std::cout << "Number of hyst evaluations: " << INEXACT_TrialAndErrorEvalCounter << std::endl;
        //        std::cout << "Required runtime: " << endTime-startTime << std::endl;

        LOG_DBG(solvestephyst_linesearch) << "=== Statistics for TrialAndErrorEval ===";
        LOG_DBG(solvestephyst_linesearch) << "Number of iteraion steps: " << k;
        LOG_DBG(solvestephyst_linesearch) << "Number of hyst evaluations: " << INEXACT_TrialAndErrorEvalCounter;
        LOG_DBG(solvestephyst_linesearch) << "Required runtime: " << endTime-startTime;
      } else {
        endTime = 0.0;
      }

      // flag indicating if at least one linesearch was performed (this is the case at this point)
      // this does not mean, that linesearch was successful!
      validLSperformed = true;

      /*
       * check if linesearch was successful, too
       */
      if((alphaClassic < 0)&&(!TrialAndErrorParameter.allowNegativeEta_)){
        validEtaFound = false;
        successClassicApproach = false;
      }

      if(successClassicApproach){
        validEtaFound = true;
      }

      // valid eta found but maybe too small!
      if(validEtaFound){
        if(abs(alphaClassic) < TrialAndErrorParameter.minEta_){
          if(alphaClassic < 0){
            alphaClassic = -TrialAndErrorParameter.minEta_;
          } else {
            alphaClassic = TrialAndErrorParameter.minEta_;
          }
        }
      }

      LSResults curMethod;
      curMethod.name = "INEXACT_TrialAndError";
      curMethod.lsFlag = LS_TRIAL_AND_ERROR;
      curMethod.success = successClassicApproach;
      curMethod.alphaFound = alphaClassic;
      curMethod.resultingRes = resOptClassic/referenceValue;
      curMethod.totalEvalTime = endTime-startTime;
      curMethod.numResidualEvals = INEXACT_TrialAndErrorEvalCounter;

      testedLSMethods.push(curMethod);

      etaReturn = alphaClassic;
    }

    delete LinesearchTimer;

    if(validLSperformed == false) {
      std::stringstream warnmsg;
      warnmsg << "Linesearch method " << lineSearch_ << " not implemented for hysteresis. ";
      warnmsg << "No linesearch was performed. Return eta = 1.0";
      WARN(warnmsg.str());
//      std::cout << warnmsg.str() << std::endl;
      return 1.0;
    }

    if(!validEtaFound) {
      std::stringstream warnmsg;
      if(currentSolutionApproach_ == SOLVE_JACOBIAN){
        // we already solved via FP Jacobian; if we return 0 as eta to trigger evaluation of FP Jacobian, we would
        // get the same problem again > loop till number of iterations exceeds; fail
        // > instead: take best found eta even if it did not pass the linsearch
        warnmsg << "Linesearch was not successful; solution was already computed via FD_Jacobian; take search direction with last found eta";
        WARN(warnmsg.str());
      } else {
//        warnmsg << "Linesearch was not successful; dismiss search direction. Return eta = 0.0; Trigger solution via FD_Jacobian";
        warnmsg << "Linesearch was not successful; dismiss search direction; trigger solution via FD_Jacobian";
        WARN(warnmsg.str());
        return 0.0;
      }
    }

    // if more than one linesearch was tested, we have to pick one (depending on choice in input.xml)
    // if DEBUG_compareLinesearches_ is true, we put out debug messages to cout
    if( (DEBUG_compareLinesearches_) || (testedLSMethods.size() > 1) ){
      if(DEBUG_compareLinesearches_){
        std::cout << "=== COMPARISON OF DIFFERENT LINESEARCHES ===" << std::endl;
      }
      LOG_TRACE(solvestephyst_linesearch) << "=== COMPARISON OF DIFFERENT LINESEARCHES ===";

      Double bestRes = 1e16;
      Double etaClosestToOne = 1e16;
      Double etaClosestToOneABS = 1e16;
      LSResults curMethod;
      LSResults bestMethod_Res;
      LSResults bestMethod_Eta;
      LSResults bestMethod_EtaABS;

      while(!testedLSMethods.empty()){
        curMethod = testedLSMethods.top();

        if(DEBUG_compareLinesearches_){
          if(curMethod.success){
            std::cout << "SUCCESSFUL == ";
            LOG_TRACE(solvestephyst_linesearch) << "SUCCESSFUL == ";
          } else {
            std::cout << "FAILED == ";
            LOG_TRACE(solvestephyst_linesearch) << "FAILED == ";
          }
          std::cout << curMethod.name << "\t: etaLS = " << curMethod.alphaFound << "\t; res = " << curMethod.resultingRes;
          std::cout << "\t; evals = " << curMethod.numResidualEvals << "\t; time = " << curMethod.totalEvalTime << std::endl;

          LOG_TRACE(solvestephyst_linesearch) << curMethod.name << "\t: etaLS = " << curMethod.alphaFound << "\t; res = " << curMethod.resultingRes;
          LOG_TRACE(solvestephyst_linesearch) << "\t; evals = " << curMethod.numResidualEvals << "\t; time = " << curMethod.totalEvalTime << std::endl;
        }

        if(curMethod.resultingRes < bestRes){
          if(curMethod.success){
            bestRes = curMethod.resultingRes;
            bestMethod_Res = curMethod;
          }
        }
        if( abs(1-curMethod.alphaFound) < abs(1-etaClosestToOne)){
          if(curMethod.success){
            etaClosestToOne = curMethod.alphaFound;
            bestMethod_Eta = curMethod;
          }
        }
        if( abs(1-abs(curMethod.alphaFound)) < abs(1-abs(etaClosestToOneABS))){
          if(curMethod.success){
            etaClosestToOneABS = curMethod.alphaFound;
            bestMethod_EtaABS = curMethod;
          }
        }

        testedLSMethods.pop();
      }

      std::string take_A = "";
      std::string take_B = "";
      std::string take_C = "";
      Double etaA,etaB,etaC;

      etaA = bestMethod_Res.alphaFound;;
      etaB = bestMethod_Eta.alphaFound;
      etaC = bestMethod_EtaABS.alphaFound;
      if(LSSelectionCriterion_ == 1){
        etaReturn = bestMethod_Res.alphaFound;
        take_A = "(SELECTED)";
        currentLinesearchApproach_ = bestMethod_Res.lsFlag;
      } else if(LSSelectionCriterion_ == 2){
        etaReturn = bestMethod_Eta.alphaFound;
        take_B = "(SELECTED)";
        currentLinesearchApproach_ = bestMethod_Eta.lsFlag;
      } else {
        etaReturn = bestMethod_EtaABS.alphaFound;
        take_C = "(SELECTED)";
        currentLinesearchApproach_ = bestMethod_EtaABS.lsFlag;
      }

      if(DEBUG_compareLinesearches_){
        std::cout << "=== A (minResidual) "<<take_A<<": eta = " << etaA << "; res = " << bestMethod_Res.resultingRes << std::endl;
        LOG_TRACE(solvestephyst_linesearch) << "=== A (minResidual) "<<take_A<<": eta = " << etaA << "; res = " << bestMethod_Res.resultingRes;

        std::cout << "=== B (eta closest to 1) "<<take_B<<": eta = " << etaB << "; res = " << bestMethod_Eta.resultingRes << std::endl;
        LOG_TRACE(solvestephyst_linesearch) << "=== B (eta closest to 1) "<<take_B<<": eta = " << etaB << "; res = " << bestMethod_Eta.resultingRes;

        std::cout << "=== C (abs(eta) closest to 1) "<<take_C<<": eta = " << etaC << "; res = " << bestMethod_EtaABS.resultingRes << std::endl;
        LOG_TRACE(solvestephyst_linesearch) << "=== C (abs(eta) closest to 1) "<<take_C<<": eta = " << etaC << "; res = " << bestMethod_EtaABS.resultingRes;
      }
    } else {
      // if only a single linesearch was performed, take it (at this point it was already checked that it was
      // successful as otherwise one of the two checks  if(validLSperformed == false) or if(!validEtaFound) would
      // have been triggered)
      LSResults lastTestedMethod;
      lastTestedMethod = testedLSMethods.top();
      currentLinesearchApproach_ = lastTestedMethod.lsFlag;
    }

    // already checked for single linesearches
//    if(abs(etaReturn)<minEta_){
//      if(etaReturn < 0){
//        etaReturn = -minEta_;
//      } else {
//        etaReturn = minEta_;
//      }
//      std::cout << "=== abs(etaSelected) < minEta (="<< minEta_ <<") -> take eta = "<<etaReturn<<" ===" << std::endl;
//      LOG_TRACE(solvestephyst_linesearch) << "=== abs(etaSelected) < minEta (="<< minEta_ <<") -> take eta = "<<etaReturn<<" ===";
//    } else {
//      std::cout << "=== abs(etaSelected) > minEta (="<< minEta_ <<") -> take eta = "<<etaReturn<<" ===" << std::endl;
//      LOG_TRACE(solvestephyst_linesearch) << "=== abs(etaSelected) > minEta (="<< minEta_ <<") -> take eta = "<<etaReturn<<" ===";
//    }
//
//    if((etaReturn < 0)&&(dismissNegativeEta_)){
//      std::cout << "=== DISMISS negative eta = "<<etaReturn<<"; return eta = 0.0 ===" << std::endl;
//      LOG_TRACE(solvestephyst_linesearch) << "=== DISMISS negative etaReturn = "<<etaReturn<<" ===";
//      return 0.0;
//    }

    return etaReturn;
  }

  Double SolveStepHyst::SetupRESVector_Linesearch(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
  SBM_Vector& targetVector, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForNonLinRHSAssembly,
          SBM_Vector& oldSolution, SBM_Vector& solutionUpdate, Double currentEta){

    LOG_TRACE(solvestephyst_residual) << "== SetupRESVector_Linesearch  ==";

    /*
     * basically the same as SetupRes but omit compututation of Keff*u and Keff*du
     * > only allowed if Keff is constant, i.e. not solution dependend in any way!
     */
    Double sysNorm = 0;

    // pass new solution to SetupRESorRHS
    // if Keff is somewhat nonlinear (even in the non-jacobian form), then compute
    // Keff*newSol inside that function
    // otherwise, sysNorm will stay 0
    // if the later is the case > build up Keff*newSol from precomputed values Keff*oldSol + eta*Keff*solIncrement
    SBM_Vector newSol;
    newSol = oldSolution;
    newSol.Add(currentEta,solutionUpdate);

    bool comingFromLinesearch = true;

//    std::cout << "SetupResVector_Linesearch" << std::endl;
//    std::cout << "Norm of targetVector at start: " << targetVector.NormL2() << std::endl;
    
    sysNorm = SetupRESorRHS(flag_for_matrix_assembly, flag_for_vector_assembly, newSol, solutionForMatrixAssembly, solutionForNonLinRHSAssembly,
      targetVector, comingFromLinesearch);

//    std::cout << "Norm of targetVector after SetupRESorRHS: " << targetVector.NormL2() << std::endl;
    
    /*
     * IMPORTANT: the following steps are ONLY allowed, if we compute the residual
     * with a linear system, i.e. all nonlinearities are contained in the rhs vectors!
     * if we have additional material nonlinearities in some regions, i.e. bh-saturation curve or so
     * we must recompute the system!
     * > This check is done in SetupRESorRHS; if no rebuild was done, sysNorm = 0 is returned
     */

    if(sysNorm == 0){
//      std::cout << "SetupRes - sysNorm = 0 > compute Keff*solution from Keff*oldSol and Keff*solUpdate" << std::endl;
      SBM_Vector systemUpdate;
      LOG_TRACE(solvestephyst_residual) << "=== Subtract K_eff*solutionVector (linesearch case, use precomputed values) ===";

      if(linSysMatrix_n_0_TIMES_solVec_n_kPrev_precomputed_ == false){
//        std::cout << "Precompute Keff*oldSolution" << std::endl;
        LOG_DBG(solvestephyst_residual) << "--- Precompute K_eff*previousSolutionVector ---";
        linSysMatrix_n_0_TIMES_solVec_n_kPrev_ = zeroVec_;

        Helper_ComputeKeffTimesSolution(oldSolution, linSysMatrix_n_0_TIMES_solVec_n_kPrev_,linearSystemBackuped_, flag_for_matrix_assembly);
        linSysMatrix_n_0_TIMES_solVec_n_kPrev_precomputed_ = true;
      }

      if(linSysMatrix_n_0_TIMES_solIncrement_n_k_precomputed_ == false){
//        std::cout << "Precompute Keff*solutionUpdate" << std::endl;
        LOG_DBG(solvestephyst_residual) << "--- Precompute K_eff*solutionUpdate ---";
        linSysMatrix_n_0_TIMES_solIncrement_n_k_ = zeroVec_;

        Helper_ComputeKeffTimesSolution(solutionUpdate, linSysMatrix_n_0_TIMES_solIncrement_n_k_,linearSystemBackuped_, flag_for_matrix_assembly);
        linSysMatrix_n_0_TIMES_solIncrement_n_k_precomputed_ = true;
      }

      systemUpdate = linSysMatrix_n_0_TIMES_solVec_n_kPrev_;
      systemUpdate.Add(currentEta,linSysMatrix_n_0_TIMES_solIncrement_n_k_);
      sysNorm = systemUpdate.NormL2();

      targetVector.Add(-1.0,systemUpdate);
    }
//    std::cout << "Norm of targetVector at end: " << targetVector.NormL2() << std::endl;
		return sysNorm;
	}

  Double SolveStepHyst::SetupRESorRHS(matrixFlag flag_for_matrix_assembly, vectorFlag flag_for_vector_assembly,
          SBM_Vector& solutionForUpdate, SBM_Vector& solutionForMatrixAssembly, SBM_Vector& solutionForRHSAssembly,
          SBM_Vector& returnVector, bool comingFromLinesearch){

    LOG_TRACE(solvestephyst_residual) << "====================================";
    LOG_TRACE(solvestephyst_residual) << "== SolveStepHyst - SetupRESorRHS  ==";
    LOG_TRACE(solvestephyst_residual) << "====================================";

//    std::cout << "norm of returnVector start, pre init: " << returnVector.NormL2() << std::endl;
    
		/*
     * 0. Init returnVector
     */
		returnVector.Init();

//    std::cout << "norm of returnVector after init: " << returnVector.NormL2() << std::endl;
    
		/*
     * 1. Setup Linear part and predictor corrector update
     */
    assert(LinRHSRequiresAssembly_ == false);
    assert(PredictorCorrectorRequiresAssembly_ == false);

    LOG_TRACE(solvestephyst_residual) << "=== Add linear rhs and predictor-corrector update to return vector ===";
    returnVector.Add(1.0,linRhsVec_n_0_);
		returnVector.Add(1.0,predictorCorrectorUpdate_n_0_);
    LOG_DBG2(solvestephyst_residual) << "---- Norm of return vector after adding linRHS and predictor-corrector: "<<returnVector.NormL2()<<" ---";
//    std::cout << "norm of returnVector after adding linRhsVec_n_0_ and predictorCorrectorUpdate_n_0_: " << returnVector.NormL2() << std::endl;
		/*
     * 2. Setup nonlin RHS with hysteresis
     * > nonlin RHS is solution dependent
     * > if current solution differes from previously used one >  reassemble
     */
    SBM_Vector diff;
    diff = lastUpdateVecForNonLinRHS_;
    diff.Add(-1.0,solutionForRHSAssembly);
    bool forceHystReassembly = false;
    bool NonLinRHSRequiresReassembly = (diff.NormL2() > 1e-16);

    LOG_DBG2(solvestephyst_residual) << "---- Difference between previously used assembly-solution and requested assembly-solution: " << diff.NormL2() << " ---";

    if(currentRESandRHSFlag_ != flag_for_vector_assembly){
      LOG_DBG(solvestephyst_residual) << vectorFlagToString(currentRESandRHSFlag_) << " to " << vectorFlagToString(flag_for_vector_assembly);
      forceHystReassembly = true;
    }
//    std::cout << " currentRESandRHSFlag_ = " << vectorFlagToString(currentRESandRHSFlag_) << std::endl;
//    std::cout << " flag_for_vector_assembly = " << vectorFlagToString(flag_for_vector_assembly) << std::endl;

    if(forceRHSevaluation_){
      LOG_DBG(solvestephyst_residual) << "--- Force assembly of NL-RHS (only this time)  ---";
    }

    LOG_TRACE(solvestephyst_residual) << "=== Add NL-RHS to return vector ===";
    if( NonLinRHSRequiresReassembly || forceHystReassembly || forceRHSevaluation_ ){
//      std::cout << "Reassemble NL RHS" << std::endl;
//      std::cout << "norm of returnVector at start of reassembly: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_ at start of reassembly: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
      LOG_DBG(solvestephyst_residual) << "--- Reassemble NL-RHS ---";
      LOG_DBG2(solvestephyst_residual) << "---- Norm of previous NL-RHSreturn vector : "<<nonLinRhsVec_n_k_.NormL2()<<" ---";

			nonLinRhsVec_n_k_.Init();

//      std::cout << "norm of returnVector after init: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_ after init: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
			algsys_->InitRHS(nonLinRhsVec_n_k_); //load storage into algsys
//      std::cout << "norm of returnVector after init2: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_ after init2: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
      // evalPurpose 1 > assemble (refers not only for matrix but also for rhs)
      UInt evaluationPurpose = 1;
      PDE_.SetFlagInCoefFncHyst("evaluationPurpose",evaluationPurpose); // set variable directly
      //			PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSHyst",timeLevelRHSHyst);

      Integer timeLevelRHSPol;
      Integer timeLevelRHSStrain;
      Integer timeLevelRHSCoupling;

      Helper_GetHystOperatorParamsFromVectorFlag(flag_for_vector_assembly, timeLevelRHSPol, timeLevelRHSStrain, timeLevelRHSCoupling);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSPol",timeLevelRHSPol);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSStrain",timeLevelRHSStrain);
      PDE_.SetFlagInCoefFncHyst("SetTimeLevelRHSCoupling",timeLevelRHSCoupling);

      int timeLevelMaterial;
      int FPMaterialTensor;
      int FPRHSCorrectionTensor;
      int timeLevelDeltaMatPol;
      int timeLevelDeltaMatStrain;
      int timeLevelDeltaMatCoupling;
      int inversionApproach;
      Helper_GetHystOperatorParamsFromMatrixFlag(flag_for_matrix_assembly, timeLevelMaterial, FPMaterialTensor, FPRHSCorrectionTensor,
              timeLevelDeltaMatPol, timeLevelDeltaMatStrain, timeLevelDeltaMatCoupling, inversionApproach);

      // important: allow reevaluation of hyst operator > might have changed
      PDE_.SetFlagInCoefFncHyst("resetReeval",true);

//      std::cout << "SetupRESorRHS - NONLIN" << std::endl;
      PDE_.SetFlagInCoefFncHyst("SetFPCorrectionFlag",FPRHSCorrectionTensor); // during each rhs assembly

			SBM_Vector solBackup;
      SBM_Vector rhsBackup;
			solBackup = solVec_;
      rhsBackup = rhsVec_;
      
//      std::cout << "norm of rhsBackup before assembleNonLinRHS: " << rhsBackup.NormL2() << std::endl;
//      SBM_Vector zeroVec;
//      zeroVec = rhsBackup;
//      zeroVec.Init();
//      std::cout << "norm of rhsBackup after defining zeroVec: " << rhsBackup.NormL2() << std::endl;
//      std::cout << "norm of zeroVec afterinit: " << zeroVec.NormL2() << std::endl;
      // AssembleNonLinRHS accesses solVec_ for assembly > set solVec_ accordingly
			solVec_ = solutionForRHSAssembly;

      /*
       * DANGEROUS! we set rhs vector to nonlinrhs
       * > if we are not setting the rhs to its correct value before solving, we get wrong results
       * > maybe edit this function, such that it takes a backup of the rhs and then restores it!s
       */
//      algsys_->InitRHS(zeroVec); // test 18.4.2020: clear RHS before assemble nonlinrhs! somehow it seems that
      // the nonlinrhs adds up > no effect
      assemble_->AssembleNonLinRHS();
            
//      std::cout << "norm of returnVector after assembly: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_ after assembly: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
			algsys_->GetRHSVal( nonLinRhsVec_n_k_ );
//            std::cout << "norm of returnVector after load from algsys: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_ after load from algsys: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
      LOG_DBG2(solvestephyst_residual) << "---- Norm of new NL-RHSreturn vector : "<<nonLinRhsVec_n_k_.NormL2()<<" ---";
//      std::cout << "norm of rhsBackup after assembleNonLinRHS: " << rhsBackup.NormL2() << std::endl;
			// reset solVec_ and rhs
			solVec_ = solBackup;
      algsys_->InitRHS(rhsBackup);

      lastUpdateVecForNonLinRHS_ = solutionForRHSAssembly;
      forceRHSevaluation_ = false;
      currentRESandRHSFlag_ = flag_for_vector_assembly;
		} else {
//      std::cout << "DO NOT Reassemble NL RHS" << std::endl;
      LOG_DBG(solvestephyst_residual) << "--- Reuse NL-RHS ---";
    }

    //    std::cout << "TimeLevels_IN_CoefFunctionHyst -- post" << std::endl;
    //    PDE_.SetFlagInCoefFncHyst("PrintTimeLevels",0);
//            std::cout << "norm of returnVector before adding nonLinRhsVec_n_k_: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
		returnVector.Add(1.0,nonLinRhsVec_n_k_);
//                std::cout << "norm of returnVector after adding nonLinRhsVec_n_k_: " << returnVector.NormL2() << std::endl;
//      std::cout << "norm of nonLinRhsVec_n_k_: " << nonLinRhsVec_n_k_.NormL2() << std::endl;
//    
    LOG_DBG2(solvestephyst_residual) << "---- Norm of return vector after adding NL-RHS: "<<returnVector.NormL2()<<" ---";

    /*
     * 3. Assemble system if required
     * > here we have to use the version with prescribed timelevels as those
     *   might differ from the ones that are required by nonLinMethod
     */
    bool reassemblyPerformed = false;
//    std::cout << "flag_for_matrix_assembly = " << matrixFlagToString(flag_for_matrix_assembly) << std::endl;
    if( (flag_for_matrix_assembly == HYSTFREE) ||
            (flag_for_matrix_assembly == HYSTFREE_GLOBALLY_SCALED_B_VERSION) ||
            (flag_for_matrix_assembly == HYSTFREE_LOCALLY_SCALED_B_VERSION) ||
            (flag_for_matrix_assembly == HYSTFREE_LOCALLY_SCALED_H_VERSION) ||
            (flag_for_matrix_assembly == HYSTFREE_SCALED_H_VERSION) ) {

      LOG_DBG(solvestephyst_residual) << "--- Assemble system if this is required; the requested system should be HYSTFREE ---";
      /*
       * Note: AssembleSystem will only trigger a reassembly if this is required (or forced), i.e. if the system is purely
       * linear, it will only be assembled once and otherwise just loaded
       */
//      std::cout << "Check if assembly is needed" << std::endl;
      reassemblyPerformed = AssembleSystem(flag_for_matrix_assembly, solutionForMatrixAssembly,false);//, iterationCounter);

      if(reassemblyPerformed){
//        std::cout << ">>>> system got reassembled" << std::endl;
        LOG_DBG2(solvestephyst_residual) << "---- system was reassembled ----";
      }
    } else {
      EXCEPTION("RESIDUAL/RHS setup has to be performed without DELTAMAT/JACOBIAN activated");
    }
		/*
     * 4. Update res/rhs according to update solution vector
     */
    SBM_Vector systemUpdate;
    Double normOfSystemUpdate;
    systemUpdate = zeroVec_;
    // note: compute KeffTimesSolution does not assemble system;
    // NOTE: solutionForUpdate always has to be a full solution vector, no increment!
    //       however, it need not be the actual solution but can also be the one from the prev. ts, it etc
//    std::cout << "res without Keff*sol = " << returnVector.NormL2() << std::endl;

    if(comingFromLinesearch==false){
//      std::cout << "Coming from LS = false" << std::endl;
      LOG_TRACE(solvestephyst_residual) << "=== Subtract K_eff*solutionVector (std. case) ===";
      Helper_ComputeKeffTimesSolution(solutionForUpdate, systemUpdate,!reassemblyPerformed, flag_for_matrix_assembly);
      normOfSystemUpdate = systemUpdate.NormL2();
      returnVector.Add(-1.0,systemUpdate);
      LOG_DBG2(solvestephyst_residual) << "---- Norm of return vector after adding K_eff*solution: "<<returnVector.NormL2()<<" ---";
    } else {
      if(reassemblyPerformed){
        LOG_TRACE(solvestephyst_residual) << "=== Subtract K_eff*solutionVector (linesearch case, full eval. due to reassembly of system) ===";
//        std::cout << "=== Subtract K_eff*solutionVector (linesearch case, full eval. due to reassembly of system) ===" << std::endl;
        Helper_ComputeKeffTimesSolution(solutionForUpdate, systemUpdate,!reassemblyPerformed, flag_for_matrix_assembly);
        normOfSystemUpdate = systemUpdate.NormL2();
        returnVector.Add(-1.0,systemUpdate);
        LOG_DBG2(solvestephyst_residual) << "---- Norm of return vector after adding K_eff*solution: "<<returnVector.NormL2()<<" ---";
      } else {
        LOG_DBG(solvestephyst_residual) << "--- Return 0 to trigger construction of K_eff*solutionVector from oldSolution and solutionUpdate (linesearch) ---";
//        std::cout << "=== DO NOT Subtract K_eff*solutionVector (linesearch case, full eval. due to reassembly of system) ===" << std::endl;
        normOfSystemUpdate = 0;
      }
    }

//    std::cout << "normOfSystemUpdate (in setupRESorRHS: " << normOfSystemUpdate << std::endl;

    return normOfSystemUpdate;
  }

  /*
   * V. Helper and output functions
   */
  void SolveStepHyst::Helper_GetParamsFromSolveFlag(solveFlag flag_for_systemSolution, matrixFlag& flag_for_systemMatrix, matrixFlag& flag_for_resEvalMatrix,
          vectorFlag& flag_for_rhsVector, vectorFlag& flag_for_resVector, bool& estimatorUseful){
    /*
     * For each solution approach (determined by input flag_for_systemSolution),
     * there are four important solution steps that require interaction with algsys and the hyst operators
     * a) assemble system that shall be solved
     * b) setup rhs that shall be solved for
     * c) assemble system for residual evaluation
     * d) evaluate hysteretic rhs vector for residual evaluation
     * > this function shall return the corresponding flags that can than be converted further to input parameter
     *    that are passed to hyst operator / algsys
     * > example:
     *    flag_for_systemSolution = SOLVE_JACOBIAN
     *    a) system shall be FD-approximation of Jacobian of residual > flag_for_systemMatrix = FDJACOBIAN
     *    b) rhs has to be the current residual > flag_for_rhsVector = CURRENT_RES
     *    c) to setup the residual, we require the hyst-free system matrix > flag_for_resEvalMatrix = HYSTFREE
     *    d) residual = current residual > flag_for_resVector = CURRENT_RES
     *
     * bool flag estimatorUseful:
     *  in case of quasi-Newton methods it might be useful to
     *  apply a fixpoint estimator step to determine if we have reached a reversal point;
     *  for fixpoint cases, it makes no sense, though
     */

    switch(flag_for_systemSolution){
      case SOLVE_NOTSET :
        EXCEPTION("No valid solution approach was set; cannot determine an appropriate set of flags");
      case SOLVE_GLOBAL_FIXPOINT_B :
        /*
         * Global Fixpoint B
         * > solve (non-)linear, hystFree system with scaled matrerial tensors (e.g. mu_FP = 2*mu_0)
         * > solve for current residual vector
         */
        flag_for_systemMatrix = HYSTFREE_GLOBALLY_SCALED_B_VERSION;
        flag_for_resEvalMatrix = HYSTFREE; //HYSTFREE_GLOBALLY_SCALED_B_VERSION;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = false;
        break;
      case SOLVE_LOCAL_FIXPOINT_B :
        /*
         * Local Fixpoint B
         * > solve (non-)linear, hystFree system with scaled matrerial tensors (e.g. mu_FP = 2*mu_0)
         * > solve for current residual vector
         */
        flag_for_systemMatrix = HYSTFREE_LOCALLY_SCALED_B_VERSION;
        flag_for_resEvalMatrix = HYSTFREE_LOCALLY_SCALED_B_VERSION;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = false;
        break;
      case SOLVE_LOCAL_FIXPOINT_B_v2 :
        /*
         * Local Fixpoint B v2
         * > solve (non-)linear, hystFree system with scaled matrerial tensors (e.g. mu_FP = 2*mu_0)
         * > solve for current residual vector
         */
        flag_for_systemMatrix = HYSTFREE_LOCALLY_SCALED_B_VERSION;
        flag_for_resEvalMatrix = HYSTFREE; // use standard un-scaled matrix for residual and rhs computation
        // by this, we do not need to correct the magnetization term
        // Proof:
        //  rot H = rot nu_0 B - rot M
        //        = rot nu_FP B - rot nu_FP B + rot nu_0 B - rot M
        //        = rot nu_FP B - rot M + rot (nu_0 - nu_FP) B = J
        // > rot nu_FP B = J + rot M - rot (nu_0 - nu_FP) B
        // > rot nu_FP B_k+1 = J + rot M(B_k) - rot (nu_0 - nu_FP) B_k  << eigentlich m��sste hier B_k+1 stehen!
        // > rot nu_FP deltaB = rot nu_FP (B_k+1 - B_k)
        //       = J + rot M(B_k) - rot (nu_0 - nu_FP) B_k - rot nu_FP B_k
        //       = J + rot M(B_k) - rot nu_0 B_k + nu_FP B_k - rot nu_FP B_k
        //       = J + rot M(B_k) - rot nu_0 B_k = std residual with Hystfree system matrix update (= rot nu_0 B_k)
        //  > auf RHS m��sste eigentlich der term rot (nu_0 - nu_FP) (B_k+1 - B_k) stehen
        //  > geht aber nicht, da B_k+1 nicht bekannt > approximation mittels B_k wie oben ist wohl die beste Option

        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = false;
        break;
      case SOLVE_GLOBAL_FIXPOINT_H :
        /*
         * Fixpoint H
         * > solve (non-)linear, hystFree system with scaled matrerial tensors (different scaling as in B-form)
         * > solve for current residual vector
         * > no local inversion of hyst operator; take last known value of H for evaluation
         *
         * Pseudo-Algorithm as indicated in
         * "Anwendung der Methode der finiten Elemente und des Vektor-Preisach-Modells zur Berechnung ebener magnetostatischer Felder in hysteresebehafteten Medien" - Kurz et al.
         *
         * 1) Use current value H_k to evaluate FORWARD preisach model (resulting in P(H_k) where P shall denote the magnetic polarization)
         * 2) Solve Fixpoint system (based on B = \mu_FP H + P_FP(H) = \mu_0 H + ( P_FP(H) + (\mu_FP - \mu_0) H) = \mu_0 H + P(H),
         *        i.e. P_FP(H) = P(H) + (\mu_FP - \mu_0) H; \nu_FP = \mu_FP^-1 )
         *      rot(\nu_FP B_k+1) = J + rot(\nu_FP P(H_k))
         * 3) Estimate H^* that is consistent with B_k+1 and P(H_k), i.e.
         *      H^* = \nu_FP (B_k+1 - P(H_k))
         * 4) Compute new H_k+1 using either linesearch or a fitting constant relaxation factor omega
         *      H_k+1 = (1-omega)*H_k + omega*H^*
         * 5) Continue at 1 till convergence
         *
         * Note:
         * a) Kurz et al. suggest the secondary souce approach, i.e. H is constructed from H_j and H_m where
         *    H_j is the solution of rot(H_j) = J
         *    and
         *    H_m is the solution of rot(H_m) = rot(M) (M being the magnetization)
         * b) For the H-Version, \mu_FP can be set to \mu_0 directly; the actual important part to ensure convergence
         *    is the choice of omega which should be between 0 and 2/\mu_r_max
         *
         * > by setting HYSTFREE_SCALED_H_VERSION we set the FP material tensor to \mu_0
         * > the first time we set this flag, we trace the hysteresis curve to retrieve \mu_r_max
         * > the flag useFPH gets set in CoefFunctionHyst which disables the inversion of the Preisach operator and
         *   insteads computes H as indicated by steps 3 and 4 above
         *
         */
        flag_for_systemMatrix = HYSTFREE_SCALED_H_VERSION;
        flag_for_resEvalMatrix = HYSTFREE; //HYSTFREE_SCALED_H_VERSION;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = false;
        break;
      case SOLVE_LOCAL_FIXPOINT_H :
        /*
         * Fixpoint H
         * > solve (non-)linear, hystFree system with scaled matrerial tensors (different scaling as in B-form)
         * > solve for current residual vector
         * > no local inversion of hyst operator; take last known value of H for evaluation
         *
         * Pseudo-Algorithm as indicated in
         * "Anwendung der Methode der finiten Elemente und des Vektor-Preisach-Modells zur Berechnung ebener magnetostatischer Felder in hysteresebehafteten Medien" - Kurz et al.
         *
         * 1) Use current value H_k to evaluate FORWARD preisach model (resulting in P(H_k) where P shall denote the magnetic polarization)
         * 2) Solve Fixpoint system (based on B = \mu_FP H + P_FP(H) = \mu_0 H + ( P_FP(H) + (\mu_FP - \mu_0) H) = \mu_0 H + P(H),
         *        i.e. P_FP(H) = P(H) + (\mu_FP - \mu_0) H; \nu_FP = \mu_FP^-1 )
         *      rot(\nu_FP B_k+1) = J + rot(\nu_FP P(H_k))
         * 3) Estimate H^* that is consistent with B_k+1 and P(H_k), i.e.
         *      H^* = \nu_FP (B_k+1 - P(H_k))
         * 4) Compute new H_k+1 using either linesearch or a fitting constant relaxation factor omega
         *      H_k+1 = (1-omega)*H_k + omega*H^*
         * 5) Continue at 1 till convergence
         *
         * Note:
         * a) Kurz et al. suggest the secondary souce approach, i.e. H is constructed from H_j and H_m where
         *    H_j is the solution of rot(H_j) = J
         *    and
         *    H_m is the solution of rot(H_m) = rot(M) (M being the magnetization)
         * b) For the H-Version, \mu_FP can be set to \mu_0 directly; the actual important part to ensure convergence
         *    is the choice of omega which should be between 0 and 2/\mu_r_max
         *
         * > by setting HYSTFREE_SCALED_H_VERSION we set the FP material tensor to \mu_0
         * > the first time we set this flag, we trace the hysteresis curve to retrieve \mu_r_max
         * > the flag useFPH gets set in CoefFunctionHyst which disables the inversion of the Preisach operator and
         *   insteads computes H as indicated by steps 3 and 4 above
         *
         */
        flag_for_systemMatrix = HYSTFREE_LOCALLY_SCALED_H_VERSION;
        flag_for_resEvalMatrix = HYSTFREE; //HYSTFREE_SCALED_H_VERSION;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = false;
        break;
      case SOLVE_JACOBIAN:
      case SOLVE_CHORD:
        /*
         * quasi-Newton schemes > system matrix approximates Jacobian of the residual
         *                       rhs = current residual
         * SOLVE_JACOBIAN
         * > use Finite Differences on element level to approximate element wise residual;
         * > for residual evaluation we use the hysteresis-free system matrix and add the contribution of the
         *       hysteresis operator via nonlinear rhs
         *
         * SOLVE_CHORD
         * > as SOLVE_JACOBIAN but Jacobian gets evaluated only once (or if forced)
         */
        flag_for_systemMatrix = FDJACOBIAN;
        flag_for_resEvalMatrix = HYSTFREE;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = true;
        break;
      case SOLVE_BGM:
      case SOLVE_BGM_ENHANCED:
        /*
         * quasi-Newton schemes > system matrix approximates Jacobian of the residual
         *                       rhs = current residual
         * SOLVE_BGM
         * > estimate Jacobian / its contribution to the solution from the previous residual and solution
         *  increments on global level
         * > to handle boundary conditions correctly and to get a fitting starting point, solve initial system using
         *  FDJacobian
         * > for residual evaluation we use the hysteresis-free system matrix and add the contribution of the
         *  hysteresis operator via nonlinear rhs
         */
        flag_for_systemMatrix = FDJACOBIAN;
        flag_for_resEvalMatrix = HYSTFREE;
        flag_for_rhsVector = CURRENT_RES;
        flag_for_resVector = CURRENT_RES;
        // note: BGM does not use the estimator during its actual iterations as we compute the global Jacobian
        //        from solution and residual increments instead of using the Hyst operator
        //        however, during the first iteration of BGM or after BGM stalls out, we perform a standard solution
        //        step via algsys using a FDJacobian (which profits from an estimator)
        estimatorUseful = true;
        break;
      case SOLVE_DELTAMAT_IT:
        /*
         * quasi-Newton schemes > system matrix approximates Jacobian of the residual
         *                       rhs = current residual
         * SOLVE_DELTAMAT_IT
         * > use delta matrix approach on element level to approximate element wise residual; delta matrix computed
         *    from current hysteresis values and values from previous iteration
         * > for residual evaluation we use the hysteresis-free system matrix and add the contribution of the
         *    hysteresis operator via nonlinear rhs
         * > on rhs we actually should solve for the residual towards the last iteration but tests showed that the
         *  current residual does better
         */
        flag_for_systemMatrix = DELTAMATJACOBIAN_LAST_IT;
        flag_for_resEvalMatrix = HYSTFREE;
        flag_for_rhsVector = CURRENT_RES; // instead of PREVIOUS_IT_RES
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = true;
        break;
      case SOLVE_DELTAMAT_TS:
        /*
         * quasi-Newton schemes > system matrix approximates Jacobian of the residual
         *                       rhs = current residual
         * SOLVE_DELTAMAT_TS
         * > use delta matrix approach on element level to approximate element wise residual; delta matrix computed
         *    from current hysteresis values and values from previous time step
         * > for residual evaluation we use the hysteresis-free system matrix and add the contribution of the
         *    hysteresis operator via nonlinear rhs
         * > on rhs we solve for the residual towards the last time step
         */
        flag_for_systemMatrix = DELTAMATJACOBIAN_LAST_TS;
        flag_for_resEvalMatrix = HYSTFREE;
        flag_for_rhsVector = PREVIOUS_TS_RES;
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = true;
        break;
      case SOLVE_DELTAMAT_TS_TOWARDS_IT:
        /*
         * quasi-Newton schemes > system matrix approximates Jacobian of the residual
         *                       rhs = current residual
         * SOLVE_DELTAMAT_TS
         * > use delta matrix approach on element level to approximate element wise residual; delta matrix computed
         *    from current hysteresis values and values from previous time step
         * > for residual evaluation we use the hysteresis-free system matrix and add the contribution of the
         *    hysteresis operator via nonlinear rhs
         * > on rhs we actually should solve for the residual towards the last iteration but tests showed that the
         *  current residual does better
         * > uses deltaMatrix as SOLVE_DELTAMAT_TS but steps towards last iterate as SOLVE_DELTAMAT_IT
         */
        flag_for_systemMatrix = DELTAMATJACOBIAN_LAST_TS;
        flag_for_resEvalMatrix = HYSTFREE;
        flag_for_rhsVector = CURRENT_RES; // instead of PREVIOUS_TS_RES
        flag_for_resVector = CURRENT_RES;
        estimatorUseful = true;
        break;
      default:
        EXCEPTION("Solution approach "<<solveFlagToString(SOLVE_NOTSET)<<" unknown");
    }
  }

  void SolveStepHyst::Helper_GetHystOperatorParamsFromMatrixFlag(matrixFlag flag_for_assembly, int& timeLevelMaterial, int& FPMaterialTensor,
          int& FPRHSCorrectionTensor, int& timeLevelDeltaMatPol, int& timeLevelDeltaMatStrain, int& timeLevelDeltaMatCoupling, int& inversionApproach){

    int timeLevelDeltaMatGeneral;
    switch(flag_for_assembly){
      case NOTSET:
        EXCEPTION("No valid flag for matrix assembly was set; cannot determine an appropriate set of time levels for hyst operator");
      case HYSTFREE:
        /*
         * (non-)linear but hysteresis free system (i.e. no approximation of Jacobian!)
         * > deltaMatrix is deactivated completely (-2)
         * > material tensors are evaluated at the current time step though (in case on non-linearities)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 0; // no additional scaling as in fixpoint case
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -2;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case FDJACOBIAN:
        /*
         * (non-)linear approximation of Jacobian build from Finite Differences on element level
         * > deltaMatrix is set to 3
         * > material tensors are evaluated at the current time step (in case on non-linearities)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 0; // no additional scaling as in fixpoint case
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = 3;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case DELTAMATJACOBIAN_LAST_TS:
        /*
         * (non-)linear approximation of Jacobian on element level using deltaMatrix approach
         * (i.e. compute approximation from current and previous hyst values; here: previous value = value from last time step)
         * > deltaMatrix is set to -1 to trigger usage of last time step values
         * > material tensors are evaluated at the current time step (in case on non-linearities)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 0; // no additional scaling as in fixpoint case
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -1;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case DELTAMATJACOBIAN_LAST_IT:
        /*
         * (non-)linear approximation of Jacobian on element level using deltaMatrix approach
         * (i.e. compute approximation from current and previous hyst values; here: previous value = value from last iteration)
         * > deltaMatrix is set to 1 to trigger usage of last iteration values
         * > material tensors are evaluated at the current time step (in case on non-linearities)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 0; // no additional scaling as in fixpoint case
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = 1;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case REUSED:
        /*
         * (non-)linear approximation of Jacobian on element level; reuse last computed Jacobian/deltaMatrix
         * > deltaMatrix is set to 2 to trigger reusebuild from Finite Differences on element level
         * > material tensors are evaluated at the current time step (in case on non-linearities)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 0; // no additional scaling as in fixpoint case
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = 2;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case HYSTFREE_GLOBALLY_SCALED_B_VERSION:
        /*
         * (non-)linear but hysteresis free system (i.e. no approximation of Jacobian!)
         * > deltaMatrix is deactivated completely (-2)
         * > material tensors are evaluated at the current time step though (in case on non-linearities)
         * > material tensors are scaled for better convergence (B-version case)
         * > scaling factor is constant for all timesteps (i.e. global)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 1;
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -2;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case HYSTFREE_LOCALLY_SCALED_B_VERSION:
        /*
         * (non-)linear but hysteresis free system (i.e. no approximation of Jacobian!)
         * > deltaMatrix is deactivated completely (-2)
         * > material tensors are evaluated at the current time step though (in case on non-linearities)
         * > material tensors are scaled for better convergence (B-version case)
         * > scaling factor is determined once every timestep from current slope (i.e. local)
         *  > quite similar to chord method
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 2;
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -2;
        inversionApproach = 0; // local inversion (default; use inversion as defined in mat.xml)
        break;
      case HYSTFREE_SCALED_H_VERSION:
        /*
         * (non-)linear but hysteresis free system (i.e. no approximation of Jacobian!)
         * > deltaMatrix is deactivated completely (-2)
         * > material tensors are evaluated at the current time step though (in case on non-linearities)
         * > material tensors are scaled for better convergence (H-version case)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 3;
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -2;
        inversionApproach = 1; // no local inversion (evaluate hyst operator using the latest value of H)
        break;
      case HYSTFREE_LOCALLY_SCALED_H_VERSION:
        /*
         * (non-)linear but hysteresis free system (i.e. no approximation of Jacobian!)
         * > deltaMatrix is deactivated completely (-2)
         * > material tensors are evaluated at the current time step though (in case on non-linearities)
         * > material tensors are scaled for better convergence (H-version case)
         */
        timeLevelMaterial = 0; // current value
        FPMaterialTensor = 4;
        FPRHSCorrectionTensor = FPMaterialTensor;
        timeLevelDeltaMatGeneral = -2;
        inversionApproach = 1; // no local inversion (evaluate hyst operator using the latest value of H)
        break;
      default:
        EXCEPTION("Matrix flag "<<matrixFlagToString(flag_for_assembly,0)<<" unknown");
    }

    timeLevelDeltaMatPol = timeLevelDeltaMatGeneral;
    // delta matrix formulation is selected in input.xml
    // a) Jacobian approximates dPol/dE or dPol/dB only
    // b) Jacobian approximates dPol/dE and dS/dE (accordingly for B)
    // c) Jacobian approximates dPol/dE, dS/dE and de/dE (accordingly for B)
    // other combinations would be possible, too (but a,b,c are the ones which make most sense)
    if(deltaMatStrain_){
      timeLevelDeltaMatStrain = timeLevelDeltaMatGeneral;
    } else {
      // deactivate specific deltaMatrix only
      timeLevelDeltaMatStrain = 0;
    }

    if(deltaMatCoupling_){
      timeLevelDeltaMatCoupling = timeLevelDeltaMatGeneral;
    } else {
      // deactivate specific deltaMatrix only
      timeLevelDeltaMatCoupling = 0;
    }
  }

  void SolveStepHyst::Helper_GetHystOperatorParamsFromVectorFlag(vectorFlag flag_for_res_or_rhs, int& timeLevelPol, int& timeLevelSrain,
          int& timeLevelCoupling){

    int timeLevelHystGeneral;
    switch(flag_for_res_or_rhs){
      case NOTSET_VEC:
        EXCEPTION("No valid flag for rhs/res assembly was set; cannot determine an appropriate set of time levels for hyst operator");
      case PREVIOUS_TS_RES:
        /*
         * setup full residual, i.e.
         * - linear rhs
         * - non-linear rhs
         * - Keff*solution (Keff must not be a deltaMatrix!)
         * - hystOperator from previous time step
         */
        timeLevelHystGeneral = -1;
        break;
      case CURRENT_RES:
        /*
         * setup full residual, i.e.
         * - linear rhs
         * - non-linear rhs
         * - Keff*solution (Keff must not be a deltaMatrix!)
         * - current state of hyst operator
         */
        timeLevelHystGeneral = 0;
        break;
      case PREVIOUS_IT_RES:
        /*
         * setup full residual, i.e.
         * - linear rhs
         * - non-linear rhs
         * - Keff*solution (Keff must not be a deltaMatrix!)
         * - hystOperator from previous iteration
         */
        timeLevelHystGeneral = 1;
        break;
      case RES_REM:
        /*
         * setup full residual, i.e.
         * - linear rhs
         * - non-linear rhs
         * - Keff*0 (Keff must not be a deltaMatrix!)
         * - hystOperator evaluated at 0 (remanence; actually not correct for magnetics as we would require coercive field there)
         */
        timeLevelHystGeneral = 0; // do not forget to set solVec to 0 and then reevaluate!
        break;
      default:
        EXCEPTION("Vector flag "<<vectorFlagToString(flag_for_res_or_rhs,0)<<" unknown");
    }

    // if delta formulation is active for strains/coupling we consider have to set the
    // hyst values for strains/coupling tensors to the same timelevel as the one for polarization
    // i.e. deltaPol towards last ts > deltaStrain towards last ts
    // if delta formulation is not active, set timelevel to 0 (to get current value)
    timeLevelPol = timeLevelHystGeneral;
    if(deltaMatStrain_){
      timeLevelSrain = timeLevelHystGeneral;
    } else {
      // deactivate specific deltaMatrix only
      timeLevelSrain = 0;
    }

    if(deltaMatCoupling_){
      timeLevelCoupling = timeLevelHystGeneral;
    } else {
      // deactivate specific deltaMatrix only
      timeLevelCoupling = 0;
    }
  }

  void SolveStepHyst::Helper_ComputeKeffTimesSolution(SBM_Vector& solutionForUpdate, SBM_Vector& returnVector,bool useBackup, matrixFlag flag_for_matrix_assembly){
    LOG_TRACE(solvestephyst_helperfuncs) << "== Helper_ComputeKeffTimesSolution (CalcKu) ==";

    if(useBackup == false){
//      std::cout << "Compute KeffTimesSolution; no backup; use current system matrix" << std::endl;
      LOG_DBG(solvestephyst_helperfuncs) << "(CalcKu) --- UseBackup == false ---";

      LOG_DBG(solvestephyst_helperfuncs) << "(CalcKu) --- 1. rebuild system ---";
      BuildSystem(true);
//      std::cout << "Current build flag: " << matrixFlagToString(currentBuildFlag_) << std::endl;

      LOG_DBG(solvestephyst_helperfuncs) << "(CalcKu) --- 2. multiply with newly built system ---";
      algsys_->ComputeSysMatTimesVector(SYSTEM,solutionForUpdate,returnVector,false);
    } else {
//      std::cout << "use backuped system for KeffTimesSolution" << std::endl;
      LOG_DBG(solvestephyst_helperfuncs) << "(CalcKu) --- UseBackup == true ---";

      LOG_DBG(solvestephyst_helperfuncs) << "(CalcKu) --- multiply with backup of HYSTFREE system ---";
      if(flag_for_matrix_assembly == HYSTFREE){
//        std::cout << "Compute KeffTimesSolution; load backup; use HYSTFREE matrix" << std::endl;
        algsys_->ComputeSysMatTimesVector(SYSTEM_HYSTFREE,solutionForUpdate,returnVector,false);
      } else {
//        std::cout << "Compute KeffTimesSolution; load backup; use SYSTEM_FIXPOINT matrix" << std::endl;
        algsys_->ComputeSysMatTimesVector(SYSTEM_FIXPOINT,solutionForUpdate,returnVector,false);
      }

    }
    LOG_DBG2(solvestephyst_helperfuncs) << "(CalcKu) --- NORM OF INPUT: " << solutionForUpdate.NormL2();
    LOG_DBG2(solvestephyst_helperfuncs) << "(CalcKu) --- NORM OF OUTPUT: " << returnVector.NormL2();
  }

  // TODO: make use of these functions
  void SolveStepHyst::Helper_CreateSystemBackup(bool backupSOL, bool backupRHS, bool backupHystValues, bool backupSYSTEM, FEMatrixType backupSLOT){
    if(backupSOL){
      // solVec_ directly points to system!
      // i.e. changing solVec_ will influence algsys; same is true for rhsVec_
      solVecBAK_ = solVec_;
    }
    if(backupRHS){
      rhsVecBAK_ = rhsVec_;
    }
    if(backupHystValues){
      PDE_.SetFlagInCoefFncHyst("SetPreviousHystValues",3);
    }
    if(backupSYSTEM){
      algsys_->BackupCurrentSystemMatrix(backupSLOT);
    }
  }

  void SolveStepHyst::Helper_LoadSystemBackup(bool loadSOL, bool loadRHS, bool loadHystValues, bool loadSYSTEM, FEMatrixType loadSLOT){
    if(loadSOL){
      solVec_ = solVecBAK_;
    }
    if(loadRHS){
      rhsVec_ = rhsVecBAK_;
    }
    if(loadHystValues){
      PDE_.SetFlagInCoefFncHyst("RestorePreviousHystValues",3);
    }
    if(loadSYSTEM){
      algsys_->LoadBackupToCurrentSystemMatrix(loadSLOT);
    }
  }

  Vector<UInt> SolveStepHyst::Helper_GetStorageIndices(){
    /*
     * Retrive storage indices for residual/increment updates
     */
    /*
     * for L-BFGS and L-BGM we hold up to hisotryLengthMax_ residual/incrmeent updates
     * in storage;
     *
     * Convention:
     *  returnVec[0] = index/position of oldest increment
     * if(historyLength_ >= historyLengthMax_)
     *  returnVec[historyLengthMax_-1] = index of latest increment
     * else
     *  returnVec[historyLength_-1] = index of latest increment
     */
    Vector<UInt> indices;
//    std::cout << "historyLength_: " << historyLength_ << std::endl;
//    std::cout << "hisotryLengthMax_: " << hisotryLengthMax_ << std::endl;
    if(historyLength_ <= hisotryLengthMax_){
      /*
       * unless the maximal number of stored increments is reached,
       * storage index and position in array coincide
       */
      indices = Vector<UInt>(historyLength_);
      for(UInt i = 0; i < historyLength_; i++){
        indices[i] = i;
      }

    } else {
      /*
       * if maximal size is exceeded, wind around and overwrite from
       * 0 on
       */
      indices = Vector<UInt>(hisotryLengthMax_);
      for(UInt i = 0; i < hisotryLengthMax_; i++){
        UInt position = (historyLength_ + i)%hisotryLengthMax_;
        indices[i] = position;
      }
    }
//    std::cout << "Indices: " << indices.ToString() << std::endl;
    return indices;
  }

  // TODO: eventually move to global SBM class
  void SolveStepHyst::Helper_MatrixTimesSBMVector(Matrix<Double>& matIn, SBM_Vector& vecIn, SBM_Vector& vecOut){

    UInt numRows, numCols;
    numRows = matIn.GetNumRows();
    numCols = matIn.GetNumCols();

    //    std::cout << "Values retrieved from matIn: " << std::endl;
    //    std::cout << "numRows x numCols: " << numRows << " x " << numCols << std::endl;
    //
    UInt size1 = 0;
    UInt size2 = 0;
    UInt numVecs1 = 0;
    UInt numVecs2 = 0;

    Vector<UInt> sizeSubvectors1 = Vector<UInt>(1);
    Vector<UInt> sizeSubvectors2 = Vector<UInt>(1);

    Vector<UInt> posOffsets1 = Vector<UInt>(1);
    Vector<UInt> posOffsets2 = Vector<UInt>(1);

    Helper_GetSizeInfoFromSBMVector(vecIn, size1, numVecs1, sizeSubvectors1, posOffsets1);

    //    std::cout << "Values retrieved from vecIn: " << std::endl;
    //    std::cout << "totalSize: " << size1 << std::endl;
    //    std::cout << "numSubvectors: " << numVecs1 << std::endl;
    //    std::cout << "sizeSubvectors: " << sizeSubvectors1 << std::endl;
    //    std::cout << "positionOffsets: " << posOffsets1 << std::endl;
    //
    Helper_GetSizeInfoFromSBMVector(vecOut, size2, numVecs2, sizeSubvectors2, posOffsets2);

    //    std::cout << "Values retrieved from vecOut: " << std::endl;
    //    std::cout << "totalSize: " << size2 << std::endl;
    //    std::cout << "numSubvectors: " << numVecs2 << std::endl;
    //    std::cout << "sizeSubvectors: " << sizeSubvectors2 << std::endl;
    //    std::cout << "positionOffsets: " << posOffsets2 << std::endl;
    //
    bool checkSize = true;

    if(checkSize){
      if(numCols != size1){
        EXCEPTION("Matrix numCols != VecIn total length");
      }
      if(numRows != size2){
        EXCEPTION("Matrix numRows != VecOut total length");
      }
    }

    Matrix<Double> tmpMat = Matrix<Double>(1,1);
    Vector<Double> tmpVec = Vector<Double>(1);

    vecOut.Init();

    for(UInt jj = 0; jj < numVecs2; jj++){
      tmpVec.Resize(sizeSubvectors2[jj]);

      for(UInt ii = 0; ii < numVecs1; ii++){
        tmpMat.Resize(sizeSubvectors2[jj],sizeSubvectors1[ii]);
        matIn.GetSubMatrix(tmpMat,posOffsets2[jj],posOffsets1[ii]);

        matIn.Mult(vecIn(ii),tmpVec);

        vecOut(jj).Add(tmpVec);

        //        std::cout << "ii / jj: " << ii << " / " << jj << std::endl;
        //        std::cout << "vecOut(jj): " << vecOut(jj).ToString() << std::endl;
      }
    }
  }

  void SolveStepHyst::Helper_GetSizeInfoFromSBMVector(SBM_Vector& vecIn, UInt& totalSize,
          UInt& numSubvectors, Vector<UInt>& sizeSubvectors, Vector<UInt>& posOffsets ){

    totalSize = 0;
    numSubvectors = vecIn.GetSize();
    sizeSubvectors.Resize(numSubvectors);
    posOffsets.Resize(numSubvectors);
    UInt curSize = 0;

    for(UInt i = 0; i < numSubvectors; i++){
      curSize = vecIn.GetPointer(i)->GetSize();
      totalSize += curSize;
      sizeSubvectors[i] = curSize;
      if(i == 0){
        posOffsets[i] = 0;
      } else {
        posOffsets[i] = sizeSubvectors[i-1]+posOffsets[i-1];
      }
    }
  }

  void SolveStepHyst::Helper_DyadicProductFromSBMVectors(Matrix<Double>& dyadicProduct, SBM_Vector& vec1,
          SBM_Vector& vec2,Double scaling,bool addIdentity){

    UInt size1 = 0;
    UInt size2 = 0;
    UInt numVecs1 = 0;
    UInt numVecs2 = 0;

    Vector<UInt> sizeSubvectors1 = Vector<UInt>(1);
    Vector<UInt> sizeSubvectors2 = Vector<UInt>(1);

    Vector<UInt> posOffsets1 = Vector<UInt>(1);
    Vector<UInt> posOffsets2 = Vector<UInt>(1);

    Helper_GetSizeInfoFromSBMVector(vec1, size1, numVecs1, sizeSubvectors1, posOffsets1);
    Helper_GetSizeInfoFromSBMVector(vec2, size2, numVecs2, sizeSubvectors2, posOffsets2);

    dyadicProduct.Resize(size1,size2);
    Matrix<Double> tmp = Matrix<Double>(1,1);

    for(UInt ii = 0; ii < numVecs1; ii++){

      for(UInt jj = 0; jj < numVecs2; jj++){

        tmp.Resize(sizeSubvectors1[ii],sizeSubvectors2[jj]);
        tmp.DyadicMult(vec1(ii),vec2(jj));

        dyadicProduct.SetSubMatrix(tmp,posOffsets1[ii],posOffsets2[jj]);
      }
    }

    if(scaling != 1.0){
      for(UInt i = 0; i < size1; i++){
        for(UInt j = 0; j < size2; j++) {
          dyadicProduct(i,j) *= scaling;
        }
      }
    }

    if(addIdentity){
      for(UInt i = 0; i < size1; i++){
        for(UInt j = 0; j < size2; j++) {
          if(i == j){
            dyadicProduct(i,j) += 1.0;
          }
        }
      }
    }
  }

  // TODO: clean up arguments
  void SolveStepHyst::Output_WriteHystIterToInfoXML(const std::string& pdeName,
          const std::string& nonLinSolveMethod,
          const std::string& linesearchMethod,
          const UInt coupledIterStep,
          const UInt solStep,
          const UInt iterationCounter,
          const std::string& convergenceLogMSG){

//          const UInt iterationCounter,
//          const UInt iterationCounterTotal,
//          const Double residualErrAbs, const Double residualErrRel, bool useRelativeRelErr,
//          const Double incrementalErrAbs, const Double incrementalErrRel, bool useRelativeIncErr,
//          double etaLineSearch, bool reset, bool failback){

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
    PtrParamNode info = iteration->GetByVal("it","value",iterationCounter,ParamNode::INSERT)->Get("info",ParamNode::APPEND);
//    if(reset){
//      info->Get("reset",ParamNode::APPEND);
//    } else if(failback){
//      info->Get("failback",ParamNode::APPEND);
//    } else {
      info->Get("method")->SetValue(nonLinSolveMethod);
      info->Get("linesearch")->SetValue(linesearchMethod);
      info->Get("convergence")->SetValue(convergenceLogMSG);
//
//      if(linesearchMethod != "none"){
//        info->Get("etaLS")->SetValue(etaLineSearch);
//      }
//      PtrParamNode increment = iteration->GetByVal("it","value",iterationCounterTotal,ParamNode::INSERT)->Get("increment",ParamNode::APPEND);
//      increment->Get("incErr")->SetValue(incrementalErrAbs);
//      increment->Get("incErrRel")->SetValue(incrementalErrRel);
//      increment->Get("incErrCriterion")->SetValue(incStopCrit_);
//      if(useRelativeIncErr){
//        increment->Get("useRelativeErr")->SetValue("True");
//      } else {
//        increment->Get("useRelativeErr")->SetValue("False");
//      }
//      PtrParamNode residual = iteration->GetByVal("it","value",iterationCounterTotal,ParamNode::INSERT)->Get("residual",ParamNode::APPEND);
//      residual->Get("resErr")->SetValue(residualErrAbs);
//      residual->Get("resErrRel")->SetValue(residualErrRel);
//      residual->Get("resErrCriterion")->SetValue(residualStopCrit_);
//      if(useRelativeIncErr){
//        residual->Get("useRelativeErr")->SetValue("True");
//      } else {
//        residual->Get("useRelativeErr")->SetValue("False");
//      }
//    }
  }

  std::string SolveStepHyst::linesearchFlagToString(linesearchFlag flag){
    std::stringstream os;
    switch(flag)
    {
      case NO_LINESEARCH : os << " no linesearch (FULL STEP)"; break;
      case LS_GOLDENSECTION : os << "GoldenSection-Linesearch (EXACT)"; break;
      case LS_QUADINTERPOL : os << "QuadraticInterpolation-Linesearch (EXACT)"; break;
      case LS_VECBASEDPOLY : os << "VectorBasedPolynomial-Linesearch (SUFF. DECREASE)"; break;
      case LS_THREEPOINTPOLY : os << "ThreePointPolynomial-Linesearch (SUFF. DECREASE)"; break;
      case LS_BACKTRACKING_ARMIJO : os << "Backtracking-Linesearch (SUFF. DECREASE, classical Armijo-rule)"; break;
      case LS_BACKTRACKING_ARMIJO_NONMONOTONIC : os << "Backtracking-Linesearch (SUFF. DECREASE, non monotonic Armijo-rule)"; break;
      case LS_BACKTRACKING_GRADFREE : os << "Backtracking-Linesearch (SUFF. DECREASE, gradientFree Armijo-rule)"; break;
      case LS_TRIAL_AND_ERROR : os << "Backtracking-Linesearch (BEST FOUND DECREASE, trial-and-error)"; break;
      default : os.setstate(std::ios_base::failbit);
    }

    return os.str();
  }

  std::string SolveStepHyst::solveFlagToString(solveFlag flag){
    std::stringstream os;
    switch(flag)
    {
      case SOLVE_NOTSET : os << "NO solution method selected"; break;
      case SOLVE_GLOBAL_FIXPOINT_B : os << "Fixpoint (B-stepping); global slope"; break;
      case SOLVE_LOCAL_FIXPOINT_B : os << "Fixpoint (B-stepping); local slope"; break;
      case SOLVE_LOCAL_FIXPOINT_B_v2 : os << "Fixpoint (B-stepping); local slope, v2"; break;
      case SOLVE_GLOBAL_FIXPOINT_H : os << "Fixpoint (H-stepping)"; break;
      case SOLVE_DELTAMAT_IT : os << "DELTA_MAT_IT, delta to last IT"; break;
      case SOLVE_DELTAMAT_TS : os << "DELTA_MAT_TS, delta to last TS"; break;
      case SOLVE_DELTAMAT_TS_TOWARDS_IT : os << "DELTA_MAT_TS/IT, delta to last TS, incr. to last IT"; break;
      case SOLVE_JACOBIAN : os << "Quasi-Newton (FD-Jacobian)"; break;
      case SOLVE_CHORD : os << "Quasi-Newton (FD-Jacobian, chord method, evaluated only at start of TS)"; break;
      case SOLVE_BGM : os << "Quasi-Newton (Broyden's method)"; break;
      case SOLVE_BGM_ENHANCED : os << "Quasi-Newton (FD-Jacobian + Broyden's method)"; break;
      default : os.setstate(std::ios_base::failbit);
    }

    return os.str();
  }

  std::string SolveStepHyst::matrixFlagToString(matrixFlag flag, int verbose){
    std::stringstream os;
    if(verbose == 0){
      switch(flag)
      {
        case NOTSET   : os << "NOTSET";    break;
        case HYSTFREE: os << "HYSTFREE"; break;
        case FDJACOBIAN : os << "FDJACOBIAN";  break;
        case DELTAMATJACOBIAN_LAST_TS  : os << "DELTAMATJACOBIAN_LAST_TS";   break;
        case DELTAMATJACOBIAN_LAST_IT  : os << "DELTAMATJACOBIAN_LAST_IT";   break;
        case REUSED  : os << "REUSED";   break;
        case HYSTFREE_GLOBALLY_SCALED_B_VERSION  : os << "HYSTFREE_GLOBALLY_SCALED_B_VERSION";   break;
        case HYSTFREE_LOCALLY_SCALED_B_VERSION  : os << "HYSTFREE_LOCALLY_SCALED_B_VERSION";   break;
        case HYSTFREE_LOCALLY_SCALED_H_VERSION  : os << "HYSTFREE_LOCALLY_SCALED_H_VERSION";   break;
        case HYSTFREE_SCALED_H_VERSION  : os << "HYSTFREE_SCALED_H_VERSION";   break;
        default    : os.setstate(std::ios_base::failbit);
      }
    } else {
      switch(flag)
      {
        case NOTSET   : os << "NOTSET > No information set";    break;
        case HYSTFREE: os << "HYSTFREE > No Jacobian approximation"; break;
        case FDJACOBIAN : os << "FDJACOBIAN > Finite Difference approximation of Jacobian";  break;
        case DELTAMATJACOBIAN_LAST_TS  : os << "DELTAMATJACOBIAN_LAST_TS > Approximation of Jacobian using input and output deltas of hyst-operator (towards last TS)";   break;
        case DELTAMATJACOBIAN_LAST_IT  : os << "DELTAMATJACOBIAN_LAST_IT > Approximation of Jacobian using input and output deltas of hyst-operator (towards last IT)";   break;
        case REUSED  : os << "REUSED > previous state of Jacobian approximation loaded";   break;
        case HYSTFREE_GLOBALLY_SCALED_B_VERSION  : os << "HYSTFREE_GLOBALLY_SCALED_B_VERSION > no Jacobian; material tensor scaled for better fixpoint convergence (B-Version)";   break;
        case HYSTFREE_LOCALLY_SCALED_B_VERSION  : os << "HYSTFREE_LOCALLY_SCALED_B_VERSION > Jacobian computed once each time step to estimate current slope for better fixpoint convergence (B-Version)";   break;
        case HYSTFREE_SCALED_H_VERSION  : os << "HYSTFREE_SCALED_H_VERSION > no Jacobian; material tensor scaled for better fixpoint convergence (H-Version)";   break;
        case HYSTFREE_LOCALLY_SCALED_H_VERSION  : os << "HYSTFREE_LOCALLY_SCALED_H_VERSION > Jacobian computed once each time step to estimate current slope for better fixpoint convergence (H-Version)";   break;
        default    : os.setstate(std::ios_base::failbit);
      }
    }

    return os.str();
  }

  std::string SolveStepHyst::vectorFlagToString(vectorFlag flag, int verbose){
    std::stringstream os;
    if(verbose == 0){
      switch(flag)
      {
        case NOTSET_VEC   : os << "NOTSET_VEC";    break;
        case PREVIOUS_TS_RES: os << "PREVIOUS_TS_RES"; break;
        case CURRENT_RES : os << "CURRENT_RES";  break;
        case PREVIOUS_IT_RES  : os << "PREVIOUS_IT_RES";   break;
        case RES_REM  : os << "RES_REM";   break;
        default    : os.setstate(std::ios_base::failbit);
      }
    } else {
      switch(flag)
      {
        case NOTSET_VEC   : os << "NOTSET_VEC; no information set";    break;
        case PREVIOUS_TS_RES: os << "PREVIOUS_TS_RES; hyst operator evaluated at PREVIOUS TS-solution; subtract solution from PREVIOUS TS > INCREMENTAL STEP TO PREVIOUS TS"; break;
        case CURRENT_RES : os << "CURRENT_RES; hyst operator evaluated at CURRENT IT-solution; subtract solution from CURRENT IT > INCREMENTAL STEP TO CURRENT IT"; break;
        case PREVIOUS_IT_RES  : os << "PREVIOUS_IT_RES; hyst operator evaluated at PREVIOUS IT-solution; subtract solution from CURRENT IT  > INCREMENTAL STEP TO CURRENT IT";   break;
        case RES_REM  : os << "RES_REM; hyst operator evaluated at ZEROVEC; subtract ZEROVEC > FULL STEP";   break;
        default    : os.setstate(std::ios_base::failbit);
      }
    }

    return os.str();
  }

  std::string SolveStepHyst::FEMatrixTypeToString(FEMatrixType matType){
    std::stringstream os;
    switch (matType){
      case NOTYPE : os << "NOTYPE"; break;
      case SYSTEM : os << "SYSTEM"; break;
      case STIFFNESS : os << "STIFFNESS"; break;
      case DAMPING : os << "DAMPING"; break;
      case CONVECTION : os << "CONVECTION"; break;
      case MASS : os << "MASS"; break;
      case AUXILIARY : os << "AUXILIARY"; break;
      case STIFFNESS_UPDATE : os << "STIFFNESS_UPDATE"; break;
      case DAMPING_UPDATE : os << "DAMPING_UPDATE"; break;
      case MASS_UPDATE : os << "MASS_UPDATE"; break;
      case SYSTEM_HYSTFREE : os << "SYSTEM_HYSTFREE"; break;
      case SYSTEM_FIXPOINT : os << "SYSTEM_FIXPOINT"; break;
      case SYSTEM_FD_JACOBIAN : os << "SYSTEM_FD_JACOBIAN"; break;
      case SYSTEM_DELTAMAT_JACOBIAN : os << "SYSTEM_DELTAMAT_JACOBIAN"; break;
      case BACKUP : os << "SYSTEM_BAK"; break;
      default: os.setstate(std::ios_base::failbit);
    }

    return os.str();
  }


} // end of namespace
