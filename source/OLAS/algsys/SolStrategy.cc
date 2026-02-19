#include "SolStrategy.hh"
namespace CoupledField {



  SolStrategy::SolStrategy( PtrParamNode param ) {
    param_ = param;
    type_ = NO_STRATEGY;
    numSolSteps_ = 1;
    curSolStep_ = 1;

    if(param->Has("matrix/estimated_row_size"))
      estimated_row_size_  = param->Get("matrix/estimated_row_size")->As<unsigned int>();

    // matrix node for static condensation
    {
      statCondMatNode_.reset(new ParamNode());
      statCondMatNode_->SetName("matrix");
      statCondMatNode_->Get("storage",ParamNode::INSERT)->SetValue("variableBlockRow");
      statCondMatNode_->Get("reordering",ParamNode::INSERT)->SetValue("noReordering");
    }

    // default no multiharmonic analysis
    isMultHarm_ = false;

    // default no 2.5d harmonic analysis
    is25dHarm_ = false;
  }

  SolStrategy::~SolStrategy() {

  }


  shared_ptr<SolStrategy> SolStrategy::Generate( PtrParamNode node) {
    shared_ptr<SolStrategy> ret;

    // ensure valid xml element
    if( !node) {
      EXCEPTION("Pointer to element <solutionStrategy> is empty!")
    }

    // check if strategy is set. If not, generate SolStrategyStd
    if( !node->HasChildren() ) {
      PtrParamNode stdNode(new ParamNode());
      stdNode->SetName("standard");
      node->AddChildNode(stdNode);
    }

    // otherwise check type of element
    PtrParamNode stratNode = node->GetChild();
    std::string strat = stratNode->GetName();
    if( strat == "standard" ) {
      ret.reset(new SolStrategyStd(stratNode));
      } else if( strat == "twoLevel") {
      ret.reset(new SolStrategyTwoLevel(stratNode));
    } else {
      EXCEPTION("Solution strategy '" << strat << "' is not known.");
    }

    return ret;
  }

  // ==========================================================================
  //  S T A N D A R D   S O L U T I O N   S T R A T E G Y
  // ==========================================================================
  
  SolStrategyStd::SolStrategyStd(PtrParamNode node ) 
  : SolStrategy(node) {
    
    type_ = STD_STRATEGY;
    numSolSteps_ = 1;
    curSolStep_ = 1;
    
    // check, if subnodes are present. if not, generate them, but
    // without default attributes
    
    // <setup> element
    if( !node->Has("setup")) {
      setupNode_.reset(new ParamNode());
      setupNode_->SetName("setup");
      param_->AddChildNode(setupNode_);
    } else {
      setupNode_ = param_->Get("setup");
    }
    
    // <matrix> element
    if( !node->Has("matrix")) {
      matrixNode_.reset(new ParamNode());
      matrixNode_->SetName("matrix");
      param_->AddChildNode(matrixNode_);
    } else {
      matrixNode_ = param_->Get("matrix");
    }

    // <eigenSolver> element
    if( !node->Has("eigenSolver")) {
        esNode_.reset(new ParamNode());
        esNode_->SetName("eigenSolver");
        param_->AddChildNode(esNode_);
      } else {
        esNode_ = param_->Get("eigenSolver");
      }
    
    // <solver> element
    if( !node->Has("solver")) {
      solverNode_.reset(new ParamNode());
      solverNode_->SetName("solver");
      param_->AddChildNode(solverNode_);
    } else {
      solverNode_ = param_->Get("solver");
    }

    // <precond> element
    if( !node->Has("precond")) {
      precondNode_.reset(new ParamNode());
      precondNode_->SetName("precond");
      param_->AddChildNode(precondNode_);
    } else {
      precondNode_ = param_->Get("precond");
    }
    
    // <nonLinear> element
    if( !node->Has("nonLinear")) {
      nonlinNode_.reset(new ParamNode());
      nonlinNode_->SetName("nonLinear");
      param_->AddChildNode(nonlinNode_);
    } else {
      nonlinNode_ = param_->Get("nonLinear");
    }
    
    // <hysteresis> element
    if( !node->Has("hysteretic")) {
      hystNode_.reset(new ParamNode());
      hystNode_->SetName("hysteretic");
      param_->AddChildNode(hystNode_);
    } else {
      hystNode_ = param_->Get("hysteretic");
    }

    // <timeStepping> element
    if( !node->Has("timeStepping")) {
      tsNode_.reset(new ParamNode());
      tsNode_->SetName("timeStepping");
      param_->AddChildNode(tsNode_);
    } else {
      tsNode_ = param_->Get("timeStepping");
    }
    
   }

  SolStrategyStd::~SolStrategyStd(){

  }

  UInt SolStrategyStd::GetNumSolSteps(){
    return 1;
  }

  void SolStrategyStd::SetActSolStep(UInt stepNum){
    // check, that solstep does not get changed
  }
  
  UInt SolStrategyStd::GetActSolStep() {
    return 1;
  }

  bool SolStrategyStd::UseStaticCondensation(){
    bool useCondens = false;
    setupNode_->GetValue("staticCondensation", useCondens, ParamNode::INSERT);
    return useCondens;
  }
  

  UInt SolStrategyStd::GetNumSBMBlocks(){
    if( UseStaticCondensation() ) {
      return 2;
    } else {
      return 1;
    }
  }
  
  bool SolStrategyStd::UseDirichletPenalty() {
    std::string usageString = "elimination";
    setupNode_->GetValue("idbcHandling", usageString, ParamNode::INSERT);
    if( usageString == "elimination")
      return false;
    else 
      return true;
  }
  
  bool SolStrategyStd::CalcConditionNumber() {
    bool calcCond = false;
    setupNode_->GetValue("calcConditionNumber", calcCond, ParamNode::INSERT);
    return calcCond;
  }
  
  PtrParamNode SolStrategyStd::GetSetupNode(){
      return setupNode_;
    }
  
  BaseMatrix::StorageType SolStrategyStd::GetStorageType(UInt blockNum) {
    
    std::string formatString;

    // return 
    if( UseStaticCondensation() && blockNum==1 ) {
      statCondMatNode_->GetValue("storage",formatString);
    } else {
      GetMatrixNode(blockNum)->GetValue("storage",formatString);
    }
    
    BaseMatrix::StorageType st = BaseMatrix::storageType.Parse(formatString);
    return st;
    
  }
  

  PtrParamNode SolStrategyStd::GetMatrixNode(UInt level){
    PtrParamNode ret;
    if( UseStaticCondensation() && level == 1 ) {
      ret = statCondMatNode_;
    } else {
      ret = matrixNode_;
    }
    return ret; 
  }

  std::string SolStrategyStd::GetSolverId(){
    std::string id = "default";
    solverNode_->GetValue("id", id, ParamNode::INSERT);
    return id;
  }
  
  std::string SolStrategyStd::GetEigenSolverId(){
    std::string id = "default";
    esNode_->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  std::string SolStrategyStd::GetPrecondId(){
    std::string id = "default";
    precondNode_->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  PtrParamNode SolStrategyStd::GetNonLinNode(){
    return nonlinNode_;
  }

  PtrParamNode SolStrategyStd::GetHystNode(){
    return hystNode_;
  }

  PtrParamNode SolStrategyStd::GetTimeSteppingNode(){
    PtrParamNode ret;
    return ret;
  }
  
  
  
  // ==========================================================================
  //  T W O - L E V E L   S O L U T I O N   S T R A T E G Y
  // ==========================================================================
  SolStrategyTwoLevel::SolStrategyTwoLevel(PtrParamNode node ) 
  : SolStrategy(node) {
    type_ = TWO_LEVEL_STRATEGY;

    // <setup> element
    if( !node->Has("setup")) {
      setupNode_.reset(new ParamNode());
      setupNode_->SetName("setup");
      param_->AddChildNode(setupNode_);
    } else {
      setupNode_ = param_->Get("setup");
    }
    
    // read in splitting node
    PtrParamNode splitNode = node->Get("splitting", ParamNode::PASS);
    if( !splitNode ) {
      splitNode.reset(new ParamNode());
      splitNode->SetName("splitting");
      node->AddChildNode(splitNode);
    }
    
    // read in <matrix>-subnodes
    matrixNodes_.Resize(2);
    // Loop over all level
    for( UInt iLevel = 0; iLevel < 2; ++iLevel ) {
      PtrParamNode levelNode = splitNode->GetByVal("level","num",iLevel+1,ParamNode::INSERT);
      matrixNodes_[iLevel] = levelNode->Get("matrix", ParamNode::INSERT);
    }

    // Obtain solutionNode
    PtrParamNode solNode = node->Get("solution");
    ParamNodeList stepNodes = solNode->GetChildren();
    numSolSteps_ = stepNodes.GetSize();
    
    
    // Resize vectors
    eigenSolverNodes_.Resize(numSolSteps_);
    solverNodes_.Resize(numSolSteps_);
    precondNodes_.Resize(numSolSteps_);
    nonLinNodes_.Resize(numSolSteps_);
    hystNodes_.Resize(numSolSteps_);
    timeStepNodes_.Resize(numSolSteps_);
    
    
    // loop over all solution steps nodes
    for( UInt i = 0; i < stepNodes.GetSize(); ++i ) {

      PtrParamNode stepNode = stepNodes[i];
      
      // <eigenSolver> element
      eigenSolverNodes_[i] = stepNode->Get("eigenSolver",ParamNode::INSERT);

      // <solver> element
      solverNodes_[i] = stepNode->Get("solver",ParamNode::INSERT);
      
      // <precond> element
      precondNodes_[i] = stepNode->Get("precond",ParamNode::INSERT);
      
      // <nonLinear> element
      nonLinNodes_[i] = stepNode->Get("nonLinear",ParamNode::INSERT);

      // <hysteretic> element
      hystNodes_[i] = stepNode->Get("hysteretic",ParamNode::INSERT);

      // <timeStepping> element
      timeStepNodes_[i] = stepNode->Get("timeStepping", ParamNode::INSERT);
    } // loop over all solution steps
    
    // in the end, set active solution step to 1st step
    curSolStep_ = 0;
    
//    std::cerr << "=====================\n"
//        << "Param Node at End of SolStrategy Setup:\n"
//        << "=====================\n";
//    param_->Dump();
        
  }
  
  SolStrategyTwoLevel::~SolStrategyTwoLevel(){

  }

  UInt SolStrategyTwoLevel::GetNumSolSteps(){
    return numSolSteps_;
  }

  void SolStrategyTwoLevel::SetActSolStep(UInt stepNum ){
    
    if( stepNum > numSolSteps_ ) {
      EXCEPTION("Can not set active solution step to " << stepNum 
                << ", as solution strategy has only " << numSolSteps_ 
                << " solution steps." );
    }
    curSolStep_ = stepNum-1;
  }

  UInt SolStrategyTwoLevel::GetActSolStep() {
    return curSolStep_+1;
  }
  
  bool SolStrategyTwoLevel::UseStaticCondensation(){
    bool useCondens = false;
    setupNode_->GetValue("staticCondensation", useCondens, ParamNode::INSERT);
    return useCondens;
  }
  
  UInt SolStrategyTwoLevel::GetNumSBMBlocks(){
    if( UseStaticCondensation() ) {
      return 3;
    } else {
      return 2;
    }
  }

  bool SolStrategyTwoLevel::UseDirichletPenalty() {
    std::string usageString = "elimination";
    setupNode_->GetValue("idbcHandling", usageString, ParamNode::INSERT);
    if( usageString == "elimination")
      return false;
    else 
      return true;
  }

  bool SolStrategyTwoLevel::CalcConditionNumber() {
    bool calcCond = false;
    setupNode_->GetValue("calcConditionNumber", calcCond, ParamNode::INSERT);
    return calcCond;
  }

  PtrParamNode SolStrategyTwoLevel::GetSetupNode(){
    return setupNode_;
  }

  BaseMatrix::StorageType SolStrategyTwoLevel::GetStorageType(UInt blockNum) {

    std::string formatString;

    // return 
    if( UseStaticCondensation() && blockNum==2 ) {
      statCondMatNode_->GetValue("storage",formatString);
    } else {
      GetMatrixNode(blockNum)->GetValue("storage",formatString);
    }

    BaseMatrix::StorageType st = BaseMatrix::storageType.Parse(formatString);
    return st;

  }


  PtrParamNode SolStrategyTwoLevel::GetMatrixNode(UInt blockNum){
    PtrParamNode ret;
    
    if( blockNum >= GetNumSBMBlocks() ) {
      EXCEPTION("Can not return matrix node for SBM block number " << blockNum+1
                << " as the system has only " << GetNumSBMBlocks() 
                << " SBM blocks." );
    }
    
    if( UseStaticCondensation() == true && 
        blockNum == 2 ) {
      ret = statCondMatNode_;
    } else {
      ret = matrixNodes_[blockNum];
    }
    return ret; 
  }

  std::string SolStrategyTwoLevel::GetSolverId(){
    std::string id = "default";
    solverNodes_[curSolStep_]->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  std::string SolStrategyTwoLevel::GetEigenSolverId(){
    std::string id = "default";
    eigenSolverNodes_[curSolStep_]->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  std::string SolStrategyTwoLevel::GetPrecondId(){
    std::string id = "default";
    precondNodes_[curSolStep_]->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  PtrParamNode SolStrategyTwoLevel::GetNonLinNode(){
    return nonLinNodes_[curSolStep_];
  }

  PtrParamNode SolStrategyTwoLevel::GetHystNode(){
    return hystNodes_[curSolStep_];
  }

  PtrParamNode SolStrategyTwoLevel::GetTimeSteppingNode(){
    PtrParamNode ret;
    return ret;
  }

  

// ************************************************************************
// ENUM INITIALIZATION
// ************************************************************************

// Definition of stratey types
static EnumTuple strategyTuples[] = {
       EnumTuple(SolStrategy::NO_STRATEGY,  "NoStrategy"),
       EnumTuple(SolStrategy::STD_STRATEGY, "StdStrategy"),
       EnumTuple(SolStrategy::TWO_LEVEL_STRATEGY,  "TwoLevelStrategy")
};
Enum<SolStrategy::StrategyType> SolStrategy::strategyType = \
     Enum<SolStrategy::StrategyType>("Solution Strategy Types",
         sizeof(strategyTuples) / sizeof(EnumTuple),
         strategyTuples);
}
