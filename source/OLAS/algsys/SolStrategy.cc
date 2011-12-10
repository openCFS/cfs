#include "SolStrategy.hh"
namespace CoupledField {



  SolStrategy::SolStrategy( PtrParamNode param ) {
    param_ = param;
    type_ = NO_STRATEGY;
    numSolSteps_ = 0;
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
      //} else if( strat == "twoLevel") {
      //ret.reset(new SolStrategyStd(stratNode));
    } else {
      EXCEPTION("Solution strategy '" << strat << "' is not known.");
    }

    return ret;
  }

  // ==========================================================================
  //  STANDARD SOLUTION STRATEGY
  // ==========================================================================
  
  
  SolStrategyStd::SolStrategyStd(PtrParamNode node ) 
  : SolStrategy(node) {
    
    type_ = STD_STRATEGY;
    
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
    
    // <exportLinSys> element
    if( !node->Has("exportLinSys")) {
      exportNode_.reset(new ParamNode());
      exportNode_->SetName("exportLinSys");
      param_->AddChildNode(exportNode_);
    } else {
      exportNode_ = param_->Get("exportLinSys");
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

  void SolStrategyStd::SetActSolStep(){
    // check, that solstep does not get changed
  }

  bool SolStrategyStd::UseStaticCondensation(){
    bool useCondens = false;
    setupNode_->GetValue("staticCondensation", useCondens, ParamNode::INSERT);
    return useCondens;
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
    
    // if we have static condensation active and this block is queried
    // just return the &= operation of all previous matrices
    
    // in all other cases, just return the matrix type set at the
    // matrix element
    std::string formatString;
    GetMatrixNode(blockNum)->GetValue("storage",formatString);
    BaseMatrix::StorageType st = BaseMatrix::storageType.Parse(formatString);
    return st;
    
  }
  

  PtrParamNode SolStrategyStd::GetMatrixNode(UInt level){
    return matrixNode_; 
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

  std::string SolStrategyStd::GetPrecondId(UInt level){
    std::string id = "default";
    precondNode_->GetValue("id", id, ParamNode::INSERT);
    return id;
  }

  PtrParamNode SolStrategyStd::GetExportLinSysNode(){
    return exportNode_;
  }

  PtrParamNode SolStrategyStd::GetNonLinNode(){
    PtrParamNode ret;
    return ret;
  }
  PtrParamNode SolStrategyStd::GetTimeSteppingNode(){
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
