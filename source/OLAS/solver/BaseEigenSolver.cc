#include "DataInOut/ParamHandling/ParamNode.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "Utils/Timer.hh"

namespace CoupledField {

  static EnumTuple eigenSolverTypeTuples[] = 
  {
    EnumTuple( BaseEigenSolver::NO_EIGENSOLVER, "noEigenSolver" ),
    EnumTuple( BaseEigenSolver::ARPACK, "arpack" ),
    EnumTuple( BaseEigenSolver::PHIST, "phist"),
    EnumTuple( BaseEigenSolver::FEAST, "feast" ),
    EnumTuple( BaseEigenSolver::PALM, "palm" ),
  };

  // unbelievable how complicated easy stuff can be written :(
  Enum<BaseEigenSolver::EigenSolverType> BaseEigenSolver::eigenSolverType = \
  Enum<BaseEigenSolver::EigenSolverType>("EigenSolver Types",  sizeof(eigenSolverTypeTuples) / sizeof(EnumTuple),   eigenSolverTypeTuples);

  void BaseEigenSolver::PostInit()
  {
    // Assert that info Node is set
    assert( info_ );

    setupTimer_ = boost::shared_ptr<Timer>(new Timer("setup_" + eigenSolverType.ToString(eigenSolverName_)));
    info_->Get(ParamNode::SUMMARY)->Get("setup/timer")->SetValue(setupTimer_);

    solveTimer_ = boost::shared_ptr<Timer>(new Timer("solve_" + eigenSolverType.ToString(eigenSolverName_)));
    info_->Get(ParamNode::SUMMARY)->Get("solve/timer")->SetValue(solveTimer_);
  }

}
