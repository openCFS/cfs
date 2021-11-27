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
    EnumTuple( BaseEigenSolver::QUADRATIC, "quadratic" ),
  };

  // unbelievable how complicated easy stuff can be written :(
  Enum<BaseEigenSolver::EigenSolverType> BaseEigenSolver::eigenSolverType = \
  Enum<BaseEigenSolver::EigenSolverType>("EigenSolver Types",  sizeof(eigenSolverTypeTuples) / sizeof(EnumTuple),   eigenSolverTypeTuples);

  void BaseEigenSolver::PostInit()
  {
    // Assert that info Node is set
    assert(info_);

    std::string name = eigenSolverType.ToString(eigenSolverName_);

    setupTimer_ = info_->Get(ParamNode::SUMMARY)->Get("setup_" + name + "/timer")->AsTimer();

    solveTimer_ = info_->Get(ParamNode::SUMMARY)->Get("solve_" + name + "/timer")->AsTimer();
  }

}
