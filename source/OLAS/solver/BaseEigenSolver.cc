#include "DataInOut/ParamHandling/ParamNode.hh"

#include "OLAS/solver/BaseEigenSolver.hh"

namespace CoupledField {

  static EnumTuple eigenSolverTypeTuples[] = 
  {
    EnumTuple( BaseEigenSolver::NO_EIGENSOLVER, "noEigenSolver" ),
    EnumTuple( BaseEigenSolver::ARPACK, "arpack" ),
  };

  Enum<BaseEigenSolver::EigenSolverType> BaseEigenSolver::eigenSolverType = \
  Enum<BaseEigenSolver::EigenSolverType>("EigenSolver Types",  sizeof(eigenSolverTypeTuples) / sizeof(EnumTuple),   eigenSolverTypeTuples);
}
