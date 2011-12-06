#ifndef OLAS_GENERATE_EIGEN_SOLVER_HH
#define OLAS_GENERATE_EIGEN_SOLVER_HH

#include "General/Environment.hh"
#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

//! \file generateEigenSolver.hh
//! This module handles generation of Eigen solver objects. It is also
//! responsible for the instantiation of the templated Eigen solvers.

namespace CoupledField {

  // forward class declarations
  class BaseEigenSolver;
  class BaseMatrix;
  class ParamNode;

  //! Generate a basic solver object

  //! This method will generate a BaseEigenSolver solver object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat    %Matrix that is preconditioned
  //! \param strat  Pointer to solution strategy
  //! \param eSolverList Pointer to <eigenSolverList> element
  //! \param solverList Pointer to <solverList> element
  //! \param precondList Pointer to <precondList> element
  //! \param eigenInfo Output element in info-tree
  BaseEigenSolver * GenerateEigenSolverObject( BaseMatrix &mat, 
                                               shared_ptr<SolStrategy> strat,
                                               PtrParamNode eSolverList,
                                               PtrParamNode solverList,
                                               PtrParamNode precondList,
                                               PtrParamNode eigenInfo );

}

#endif
