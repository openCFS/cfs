// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_GENERATE_EIGEN_SOLVER_HH
#define OLAS_GENERATE_EIGEN_SOLVER_HH

#include "General/environment.hh"
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
  //! \param solver Type of desired solver
  //! \param xml    Pointer to ParamNode of <solver> section
  //! \param params Pointer to a parameter object with steering parameters
  //!               for the solver that is to be generated
  //! \param report Pointer to report object into which the generated solver
  //!               should write its solutiopn report.
  BaseEigenSolver * GenerateEigenSolverObject( BaseMatrix &mat, 
                                               PtrParamNode xml,
                                               PtrParamNode   eigenInfo );

}

#endif
