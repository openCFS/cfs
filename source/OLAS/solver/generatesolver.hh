// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_GENERATESOLVER_HH
#define OLAS_GENERATESOLVER_HH

#include "General/environment.hh"

//! \file generatesolver.hh
//! This module handles generation of solver objects. It is also responsible
//! for the instantiation of the templated solvers.

namespace CoupledField {

  class BaseSolver;
  class BaseMatrix;
  class ParamNode;
  class OLAS_Params;
  class OLAS_Report;
  
  //! Generate a basic solver object

  //! This method will generate a BaseSolver solver object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat    %Matrix that is preconditioned
  //! \param solver Type of desired solver
  //! \param xml    Pointer to ParamNode of <solver> section    
  //! \param params Pointer to a parameter object with steering parameters
  //!               for the solver that is to be generated
  //! \param report Pointer to report object into which the generated solver
  //!               should write its solutiopn report.
  BaseSolver* GenerateSolverObject( const BaseMatrix &mat, SolverType solver,
                                    ParamNode* xml, OLAS_Params *params, OLAS_Report *report );

}

#endif
