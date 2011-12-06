#ifndef OLAS_GENERATESOLVER_HH
#define OLAS_GENERATESOLVER_HH

#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/BaseMatrix.hh"
#include "OLAS/solver/BaseSolver.hh"

//! \file generatesolver.hh
//! This module handles generation of solver objects. It is also responsible
//! for the instantiation of the templated solvers.

namespace CoupledField {

  class BaseSolver;
  class BaseMatrix;
  class SolStrategy;
  
  //! Generate a basic solver object

  //! This method will generate a BaseSolver solver object that fits
  //! to the input matrix and return a pointer to that object.
  //! \param mat      %Matrix that is preconditioned
  //! \param strat    Pointer to solution strategy object
  //! \param xml      Pointer to ParamNode of <solverList> section
  //! \param olasInfo Base below "OLAS" in info.xml
  BaseSolver* GenerateSolverObject( const BaseMatrix &mat,
                                    shared_ptr<SolStrategy> strat,
                                    PtrParamNode xml, 
                                    PtrParamNode olasInfo );

  //! Return list of compatible matrix types

  //! This method returns a list of all matrix storage types the solver can 
  //! handle. 
  //! \note In case the solver can handle any type of sparse matrix (i.e.
  //! Krylov based solvers, requiring just matrix-vector operations),
  //! the return set will be empty!
  std::set<BaseMatrix::StorageType> 
  GetSolverCompatMatrixFormats(BaseSolver::SolverType );

  //! Return if preconditioner is capable of solving a SBM system
  bool IsSolverSBMCapable(BaseSolver::SolverType );

}

#endif
