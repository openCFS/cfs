// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MatVec/vector.hh"

#include "OLAS/precond/baseprecond.hh"
#include "OLAS/solver/diagsolver.hh"

namespace CoupledField {


  // **************
  //   Destructor
  // **************
  template<typename T>
  DiagSolver<T>::~DiagSolver() {
  }


  // ****************
  //   Solve method
  // ****************
  template<typename T>
  void DiagSolver<T>::Solve( const BaseMatrix &sysmat,
                             const BasePrecond &precond,
                             const BaseVector &rhs, BaseVector &sol, InfoNode* analysis_step ) {

    // Tracing information
    (*cla) << "### Solver for diagonal system  matrix" << std::endl;

    // just apply a jacobi-preconditioner
    if ( precond.GetPrecondType() == BasePrecond::JACOBI ) {
      precond.Apply( sysmat, rhs, sol );
    }
    else {
      EXCEPTION("Diagonal solver needs Jacobi-preconditioner");
    }

  }

// Explicit template instantiation
#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class DiagSolver<Double>;
  template class DiagSolver<Complex>;
#endif
  
}
