// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "solver/diagsolver.hh"

namespace OLAS {


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
				   const BaseVector &rhs, BaseVector &sol ) {

    // Tracing information
    (*cla) << "### Solver for diagonal system  matrix" << std::endl;

    // just apply a jacobi-preconditioner
    if ( precond.GetPrecondType() == JACOBI ) {
      precond.Apply( sysmat, rhs, sol );
    }
    else {
      Error("Diagobal solver needs Jacobi-preconditioner"__FILE__,__LINE__);
    }

  }

}
