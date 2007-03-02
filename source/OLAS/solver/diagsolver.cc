#include "solver/diagsolver.hh"

namespace OLAS {


  // **************
  //   Destructor
  // **************
  template<typename T>
  DiagSolver<T>::~DiagSolver() {
    ENTER_FCN("DiagSolver::~DiagSolver");
  }


  // ****************
  //   Solve method
  // ****************
  template<typename T>
  void DiagSolver<T>::Solve( const BaseMatrix &sysmat,
				   const BasePrecond &precond,
				   const BaseVector &rhs, BaseVector &sol ) {

    // Tracing information
    ENTER_FCN("DiagSolver::Solve");
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
