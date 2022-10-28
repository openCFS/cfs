#include "MatVec/Vector.hh"

#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/DiagSolver.hh"
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
                             const BaseVector &rhs, BaseVector &sol ) {

    // just apply a jacobi-preconditioner
    if ( ptPrecond_->GetPrecondType() == BasePrecond::JACOBI ) {
      ptPrecond_->Apply( sysmat, rhs, sol );
    }
    else {
      EXCEPTION("Diagonal solver needs Jacobi-preconditioner");
    }

  }

  // Explicit template instantiation
  template class DiagSolver<Double>;
  template class DiagSolver<Complex>;
  
}
