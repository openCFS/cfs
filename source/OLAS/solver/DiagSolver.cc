#include "MatVec/Vector.hh"

#include "OLAS/precond/BasePrecond.hh"
#include "OLAS/solver/DiagSolver.hh"
#include "Utils/Timer.hh"
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

    // Tracing information
    (*cla) << "### Solver for diagonal system  matrix" << std::endl;

    // just apply a jacobi-preconditioner
    if ( ptPrecond_->GetPrecondType() == BasePrecond::JACOBI ) {
      ptPrecond_->GetPrecondTimer()->Start();
      ptPrecond_->Apply( sysmat, rhs, sol );
      ptPrecond_->GetPrecondTimer()->Stop();
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
