#ifndef OLAS_COCR_HH
#define OLAS_COCR_HH

#include "BaseSolver.hh"

namespace CoupledField {

  //! Conjugate Orthogonal Conjugate Residual (COCR) method
  //!
  //! Iterative Krylov solver for complex symmetric linear systems Ax = b
  //! (A = A^T but in general A != A^H). COCR is the complex-symmetric
  //! analogue of CR; for real symmetric A it degenerates to ordinary CR.
  //!
  //! Reference:
  //!   [1] T. Sogabe, S.-L. Zhang,
  //!   "A COCR method for solving complex symmetric linear systems",
  //!   J. Comp. Appl. Math. 199 (2007) 297-303.
  //!
  //! Crucial implementation detail: every dot product in COCR is the
  //! BILINEAR form x^T y (no complex conjugation), realised via
  //! Vector<T>::operator* (see Vector.hh). The Hermitian inner product
  //! BaseVector::Inner() would conjugate the left operand and is NOT what
  //! the algorithm requires.
  //!
  //! Preconditioning: the algorithm uses the standard openCFS interface
  //! ptPrecond_->Apply(sysMat, r, z) which solves M z = r
  //!
  //! XML parameters read from the solver node:
  //!   - stoppingRule       : absNorm | relNormRHS | relNormRes0
  //!   - tol                : tolerance for the stopping criterion (default 1e-4)
  //!   - maxIter            : maximum number of iterations (default 2000)
  //!   - consoleConvergence : print ||r|| at every iteration (default no)
  //!
  //! Report keys written to info node:
  //!   - numIter, finalNorm, solutionIsOkay
  //!   - statistics/avgIterations, statistics/avgResReduction
  //!
  //! \note This first cut supports StdMatrix (single-field) only. Passing
  //! an SBM_Matrix raises an EXCEPTION.
  template<typename T>
  class COCRSolver : public BaseIterativeSolver {

  public:

    COCRSolver( PtrParamNode solverNode, PtrParamNode olasInfo );
    ~COCRSolver();

    //! Allocate / resize auxiliary vectors for the given system matrix
    void Setup( BaseMatrix &sysMat ) override;

    //! Not implemented for this solver
    void SetNewMatrixPattern() {
      EXCEPTION( "SetNewMatrixPattern not implemented for COCRSolver!" );
    }

    //! Solve sysMat * sol = rhs using preconditioned COCR
    void Solve( const BaseMatrix &sysMat, const BaseVector &rhs, BaseVector &sol );

    SolverType GetSolverType() { return COCR; }

  private:

    // Default values for parameters; will be overwritten if specified in XML
    //! Maximum iteration count
    Integer maxIter_ = 5000;
    //! Stopping tolerance
    Double  eps_ = 1e-4;
    //! Print ||r|| each iteration to console
    bool    consoleConvergence_ = false;

    //! Auxiliary vectors (six, following the LIS variant)
    BaseVector *r_ = nullptr;     //!< current residual r = b - A x
    BaseVector *p_ = nullptr;     //!< search direction
    BaseVector *q_ = nullptr;     //!< q = A p
    BaseVector *z_ = nullptr;     //!< preconditioned-residual recurrence z = M^-1 r
    BaseVector *qtld_ = nullptr;  //!< qtld = M^-1 q (preconditioner output)
    BaseVector *az_ = nullptr;    //!< az = A z

    //! Cached vector length; triggers re-allocation when problem size grows
    UInt vectorLength_ = 0;

    //! Disabled default constructor
    COCRSolver() {};
  };

}

#endif // OLAS_COCR_HH
