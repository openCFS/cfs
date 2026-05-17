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
  //!   [2] H. A. van der Vorst, Q. Ye,
  //!   "Residual Replacement Strategies for Krylov Subspace Iterative
  //!   Methods for the Convergence of True Residuals",
  //!   SIAM J. Sci. Comput. 22(3) (2000) 835-852.
  //!   (for the residual-replacement strategy used in recomputeTrueResEvery)
  //!
  //! The implementation follows the practical preconditioned variant used
  //! in the LIS library (lis_cocr), which is mathematically equivalent to
  //! Algorithm 1 of the paper but maintains a separate z-recurrence so
  //! that each iteration costs exactly one matrix-vector product, one
  //! preconditioner solve, three bilinear dot products and a handful of
  //! AXPYs.
  //!
  //! Crucial implementation detail: every dot product in COCR is the
  //! BILINEAR form x^T y (no complex conjugation), realised via
  //! Vector<T>::operator* (see Vector.hh). The Hermitian inner product
  //! BaseVector::Inner() would conjugate the left operand and is NOT what
  //! the algorithm requires.
  //!
  //! Preconditioning: the algorithm uses the standard openCFS interface
  //! ptPrecond_->Apply(sysMat, r, z) which solves M z = r. Compatible
  //! preconditioners include ILDL0, the full ILDL family, and any direct
  //! solver (e.g. PARDISO) applied as a preconditioner. For symmetric
  //! preconditioners M = L D L^T this gives the correct symmetric COCR.
  //!
  //! XML parameters read from the solver node:
  //!   - stoppingRule       : absNorm | relNormRHS | relNormRes0
  //!   - tol                : tolerance for the stopping criterion (default 1e-6)
  //!   - maxIter            : maximum number of iterations (default 1000)
  //!   - resDirectly        : recompute true residual every N iters
  //!                          (0 = never, default 0)
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
    void Setup( BaseMatrix &sysMat );

    //! Not implemented for this solver
    void SetNewMatrixPattern() {
      EXCEPTION( "SetNewMatrixPattern not implemented for COCRSolver!" );
    }

    //! Solve sysMat * sol = rhs using preconditioned COCR
    void Solve( const BaseMatrix &sysMat,
                const BaseVector &rhs, BaseVector &sol );

    SolverType GetSolverType() { return COCR; }

  private:

    //! Auxiliary vectors (six, following the LIS variant)
    BaseVector *r_ = nullptr;     //!< current residual r = b - A x
    BaseVector *p_ = nullptr;     //!< search direction
    BaseVector *q_ = nullptr;     //!< q = A p
    BaseVector *z_ = nullptr;     //!< preconditioned-residual recurrence z = M^-1 r
    BaseVector *qtld_ = nullptr;  //!< qtld = M^-1 q (preconditioner output)
    BaseVector *az_ = nullptr;    //!< az = A z

    //! Cached vector length; triggers re-allocation when problem size grows
    UInt vectorLength_ = 0;

    //! Internal allocator usable from Solve (Setup uses non-const sysMat)
    void PrivateSetup( const BaseMatrix &sysMat );

    //! Bilinear (non-Hermitian) dot product: returns x^T y with no
    //! conjugation, using Vector<T>::operator*. Both x and y must be
    //! concrete Vector<T> instances (the dynamic_cast asserts this).
    T BilinearInner( const BaseVector &x, const BaseVector &y ) const;

    //! Disabled default constructor
    COCRSolver() {};
  };

}

#endif // OLAS_COCR_HH
