#ifndef OLAS_SYMMLQ_HH
#define OLAS_SYMMLQ_HH

#include "OLAS/utils/math/GivensRotation.hh"

#include "BaseSolver.hh"

namespace CoupledField {

  // Forward declaration of classes
  // class GivenRotation;

  //! Implementation of the Symmetric LQ Algorithm (SYMMLQ)

  //! This class implements the Symmetric LQ Algorithm (SYMMLQ) for the
  //! iterative solution of a linear system of equations \f$Ax=b\f$
  //! with a sparse and symmetric/Hermitean matrix
  //! \f[
  //! A \in \mathcal{C}^{(n\times n)} \enspace,
  //! \f]
  //! which need not be positive definite.
  //! Starting from an initial guess \f$x_0\f$ the SYMMLQ method iteratively
  //! builds a basis of the Krylov subspace
  //! \f[
  //! K^q(A;r_0) = \{ r_0, Ar_0, A^2 r_0, A^3r_0, \ldots, A^{q-1}r_0\}
  //! \enspace.
  //! \f]
  //! In each iteration step a new approximate solution
  //! \f[
  //! x_k\in x_0 + A K^{k}(A;r_0)
  //! \f]
  //! is determined such that the Euclidean norm of the new error
  //! \f$\|A^{-1}b - x_k\|_2\f$ is minimal.\n\n
  //! Since the matrix \f$A\f$ is symmetric/Hermitain the construction of an
  //! orthonormal/unitary basis of the Krylov subspace can be achieved by
  //! a 3-term recurrence using the Lanczos method.\n\n
  //! This class represents a straightforward implementation of SYMMLQ. It uses
  //! the modified Gram-Schmidt Orthogonalisation procedure for the generation
  //! of the basis of the underlying Krylov subspace and employes successive
  //! Givens rotations for solving the resulting least-squares problem.
  //! For further details
  //! see <b>Henk van der Vorst, "Iterative Krylov Methods for Large Linear
  //! Systems", Cambridge Monographs on Applied and Computational
  //! Mathematics</b>.\n\n
  //! The behaviour of the HypreSYMMLQSolver can be tuned via the following
  //! parameters which are read from the parameter object myParams_:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for Hypre SYMMLQ solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>SYMMLQ_maxIter</td>
  //!       <td align="center">&gt; = 1</td>
  //!       <td align="center">1</td>
  //!       <td>Maximal number of iteration steps allowed for the
  //!         SYMMLQ(m) method. If this is set to 1, then we are
  //!         actually performing a full SYMMLQ iteration and the
  //!         maximal number of iterations is determined by
  //!         the value of SYMMLQ_maxKrylovDim.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>SYMMLQ_epsilon</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">1e-6</td>
  //!       <td>The value of this parameter is used for the convergence
  //!         criterion in the stopping test of SYMMLQ. The stopping rule
  //!         is as follows
  //!         \f[
  //!         \|r^{(k)}\|  \leq \varepsilon \cdot \|b\|
  //!         \f]
  //!         where \f$r^{(k)}\f$ is the residual of the k-th step
  //!         and \f$b\f$ is the right hand side. If \f$\|b\|=0\f$,
  //!         then it is replaced by \f$\|r^{(0)}\|\f$.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>SYMMLQ_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>This parameter is used to control the verbosity of the method.
  //!           If it is set to 'yes', then information on the SYMMLQ
  //!           convergence, like e.g. the norm of the residual per iteration
  //!           step will be logged to the standard %OLAS report stream (*cla) ---removed logging---.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \n\n The above parameters are read from the parameter object myParams_
  //! every time a solve of the linear system is triggered. Once the system has
  //! been solved, the HypreSYMMLQSolver writes three values into the report
  //! object myReport_\n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="2" align="center">
  //!         <b>Values the SYMMLQ solver writes into myReport_</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Report</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>numIter</td>
  //!       <td>This is the number of iterations steps the SYMMLQ method has
  //!         performed.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>finalNorm</td>
  //!       <td>This is final norm of the relative residual:
  //!         \f$\displaystyle\frac{\|r^{(k)}\|_2}{\|b\|_2}\enspace.\f$
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \note
  //! - Currently the SYMMLQSolver class only supports problem matrices that
  //!   are instances of the BNMatrix class and have scalar entries of type
  //!   Double or Complex. See also the description of TestMatrixType() below.
  template<typename T> class SYMMLQSolver : public BaseIterativeSolver {

  public:

    //! Constructor

    //! This constructor initialises the pointers to the communication objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    SYMMLQSolver( PtrParamNode solverNode, PtrParamNode olasInfo );

    //! Default Destructor

    //! The default destructor handles de-allocation of memory of the
    //! internal data structures.
    ~SYMMLQSolver();

    //! Prepare internal data structures

    //! The implementation of the Setup method for the SYMMLQ solver performs
    //! mainly memory management tasks. It will prepare the internal data
    //! structures. If this is the first call to Setup or if maxKrylovDim
    //! has changed, it will allocate new memory, freeing previously allocated
    //! memory if necessary.
    void Setup( BaseMatrix &sysMat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Solve a linear system using SYMMLQ.

    //! Solve the linear system \f$Ax=b\f$ with the square matrix \f$A\f$ and
    //! right-hand side \f$b\f$ for \f$x\f$ using right preconditioned
    //! SYMMLQ(m).
    //! \param sysMat  system matrix \f$A\f$
    //! \param precond preconditioner to be used for right preconditioning
    //! \param rhs     right-hand side vector \f$b\f$
    //! \param sol     on input initial guess for the solution \f$x\f$, on
    //!                exit approximate solution
    void Solve( const BaseMatrix &sysMat, const BasePrecond &precond,
		const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return SYMMLQ
    SolverType GetSolverType() {
      return SYMMLQ;
    }

  private:

    //! Right-hand side vector of least-squares problem

    //! This auxillary vector is used in the solution of the least squares
    //! problem to store the right-hand side of the problem. It is defined
    //! as
    //! \f[ b = \|r^{(0)}\|_2 e_1 \f]
    //! where \f$e_1\f$ is the first canonical basis vector and \f$r^{(0)}\f$
    //! is the residual of the initial guess.
    Vector<T> bVec_;

    //! Auxilliary vector for applying preconditioner

    //! This auxillary vector is used in the application of the preconditioner.
    //! In the right-preconditioned SYMMLQ method we solve a linear system
    //! with a system matrix \f$\hat{A} = A M^{-1}\f$ where \f$A\f$ is the
    //! original system matrix and \f$M^{-1}\f$ is the preconditioner. This
    //! vector is used to store the result of applying the preconditioner to a
    //! vector \f$v\f$, i.e. \f$\mbox{pVec\_} = M^{-1}v\f$.
    Vector<T> pVec_;

    //! Maximal dimension of the generated Krylov subspace

    //! This integer is used to store the maximal dimension \f$m\f$ of the
    //! Krylov subspace \f$K^m(A;r_0)\f$ generated in the solution of the
    //! linear system. The value is obtained from the SYMMLQ_maxKrylovDim entry
    //! of the myParams_ object.
    UInt maxKrylovDim_;

    //! Dimension of the problem

    //! This attribute is used to keep track of the dimension of the problem,
    //! i.e. the dimension \f$n\f$ of the linear system,
    //! \f$A\in\mathcal{C}^{n\times n}\f$. The value of the attribute is set
    //! in Setup() and checked for changes in subsequent calls to Setup(),
    //! since a change indicates that the internal data structures must be
    //! re-sized.
    UInt problemDim_;

    //! "c" coefficients of Givens rotation

    //! We use a Double array to store the "c" coefficients of the Givens
    //! rotation. We have decided to use this light-weight data structure
    //! instead of a Vector<Double>.
    Double* c_;

    //! "s" coefficients of Givens rotation

    //! We use a Double/Complex array to store the "s" coefficients of the
    //! Givens rotation. We have decided to use this light-weight data
    //! structure instead of a Vector<Double/Complex>.
    T* s_;

    //! Data structure to store the upper Hessenberg matrix H

    //! This data structure is used to store the upper Hessenberg matrix H
    //! and the R part of its factorisation.
    //! \f[
    //! H = \left( \begin{array}{cccccc}
    //! h_{1,1} & h_{1,2} & \ldots & \ldots & h_{1,i-1} & h_{i,i}   \\
    //! h_{2,1} & h_{2,2} & \ldots & \ldots & h_{2,i-1} & h_{2,i}   \\
    //!    0    & h_{3,2} & \ddots &        &  \vdots   & \vdots    \\
    //!  \vdots & \ddots  & \ddots & \ddots &  \vdots   & \vdots    \\
    //!  \vdots &         & \ddots & \ddots & h_{i,i-1} & h_{i,i}   \\
    //!  \vdots &         &        & \ddots & h_{i,i-1} & h_{i,i}   \\
    //!    0    & \ldots  & \ldots & \ldots &    0      & h_{i+1,i} \\
    //! \end{array}\right)
    //! \in \mathcal{R}^{(i+1) \times i }
    //! \f]
    //! In the QR factorisation part of the algorithm this matrix is
    //! over-written (partly) by the R part of its factorisation. Since
    //! both H and R are built up incrementally, we have at the end of the
    //! i-th step of the inner SYMMLQ(m) loop
    //! \f[
    //! \tilde{H} = \left( \begin{array}{cccccc}
    //! h_{1,1} & h_{1,2} & \ldots & \ldots & h_{1,i-1} & h_{i,i}   \\
    //! h_{2,1} & r_{2,2} & \ldots & \ldots & r_{2,i-1} & r_{2,i}   \\
    //!    0    & h_{3,2} & \ddots &        &  \vdots   & \vdots    \\
    //!  \vdots & \ddots  & \ddots & \ddots &  \vdots   & \vdots    \\
    //!  \vdots &         & \ddots & \ddots & r_{i,i-1} & r_{i,i}   \\
    //!  \vdots &         &        & \ddots & h_{i,i-1} & r_{i,i}   \\
    //!    0    & \ldots  & \ldots & \ldots &    0      & h_{i+1,i} \\
    //! \end{array}\right)
    //! \f]
    //! Note that \f$h_{1,k} = r_{1,k}\f$ for all indices k.
    //! Technically we use a two-level data layout with an arrays of pointers
    //! to arrays of decreasing length to represent H. The maximal size of
    //! this data structure can be determined a priori from the maximal
    //! dimension of the Krylov subspace. The data layout is such that
    //! hMat_[i][k] stores \f$h_{i,k}\f$.
    T** hMat_;

    //! Number of rows of upper Hessenberg matrix

    //! This variable is used to keep track on the length of the management
    //! array in the internal hMat_ data structures and corresponds to the
    //! number of rows of the upper Hessenberg matrix, i.e. it is equal to
    //! the maximal dimension of the Krylov subspace plus one. We use a
    //! separate variable to store this for safety reasons, since maxKrylovDim_
    //! is re-read from myParams_ each time Setup() or Solve() is called.
    UInt hMatNumRows_;

    //! Data structure to store the basis vectors spanning the Krylov subspace

    //! This array is used to store the basis vectors \f$v_j, j=1,\ldots k\f$
    //! computed via the modified Gram-Schmidt procedure, which span the
    //! Krylov subspace \f$K^k(A;r_0)\in\mathcal{K}^n\f$. Here
    //! \f$\mathcal{K} = \mathcal{R}\f$ or \f$\mathcal{C}\f$ depending on the
    //! problem. The array has enough space to store maxKrylovDim basis vectors
    //! and vMat_[i] is the i-th basis vector.
    Vector<T>* vMat_;

    //! Maximal number of basis vectors

    //! This variable is used to keep track on the length of the management
    //! array in the internal vMat_ data structures and corresponds to the
    //! maximal number of basis vectors for spanning the Krylov subspace, i.e.
    //! maxKrylovDim. We use a separate variable to store this for safety
    //! reasons, since maxKrylovDim_ is re-read from myParams_ each time
    //! Setup() or Solve() is called.
    UInt vMatNumCols_;

    //! Object for computing the coefficients of the Givens rotation
    GivensRotation *givens_;

    //! Core method for performing the SYMMLQ iteration

    //! This method implements the core of the SYMMLQ iteration. This core
    //! consists of generating the basis of the Krylov subspace up to the
    //! maximal allowed dimension (if required), determination of the
    //! associated upper Hessenberg matrix and solution of the related
    //! least-squares problem.\n\n
    //! This method is called by Solve() to perform the inner loop of the
    //! SYMMLQ(m) iteration, while Solve() itself handles the outer loop and
    //! the computation of the approximate solutions from the intermediate
    //! information.
    void InnerLoop( const BaseMatrix &sysMat, const BasePrecond &precond,
		    const BaseVector &rhs, BaseVector &sol, Double beta,
		    Double threshold, bool &approxIsGood, UInt &stepCount,
		    UInt globNum );

    //! Check if call to Setup() is required.

    //! Auxillary method for testing, whether the internal data structures
    //! must be re-sized. If this is the case a call to Setup() must be
    //! performed. The method will adapt the values of problemDim_ and
    //! maxKrylovDim_, if the respective dimensions have changed.
    bool MustPerformSetup( const BaseMatrix &sysMat, bool &newKdim,
			   bool &newPDim );

    //! Output usage information

    //! Currently the SYMMLQSolver class only supports problem matrices that
    //! are instances of the BNMatrix class and have scalar entries of type
    //! Double or Complex. This method is used to test the type of the input
    //! matrix for correctness. If the matrix fails the test the method will
    //! issue a corresponding error message.
    void TestMatrixType( const BaseMatrix &sysMat ) const;

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    SYMMLQSolver(){};

    //! Prepare internal data structures

    //! The implementation of the Setup method for the SYMMLQ solver performs
    //! mainly memory management tasks. It will prepare the internal data
    //! structures. If this is the first call to Setup or if maxKrylovDim
    //! has changed, it will allocate new memory, freeing previously allocated
    //! memory if necessary.
    //! \note This version of Setup() has an interface that differes from that
    //! of the inherited/over-written Setup() by a const. It is thus callable
    //! by the Solve() method. This is the actual implementation of Setup()
    //! for this class. The public Setup() falls back on this.
    void PrivateSetup( const BaseMatrix &sysMat );

  };

}

#endif
