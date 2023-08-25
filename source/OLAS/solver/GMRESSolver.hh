#ifndef OLAS_GMRES_HH
#define OLAS_GMRES_HH

#include <vector>



#include "OLAS/utils/math/GivensRotation.hh"


#include "BaseSolver.hh"

namespace CoupledField {

  class BasePrecond;

  // Forward declaration of classes
  // class GivenRotation;

  //! Implementation of the Generalised Minimal Residual Algorithm

  //! This class implements the Generalised Minimal Residual Algorithm (GMRES)
  //! for the iterative solution of a linear system of equations \f$Ax=b\f$
  //! with a sparse and not necessarily symmetric/Hermitean matrix
  //! \f[
  //! A \in \mathcal{C}^{(n\times n)} \enspace.
  //! \f]
  //! Starting from an initial guess \f$x_0\f$ the GMRES method iteratively
  //! builds a basis of the Krylov subspace
  //! \f[
  //! K^q(A;r_0) = \{ r_0, Ar_0, A^2 r_0, A^3r_0, \ldots, A^{q-1}r_0\}
  //! \enspace.
  //! \f]
  //! In each iteration step a new approximate solution
  //! \f[
  //! x_k\in x_0 + K^{k-1}(A;r_0)
  //! \f]
  //! is determined such that the Euclidean norm of the new residual
  //! \f$\|b-Ax_k\|_2\f$ is minimal. Since each new iteration step is more
  //! costly than the previous one, due to need to orthogonalise the new
  //! basis vector of the Krylov subspace against all previous ones, the
  //! GMRES method is often used in a re-started version. This is denoted
  //! as GMRES(m) and means that after computing m new approximations and
  //! reaching a maxmimal dimension (m+1) of the Krylov subspace, the subspace
  //! information is thrown away and the current approximation used as new
  //! initial guess for the next step of the GMRES(m) method.\n\n
  //! This class represents a straightforward implementation of GMRES. It uses
  //! the modified Gram-Schmidt Orthogonalisation procedure for the generation
  //! of the basis of the underlying Krylov subspace and employs successive
  //! Givens rotations for solving the resulting least-squares problem.
  //! For further details
  //! see <b>Henk van der Vorst, "Iterative Krylov Methods for Large Linear
  //! Systems", Cambridge Monographs on Applied and Computational
  //! Mathematics</b>.\n\n
  //! The GMRES method can be used in combination with a preconditioner. This
  //! class implements so called right-preconditioning. In this approach the
  //! linear system is replaced by an equivalent one in the following fashion
  //! \f[
  //! Ax = b \longleftrightarrow A M^{-1} Mx = b
  //! \f]
  //! where \f$M\f$ is the matrix representing the preconditioner.\n\n
  //! The behaviour of the %GMRESSolver can be tuned via the following
  //! parameters which are read from the parameter object myParams_:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for GMRES solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>GMRES_maxKrylovDim</td>
  //!       <td align="center">&gt; = 1</td>
  //!       <td align="center">50</td>
  //!       <td>This gives the upper limit on the dimension of the Krylov
  //!         subspace generated in each step of the GMRES(m) method before
  //!         the iteration is re-started. Thus, it is equal to (m-1), since
  //!         the residual of the initial guess constitutes the first basis
  //!         vector we can compute (m-1) new approximations before the
  //!         maximal dimension is attained.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>GMRES_maxIter</td>
  //!       <td align="center">&gt; = 1</td>
  //!       <td align="center">1</td>
  //!       <td>Maximal number of iteration steps allowed for the
  //!         GMRES(m) method. If this is set to 1, then we are
  //!         actually performing a full GMRES iteration and the
  //!         maximal number of iterations is determined by
  //!         the value of GMRES_maxKrylovDim.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>GMRES_epsilon</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">1e-6</td>
  //!       <td>The value of this parameter is used for the convergence
  //!         criterion in the stopping test of GMRES. The stopping rule
  //!         is as follows
  //!         \f[
  //!         \|r^{(k)}\|_2 \leq \varepsilon \cdot \|b\|_2
  //!         \f]
  //!         where \f$r^{(k)}\f$ is the residual of the k-th step
  //!         and \f$b\f$ is the right hand side. If \f$\|b\|_2=0\f$,
  //!         or in the case that the penalty formulation is used to treat
  //!         inhomogeneous Dirichlet boundary conditions, then \f$\|b\|_2\f$
  //!         is replaced by \f$\|r^{(0)}\|_2\f$.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>GMRES_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>This parameter is used to control the verbosity of the method.
  //!           If it is set to 'yes', then information on the GMRES
  //!           convergence, like e.g. the norm of the residual per iteration
  //!           step will be logged to the standard %OLAS report stream
  //!           (*cla) ---removed logging---.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \n\n The above parameters are read from the parameter object myParams_
  //! every time a solve of the linear system is triggered. Once the system
  //! has been solved, the %GMRESSolver writes four values into the report
  //! object myReport_\n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="2" align="center">
  //!         <b>Values the GMRES solver writes into myReport_</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Report</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>numIter</td>
  //!       <td>This is the number of iterations steps the GMRES method has
  //!         performed. For GMRES(m) this is the number of outer iterations,
  //!         while for full GMRES it is simply the total number of
  //!         iterations.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>numIterGlobal</td>
  //!       <td>This is the global number of iterations steps, i.e. the total
  //!           number of approximations computed over all inner iterations of
  //!           the GMRES(m) method. For full GMRES this is identical to
  //!           numIter.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>finalNorm</td>
  //!       <td>This is the Euclidean norm of the final residual:
  //!         \f$\displaystyle \left\|r^{(k)}\right\|_2 =
  //!         \left\|b - Ax^{(k)}\right\|_2\enspace.\f$
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>solutionIsOkay</td>
  //!       <td>This boolean is set to true, if the iteration was terminated
  //!           because the stopping criterion was fulfilled and to false
  //!           otherwise (e.g. if false convergence was detected).
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \note
  //! - The GMRES method is susceptible to the phenomenon of <b>false
  //!   convergence</b>. The meaning of this is the following. The residual
  //!   norm used for the stopping test is not computed explicitely in each
  //!   step of the full GMRES method or in each inner iteration of the
  //!   GMRES(m) method. This would require setting up the approximate
  //!   explicitely, which is quite costly, and involve an additional
  //!   matrix-vector multiplication. Instead the residual norm is obtained
  //!   as \f$\tilde{r} = \min\|b-Ax_k\|_2\f$ which is readily available in
  //!   the algorithm.
  //!   It is now possible that the actual residual is much worse than
  //!   \f$\tilde{r}\f$. For examples of this and more details on the
  //!   influence of rounding errors on GMRES see e.g.\n\n
  //!   <em>G. Sleijpen, H. van der Vorst and J. Modersitzki, Differences in
  //!   the effects of rounding errors in Krylov solvers for symmetric
  //!   indefinite linear systems, SIAM Journal on Matrix Analysis and
  //!   Applications, 2000, 22(3), 726-751.</em>\n\n
  //!   This situation can only be detected after the inner
  //!   loop was terminated. If it occurs, the Solve() method will log this
  //!   to the standard %OLAS report stream. In the case of GMRES(m) Solve()
  //!   will then start with the next outer iteration.
  template<typename T> class GMRESSolver : public BaseIterativeSolver {

  public:

    //! Constructor

    //! This constructor initialises the pointers to the communication
    //! objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    GMRESSolver( PtrParamNode solverNode, PtrParamNode olasInfo );

    //! Destructor

    //! The destructor handles de-allocation of memory of the
    //! internal data structures.
    ~GMRESSolver();

    //! Prepare internal data structures

    //! The implementation of the up method for the GMRES solver performs
    //! mainly memory management tasks. It will prepare the internal data
    //! structures. If this is the first call to Setup or if maxKrylovDim_
    //! has changed, it will allocate new memory, freeing previously allocated
    //! memory if necessary.
    void Setup( BaseMatrix &sysMat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};
    
    //! Set preconditioner object at solver

    //! Solve a linear system using GMRES.

    //! Solve the linear system \f$Ax=b\f$ with the square matrix \f$A\f$ and
    //! right-hand side \f$b\f$ for \f$x\f$ using right preconditioned
    //! GMRES(m).
    //! \param sysMat  system matrix \f$A\f$
    //! \param rhs     right-hand side vector \f$b\f$
    //! \param sol     on input initial guess for the solution \f$x\f$, on
    //!                exit approximate solution
    void Solve( const BaseMatrix &sysMat, 
    const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return GMRES
    SolverType GetSolverType() {
      return GMRES;
    }

  private:

    bool consoleConvergence_;

    //! Right-hand side vector of least-squares problem

    //! This auxilliary vector is used in the solution of the least squares
    //! problem to store the right-hand side of the problem. It is defined
    //! as
    //! \f[ b = \|r^{(0)}\|_2 e_1 \f]
    //! where \f$e_1\f$ is the first canonical basis vector and \f$r^{(0)}\f$
    //! is the residual of the initial guess.
    Vector<T> bVec_;

    //! Auxilliary vector for applying preconditioner

    //! This auxillary vector is used in the application of the
    //! preconditioner. In the right-preconditioned GMRES method we solve a
    //! linear system with a system matrix \f$\hat{A} = A M^{-1}\f$ where
    //! \f$A\f$ is the original system matrix and \f$M^{-1}\f$ is the
    //! preconditioner. This vector is used to store the result of applying
    //! the preconditioner to a vector \f$v\f$, i.e.
    //! \f$\mbox{pVec\_} = M^{-1}v\f$.
    BaseVector *pVec_;

    //! Maximal dimension of the generated Krylov subspace

    //! This integer is used to store the maximal dimension \f$m\f$ of the
    //! Krylov subspace \f$K^m(A;r_0)\f$ generated in the solution of the
    //! linear system. The value is obtained from the GMRES_maxKrylovDim entry
    //! of the myParams_ object each time Setup() or Solve() is called.
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
    //! h_{1,1} & h_{1,2} & \ldots & \ldots &  h_{1,i-1}  &  h_{1,i}  \\
    //! h_{2,1} & h_{2,2} & \ldots & \ldots &  h_{2,i-1}  &  h_{2,i}  \\
    //!    0    & h_{3,2} & \ddots &        &   \vdots    &  \vdots   \\
    //!  \vdots & \ddots  & \ddots & \ddots &   \vdots    &  \vdots   \\
    //!  \vdots &         & \ddots & \ddots & h_{i-1,i-1} & h_{i-1,i} \\
    //!  \vdots &         &        & \ddots &  h_{i,i-1}  & h_{i,i}   \\
    //!    0    & \ldots  & \ldots & \ldots &      0      & h_{i+1,i} \\
    //! \end{array}\right)
    //! \in \mathcal{R}^{(i+1) \times i }
    //! \f]
    //! In the QR factorisation part of the algorithm this matrix is
    //! over-written (partly) by the R part of its factorisation. Since
    //! both H and R are built up incrementally, we have at the end of the
    //! i-th step of the inner GMRES(m) loop
    //! \f[
    //! \tilde{H} = \left( \begin{array}{cccccc}
    //! h_{1,1} & h_{1,2} & \ldots & \ldots &  h_{1,i-1}  & h_{1,i}   \\
    //! h_{2,1} & r_{2,2} & \ldots & \ldots &  r_{2,i-1}  & r_{2,i}   \\
    //!    0    & h_{3,2} & \ddots &        &   \vdots    & \vdots    \\
    //!  \vdots & \ddots  & \ddots & \ddots &   \vdots    & \vdots    \\
    //!  \vdots &         & \ddots & \ddots & r_{i-1,i-1} & r_{i-1,i} \\
    //!  \vdots &         &        & \ddots &  h_{i,i-1}  & r_{i,i}   \\
    //!    0    & \ldots  & \ldots & \ldots &      0      & h_{i+1,i} \\
    //! \end{array}\right)
    //! \f]
    //! Note that \f$h_{1,k} = r_{1,k}\f$ for all indices k.
    //! Technically we use a two-level data layout with an array of pointers
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
    //! separate variable to store this for safety reasons, since
    //! maxKrylovDim_ is re-read from myParams_ each time Setup() or Solve()
    //! is called.
    UInt hMatNumRows_;

    //! Data structure to store the basis vectors spanning the Krylov subspace

    //! This array is used to store the basis vectors \f$v_j, j=1,\ldots k\f$
    //! computed via the modified Gram-Schmidt procedure, which span the
    //! Krylov subspace \f$K^k(A;r_0)\in\mathcal{K}^n\f$. Here
    //! \f$\mathcal{K} = \mathcal{R}\f$ or \f$\mathcal{C}\f$ depending on the
    //! problem. The array has enough space to store maxKrylovDim_ basis
    //! vectors and vMat_[i] is the i-th basis vector.
    std::vector<BaseVector*> vMat_;

    //! Object for computing the coefficients of the Givens rotation
    GivensRotation *givens_;

    //! Core method for performing the GMRES iteration

    //! This method implements the core of the GMRES iteration. This core
    //! consists of generating the basis of the Krylov subspace up to the
    //! maximal allowed dimension (if required), determination of the
    //! associated upper Hessenberg matrix and solution of the related
    //! least-squares problem.\n\n
    //! This method is called by Solve() to perform the inner loop of the
    //! GMRES(m) iteration, while Solve() itself handles the outer loop and
    //! the computation of the approximate solutions from the intermediate
    //! information.
    //! \param sysMat       matrix of the linear system to be solved
    //! \param precond      preconditioner for solution process
    //! \param rhs          right-hand side vector of linear system to be
    //!                     solved
    //! \param sol          on input the initial guess for this inner loop,
    //!                     on output the computed approximate solution
    //! \param beta         Euclidean norm of the residual of the inital guess
    //! \param threshold    threshold for stopping test; should include the
    //!                     scaling factor scalFac (see below) already
    //! \param approxIsGood on return indicates whether the computed
    //!                     approximate satisfies the stopping criterion
    //! \param stepCount    on output the number of iterations performed in
    //!                     the inner loop
    //! \param globNum      on input the number of global iterations steps
    //!                     performed so far
    void InnerLoop( const BaseMatrix &sysMat, const BasePrecond &precond,
        const BaseVector &rhs, BaseVector &sol, Double beta,
        Double threshold, bool &approxIsGood, UInt &stepCount,
        UInt globNum );

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    GMRESSolver() {
      EXCEPTION( "Default constructor of class GMRESSolver is not allowed!" );
    }

    //! Copy Constructor

    //! The copy constructor is private in order to disallow its use. If we
    //! ever need a copy constructor we better implement it ourselves and do
    //! not let the compiler do it.
    GMRESSolver( const GMRESSolver &src ) {
      EXCEPTION( "Copy constructor of class GMRESSolver is not allowed!" );
    }

    //! Prepare internal data structures

    //! The implementation of the Setup method for the GMRES solver performs
    //! mainly memory management tasks. It will prepare the internal data
    //! structures. If this is the first call to Setup or if maxKrylovDim_
    //! has changed, it will allocate new memory, freeing previously allocated
    //! memory if necessary.
    //! \note This version of Setup() has an interface that differs from that
    //! of the inherited/over-written Setup() by a const. It is thus callable
    //! by the Solve() method. This is the actual implementation of Setup()
    //! for this class. The public Setup() falls back on this.
    void PrivateSetup( const BaseMatrix &sysMat );

    //! Allocate memory to store upper Hessenberg matrix

    //! This method can be used to allocate memory and initialise the internal
    //! data structures for storing an upper Hessenberg matrix. If the matrix
    //! hMat_ has already been allocated before, this method will delete the
    //! old matrix and allocate a new one, only if this is necessary, i.e.
    //! if the maximal dimension of the Krylov subspace has increased. Thus,
    //! it is necessary that maxKrylovDim_ is set to the new value before
    //! calling this method!
    void AllocateHessenbergMatrix();

    //! Frees memory allocated to store upper Hessenberg matrix

    //! Calling this method triggers the de-allocation of the internal data
    //! structures allocated by a call to AllocateHessenbergMatrix() in order
    //! to store an upper Hessenberg matrix. It is safe to call this method,
    //! even if no matrix has yet been allocated.
    void DeleteHessenbergMatrix();

    //! Allocate memory to store base vectors of Krylov subspace

    //! This method can be used to allocate memory and initialise the internal
    //! data structures for storing the base vectors of the Krylov subspace.
    //! If a basis had already been allocated before, this method will delete
    //! the old base vectors and allocate new ones, only if this is necessary,
    //! i.e. if the maximal dimension of the Krylov subspace has increased or
    //! if the dimension of the base vectors has changed. Thus, it is
    //! necessary that maxKrylovDim_ and problemDim_ are set to the correct
    //! new values before calling this method!
    //! \note If the system matrix is an SBM_Matrix, then the base vectors
    //! are in fact vectors of the SBM_Vector class. If the order of an
    //! SBM_Matrix/SBM_Vector remains the same, but the internal layout
    //! changes, this will not be detected by this method!
    void AllocateKrylovBasis( const BaseMatrix &sysMat );

    //! Free memory allocated to store base vectors of Krylov subspace

    //! Calling this method triggers de-allocation of the internal data
    //! structures allocated by a call to AllocateKrylovBasis() in order to
    //! store the base vectors of the Krylov subspace. It is safe to call this
    //! method, even if no Krylov basis has yet been allocated.
    void DeleteKrylovBasis();

  };

}

#endif
