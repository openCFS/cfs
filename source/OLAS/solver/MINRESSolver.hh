#ifndef OLAS_MINRES_HH
#define OLAS_MINRES_HH



#include "OLAS/utils/math/GivensRotation.hh"
#include "OLAS/utils/math/LanczosMethod.hh"

#include "BaseSolver.hh"

namespace CoupledField {


  //! Implementation of the Minimal Residual (MINRES) algorithm

  //! This class implements the Minimal Residual (MINRES) algorithm for the
  //! iterative solution of a linear system of equations \f$Ax=b\f$
  //! with a sparse and symmetric/Hermitean matrix
  //! \f$ A \in \mathcal{C}^{(n\times n)}\f$,
  //! which need not be positive definite.
  //! Starting from an initial guess \f$x_0\f$ the MINRES method iteratively
  //! builds an orthonormal (unitary) basis of the Krylov subspace
  //! \f[
  //! K^q(A;r_0) = \{ r_0, Ar_0, A^2 r_0, A^3r_0, \ldots, A^{q-1}r_0\}
  //! \enspace.
  //! \f]
  //! In each iteration step a new approximate solution
  //! \f[
  //! x_k\in x_0 + K^{k-1}(A;r_0)
  //! \f]
  //! is determined such that the Euclidean norm of the new residual
  //! \f$\|b-Ax_k\|_2\f$ is minimal. If we denote by \f$Q_k\f$ the matrix
  //! consisting of the first k basis vectors then the minimisation problem
  //! can be transformed into
  //! \f[
  //! \min_{x_k}\|b-Ax_k\|_2 = \min_{v\in\mathcal{K}^n}
  //! \|b - A(x_0 + Q_k v)\|_2 = \min_{v\in\mathcal{K}^k}\|r_0 - AQ_k v\|_2
  //! \enspace.
  //!\f]
  //! Now one can use the following two facts to further transform the
  //! minimisation problem. On the one hand we have \f$A Q_k = Q_{k+1} T_k\f$
  //! where, due to the symmetry of \f$A\f$, \f$T_k\f$ is a tri-diagonal
  //! matrix
  //! \f[
  //! T_k = \left(
  //! \begin{array}{cccccc}
  //! \alpha_1 & \beta_1  &    0    & \ldots    &  \ldots     &      0     \\
  //! \beta_1  & \alpha_2 & \beta_2 &           &             &  \vdots    \\
  //!     0    & \ddots   & \ddots  & \ddots    &             &  \vdots    \\
  //! \vdots   &          & \ddots  & \ddots    &  \ddots     &      0     \\
  //! \vdots   &          &         & \ddots    &  \ddots     & \beta_{k-1}\\
  //! \vdots   &          &         &           & \beta_{k-1} &  \alpha_k  \\
  //!     0    & \ldots   & \ldots  & \ldots    &     0       &  \beta_k   \\
  //! \end{array}\right) \in\mathcal{K}^{(k+1)\times k}
  //! \enspace.
  //! \f]
  //! On the other hand we choose the scaled initial residual as the first
  //! basis vector, i.e. \f$q_1 = r_0/\|r_0\|_2\f$. Thus, the minimisation
  //! problem can be further transformed into
  //! \f[
  //! \min_{x_k}\|b-Ax_k\|_2 =
  //! \min_{v\in\mathcal{K}^k}\| \varrho Q_{k+1} e_1 - Q_{k+1} T_k v\|_2 =
  //! \min_{v\in\mathcal{K}^k}\| \varrho e_1 - T_k v\|_2
  //! \enspace,
  //!\f]
  //! where \f$\varrho=\|r_0\|_2\f$.
  //! In MINRES this least-squares problem is solved by computing a QR
  //! decomposition of \f$T_k\f$ using Givens rotations. This yields
  //! \f$T_k = G_k R_k = \hat{G}_k \hat{R}_k\f$ with
  //! \f$R_k\in\mathcal{K}^{(k+1)\times k}\f$ an upper diagonal matrix of
  //! bandwidth 2 and \f$G_k\f$ coming from the Givens rotations. The matrices
  //! with a hat on top denote the "thin" version of the QR decomposition.
  //! We now obtain
  //! \f[
  //! \min_{x_k}\|b-Ax_k\|_2 =
  //! \min_{v\in\mathcal{K}^k}\| \varrho e_1 - G_k R_k v\|_2
  //! \enspace.
  //! \f]
  //! This least-squares problem can be solved by solving the associated
  //! normal equations
  //! \f[
  //! G_k^H G_k R_k v = R_k v = G_k^H \varrho e_1
  //! \enspace.
  //! \f]
  //! Thus, we can write the next iterate as
  //! \f[
  //! x_k = x_0 + Q_k v = x_0 + Q_k \hat{R}^{-1}_k \hat{G}^H_k \varrho e_1
  //!     = x_0 + W_k \left( \hat{G}^H_k \varrho e_1 \right)
  //! \enspace.
  //! \f]
  //! A little more detailed analysis shows that
  //! - due to the symmetry of the matrix \f$A\f$ the basis vectors can be
  //!   computed by a 3-term recurrence (Lanczos relation)
  //! - the QR decomposition can, with the help of the Givens rotations,
  //!   be computed incrementally
  //! - the approximates can be computed as \f$x_k = x_{k-1} + u_k\f$
  //!   where the update vectors \f$u_k\f$ can also be computed by a 3-term
  //!   recurrence relation. The special structure of \f$W_k\f$ yields
  //!   \f$w_k =\frac{1}{l_0}\left( q_k - l_2 w_{k-2} - l_1 w_{k-1}\right)\f$.
  //!   Here \f$l_2, l_1 \mbox{ and } l_0\f$ denote the last three non-zero
  //!   entries in the final column of \f$\hat{R}_k\f$.
  //!
  //! \par References:
  //!
  //! For further details see e.g.
  //! -# <em>Henk van der Vorst, "Iterative Krylov Methods for Large Linear
  //!    Systems", Cambridge Monographs on Applied and Computational
  //!    Mathematics</em>\n
  //! -# <em>Gerard L. G. Sleijpen, Henk A. van der Vorst and Jan Modersitzki,
  //!    "Differences in the effects of rounding errors in Krylov solvers for
  //!    symmetric indefinite linear systems", SIAM Journal on Matrix Analysis
  //!    and Applications, 2000, 22(3), 726-751</em>
  //!
  //! \par Preconditioning:
  //!
  //! Preconditioning iterative solvers for indefinite symmetric problems is
  //! a difficult task. By default, the left- or right-preconditioned linear
  //! system does not have a symmetric problem matrix. On the other hand the
  //! typical preconditioning approaches do not yield a symmetric positive
  //! definite preconditioner, if \f$A\f$ is symmetric, but indefinite. Thus,
  //! it is not possible to simply switch over to another inner product as
  //! one does in CG.\n
  //! This implementation currently supports left-preconditioning, i.e.
  //! we consider the linear system \f$M^{-1}Ax=M^{-1}b\f$. This is limited to
  //! preconditioners for which \f$M^{-1}A=AM^{-T}\f$, e.g. weighted Jaboci or
  //! polynomial preconditioners.
  //! An alternative would be to apply preconditioning from both sides, i.e.
  //! to consider \f$M^{-1}AM^{-T}M^Tx=M^{-1}b\f$. However, this concept is
  //! not implemented.
  //!
  //! \par Steering parameters and report:
  //!
  //! The behaviour of the MINRESSolver can be tuned via the following
  //! parameters which are read from the parameter object myParams_:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for MINRES solver</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>MINRES_maxIter</td>
  //!       <td align="center">&gt; = 1</td>
  //!       <td align="center">50</td>
  //!       <td>Maximal number of iteration steps allowed for the
  //!         MINRES method.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>MINRES_epsilon</td>
  //!       <td align="center">&gt; 0</td>
  //!       <td align="center">1e-6</td>
  //!       <td>The value of this parameter is used for the convergence
  //!         criterion in the stopping test of MINRES. The stopping test
  //!         is as follows
  //!         \f[
  //!         \|r^{(k)}\|_2  \leq \mbox{tol} := \varepsilon \cdot \tau
  //!         \f]
  //!         where tol is computed by #ComputeThreshold()
  //!         depending on the chosen stopping rule.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>MINRES_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>This parameter is used to control the verbosity of the method.
  //!           If it is set to 'yes', then information on the MINRES
  //!           convergence, like e.g. the norm of the residual per iteration
  //!           step will be logged to the standard %OLAS report stream
  //!           (*cla) ---removed logging---.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //! \n\n The above parameters are read from the parameter object myParams_
  //! every time a solve of the linear system is triggered. Once the system
  //! has been solved, the MINRESSolver writes two values into the report
  //! object myReport_\n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="2" align="center">
  //!         <b>Values the MINRES solver writes into myReport_</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Report</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>numIter</td>
  //!       <td>This is the number of iterations steps the MINRES method has
  //!         performed.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>finalNorm</td>
  //!       <td>This is the final residual norm \f$\|r^{(k)}\|_2\f$.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \attention\n
  //! The MINRESSolver class formally supports as problem matrices objects
  //! of type BaseMatrix. However, the algorithm will only work for those
  //! children of BaseMatrix which have scalar entries of type Double or
  //! Complex. This is tested by the current implementation. Note also that
  //! the vector class associated to the respective matrix class must overload
  //! the assignment operator = for the implementation to work.
  //! See TestMatrixType() for further details.
  //!
  //! \attention\n
  //! The MINRES method is susceptible to the phenomenon of <b>false
  //! convergence</b>. The meaning of this is the following. The residual norm
  //! for the stopping test is not computed explicitely in each step, but is
  //! determined from the minimum of the least-squares problem which is
  //! readily available in the algorithm. It is now possible that the actual
  //! residual is much worse than this value. For examples of this and more
  //! details on the influence of rounding errors on MINRES see e.g.
  //! reference 2.
  template<typename T> class MINRESSolver : public BaseIterativeSolver {

  public:

    //! Constructor

    //! This constructor initialises the pointers to the communication
    //! objects.
    //! \param myParams pointer to parameter object with steering parameters
    //!                 for this solver
    //! \param myReport pointer to report object for storing general
    //!                 information on solution process
    MINRESSolver( PtrParamNode solverNode, PtrParamNode olasInfo );

    //! Default Destructor

    //! The default destructor handles de-allocation of memory of the
    //! internal data structures.
    ~MINRESSolver();

    //! Prepare internal data structures

    //! The implementation of the Setup method for the MINRES solver performs
    //! mainly memory management tasks. It will take care of (re-)allocation
    //! of the auxilliary vectors.
    void Setup( BaseMatrix &sysMat );

    //! Dummy method: Notify the solver that a new matrix pattern has been set
    void SetNewMatrixPattern() {EXCEPTION("SetNewMatrixPattern not implemented!");};

    //! Solve a linear system using MINRES.

    //! Solve the linear system \f$Ax=b\f$ with the square matrix \f$A\f$ and
    //! right-hand side \f$b\f$ for \f$x\f$ using right preconditioned
    //! MINRES(m).
    //! \param sysMat  system matrix \f$A\f$
    //! \param rhs     right-hand side vector \f$b\f$
    //! \param sol     on input initial guess for the solution \f$x\f$, on
    //!                exit approximate solution
    void Solve( const BaseMatrix &sysMat,
                const BaseVector &rhs, BaseVector &sol );

    //! Query type of this solver.

    //! This method can be used to query the type of this solver. The answer
    //! is encoded as a value of the enumeration data type SolverType.
    //! \return MINRES
    SolverType GetSolverType() {
      return MINRES;
    }

  private:

    //@{
    //! Auxillary vector for 3-term recurrence of update vectors.

    //! This vector is used in the 3-term recurrence for computing the
    //! update vector \f$w_k\f$. In the k-th step, i.e. when we compute
    //! the k-th approximate solution \f$x_k\f$, we have the following
    //! relationships w2_ = \f$w_{k-2}\f$, w1_ = \f$w_{k-1}\f$ and
    //! w0_ = \f$w_k\f$
    BaseVector *w2_, *w1_, *w0_;
    //@}

    //@{
    //! Auxillary vector for 3-term recurrence of basis vectors.

    //! This vector is used in the 3-term recurrence for computing the
    //! basis vector \f$q_{k+1}\f$. In the k-th step, i.e. when we compute
    //! the k-th approximate solution \f$x_k\f$, we have the following
    //! relationships q1_ = \f$q_{k-1}\f$, q0_ = \f$q_k\f$ and
    //! qN_ = \f$q_{k+1}\f$
    BaseVector *q1_, *q0_, *qN_;
    //@}

    //! Auxillary vector for application of preconditioner

    //! This vector is used as an intermediate vector in the application of
    //! the preconditioner.
    BaseVector *pV_;

    //! Right-hand side vector of least-squares problem

    //! This auxillary vector is used in the solution of the least squares
    //! problem to store the right-hand side of the problem. It is defined
    //! as
    //! \f[ b = \|r^{(0)}\|_2 e_1 \f]
    //! where \f$e_1\f$ is the first canonical basis vector and \f$r^{(0)}\f$
    //! is the residual of the initial guess.
    BaseVector *bV_;

    //! Object for computing the coefficients of the Givens rotation
    GivensRotation *givens_;

    //! Object for computing the Lanczos vectors
    LanczosMethod *lanczos_;

    //! Length of auxilliary vectors;

    //! This attribute is used to keep track of the length of the
    //! internal vector objects used in the 3-term recurrences. This length
    //! is identical to the dimension of the problem or larger.
    UInt vectorLength_;

    //! Default Constructor

    //! The default constructor is private in order to disallow its use. The
    //! reason is that we need pointers to the parameter and report objects to
    //! perform correct setup of an object of this class.
    MINRESSolver(){};

    //! Prepare internal data structures

    //! The implementation of the Setup method for the MINRES solver performs
    //! mainly memory management tasks. It will take care of (re-)allocation
    //! of the auxilliary vectors.
    //! \note This version of Setup() has an interface that differes from that
    //! of the inherited/over-written Setup() by a const. It is thus callable
    //! by the Solve() method. This is the actual implementation of Setup()
    //! for this class. The public Setup() falls back on this.
    void PrivateSetup( const BaseMatrix &sysMat );

    //! Test suitability of matrix type for MINRES solver

    //! Formally the MINRESSolver class supports all children of the
    //! BaseMatrix class as problem matrices. There is a caveat, however.
    //! The entries of the matrix on the scalar level must be Double or
    //! Complex and are not allowed to be tiny matrices. The latter would in
    //! principal also be possible, but would require us to have a possibility
    //! to generate for such a matrix a corresponding vector with scalar
    //! entries of dimension one. This is no principal problem, but we do not
    //! yet support it.
    void TestMatrixType( const BaseMatrix &sysMat ) const;

  };

}

#endif
