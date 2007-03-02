#ifndef BASE_LANCZOS_HH
#define BASE_LANCZOS_HH

namespace OLAS {

  //! For every symmetric (Hermitean) matrix \f$A\in \mathcal{K}^{n\times n}\f$
  //! (\f$\mathcal{K} = \mathcal{R}\mbox{ or }\mathcal{C}\f$) one can construct
  //! an orthonormal (unitary) matrix \f$Q\f$ such that \f$Q^TAQ = T\f$ with
  //! a tridiagonal matrix
  //! \f[
  //! T = \left(
  //! \begin{array}{cccccc}
  //! \alpha_1 & \beta_1  &    0    &   \ldots    &  \ldots     &      0     \\
  //! \beta_1  & \alpha_2 & \beta_2 &             &             &  \vdots    \\
  //!     0    & \ddots   & \ddots  &   \ddots    &             &  \vdots    \\
  //! \vdots   &          & \ddots  &   \ddots    &  \ddots     &      0     \\
  //! \vdots   &          &         &   \ddots    &  \ddots     & \beta_{n-1}\\
  //!     0    & \ldots   & \ldots  &      0      & \beta_{n-1} &  \alpha_n  \\
  //! \end{array}\right)
  //! \enspace.
  //! \f]
  //! From the identity \f$AQ = QT\f$ one obtains for the columns \f$q_k\f$
  //! of \f$Q\f$
  //! \f[
  //! \beta_k q_{k+1} = A q_k - \beta_{k-1} q_{k-1} - \alpha_k q_k
  //! \enspace,
  //! \f]
  //! with \f$\beta_0q_0 = 0\f$. When \f$r_k\f$ is defined as the right hand
  //! side of the above equation, then the entries of \f$T\f$ are given by
  //! \f$\alpha_k = q_k^\ast A q_k\f$ and \f$\beta_k = \pm\|r_k\|_2\f$.
  //! Without loss of generality one can choose the \f$\beta_k\f$ to be
  //! positive.
  //! The Lanczos algorithm is a method to compute the basis \f$Q\f$ and the
  //! related tridiagonal matrix \f$T\f$ based on the 3-term recurrence
  //! defined above. The vectors \f$q_k\f$ are often denotes as Lanczos
  //! vectors. The Lanczos algorithm is the core
  //! part of several other methods for solving linear system and computing
  //! eigenvalues, like e.g. MINRES or SYMMLQ. The algorithm is not uniquely
  //! defined, however, and there exist different flavours of implementing it.
  //! \n\n
  //! This class implements one version of the Lanczos algorithm in order to
  //! provide a common implementation for the different classes like
  //! MINRESSolver or SYMMLQSolver which rely on the Lanczos process.
  class LanczosMethod {

    public:

    //@{
    //! Compute next orthonormal basis vector via 3-term recurrence

    //! This method computes the next basis vector, i.e. another column
    //! vector of the matrix \f$Q\f$ from the previous ones based on a 3-term
    //! recurrence. This implementation is based on the observation that
    //! \f$\alpha_k = q_k^\ast\left(Aq_k - \beta_{k-1}q_{k-1}\right)\f$ and
    //! proceeds in the following steps
    //! \f[
    //! \begin{array}{r@{\:\longleftarrow\:}l}
    //! \displaystyle q_3      & \displaystyle Aq_2 - \beta_1 q_1  \\[2ex]
    //! \displaystyle \alpha_2 & \displaystyle q_2^\ast q_3        \\[2ex]
    //! \displaystyle r_2      & \displaystyle q_3 - \alpha_2 q_2  \\[2ex]
    //! \displaystyle \beta_2  & \displaystyle \|r_2\|_2           \\[2ex]
    //! \displaystyle q_3      & \displaystyle \frac{r_2}{\beta_2} \\
    //! \end{array}
    //! \f]
    //! \note If a preconditioner is supplied, then we replace the matrix
    //!       \f$A\f$ in the above formulas by \f$\tilde{A}=M^{-1}A\f$.
    //! \param mat     (input)  matrix \f$A\f$
    //! \param precond (input)  preconditioner for \f$A\f$
    //! \param q1      (input)  last but one basis vector
    //! \param q2      (input)  last basis vector
    //! \param q3      (output) next basis vector
    //! \param aux     (input)  auxilliary vector for applying preconditioner
    //! \param alpha2  (output) \f$\alpha_2 = q_2^\ast Aq_2\f$
    //! \param beta1   (input)  \f$\|r_1\|_2\f$
    //! \param beta2   (output) \f$\|r_2\|_2\f$
    void CompNextVector( const BaseMatrix &mat, const BasePrecond &precond,
			 const BaseVector &q1, const BaseVector &q2,
			 BaseVector &q3, BaseVector &aux, Double &alpha2,
			 const Double &beta1, Double &beta2 ) {

      ENTER_IFCN( "LanczosMethod::CompNextVector" );

      mat.Mult( q2, aux );
      precond.Apply( mat, aux, q3 );
      q3.Add( - beta1, q1 );
      q2.Inner( q3, alpha2 );
      q3.Add( - alpha2, q2 );
      beta2 = q3.NormEuclid();
      q3.ScalarDiv( beta2 );

    };

    void CompNextVector( const BaseMatrix &mat, const BasePrecond &precond,
			 const BaseVector &q1, const BaseVector &q2,
			 BaseVector &q3, BaseVector &aux, Complex &alpha2,
			 const Complex &beta1, Complex &beta2 ) {

      ENTER_IFCN( "LanczosMethod::CompNextVector" );

      mat.Mult( q2, aux );
      precond.Apply( mat, aux, q3 );
      q3.Add( - beta1, q1 );
      q2.Inner( q3, alpha2 );
      q3.Add( - alpha2, q2 );
      beta2 = q3.NormEuclid();
      q3.ScalarDiv( beta2 );

    };

    //@}

  };

}

#endif
