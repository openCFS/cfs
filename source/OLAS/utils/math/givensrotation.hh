// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef GIVENS_ROTATION_HH
#define GIVENS_ROTATION_HH

#include "utils/defs.hh"
#include "utils/environment.hh"

namespace OLAS {

  //! Class for performing a Givens rotation.

  //! This class offers methods for computing the coefficients of a real or
  //! complex valued Givens rotation and performing that rotation. A Givens
  //! rotation matrix \f$G(c,s)\f$ is defined as the matrix that when
  //! multiplied with a prescribed vector zeros its second component:
  //! \f[
  //! G(c,s) \left(\begin{array}{c} f \\ g \end{array}\right) =
  //! \left(\begin{array}{rr} c & s \\ -\bar{s} & \bar{c} \end{array}\right)
  //! \left(\begin{array}{c} f \\ g \end{array}\right) =
  //! \left(\begin{array}{c} r \\ 0 \end{array}\right)
  //! \f]
  //! under the side-constraint that \f$G(c,s)\f$ is a unitary transformation,
  //! i.e. \f$G(c,s)^H G(c,s) = \mbox{Id}\f$. The coefficients \f$c\f$ and
  //! \f$s\f$ that solve this problem are not uniquely determined and there
  //! exist different flavours for computing a solution pair (c,s). This is
  //! especially true for the complex case. Take e.g. the case that \f$f\f$ and
  //! \f$g\f$ are real and positive. Then the widely accepted convention is to
  //! set
  //! \f[
  //! \begin{array}{rcl}
  //! c &=& \displaystyle\frac{f}{\sqrt{f^2 + g^2}}\\[4ex]
  //! s &=& \displaystyle\frac{g}{\sqrt{f^2 + g^2}}\\[4ex]
  //! r &=& \displaystyle\sqrt{f^2 + g^2}\enspace.
  //! \end{array}
  //! \f]
  //! However, the negatives of \f$c\f$, \f$s\f$ and \f$r\f$ also solve the
  //! problem. This class offers two possibilities for computing the triple
  //! (c,s,r)
  //! -# By setting interface_ to LAPACK one can use an implementation of the
  //!    LAPACK library for computing the Givens rotation. In this case the
  //!    fitting one of the four LAPACK routines SLARTG, DLARTG, CLARTG or
  //!    ZLARTG will be used for computing (c,s,r). The details of the style
  //!    of computing (c,s,r) then depend on the way the given LAPACK library
  //!    implements the corresponding routine.\n\n
  //! -# By setting interface_ to %OLAS one can use our private implementation
  //!    of the Givens rotation. This implementation is based on reference [1]
  //!    below. In the complex case it implements <b>Algorithm 3</b> from [1]
  //!    omitting, however, the scaling operations introduced in [1] to
  //!    guarantee numerical stability. In the real case it implements
  //!    <b>Algorithm 9</b> from [1]. This is a fast version without scaling.
  //!    The implementation also currently does not implement any explicit
  //!    treatment of exceptions. The behaviour for input that is NaN or
  //!    \f$\pm\f$Inf is not well defined.
  //!
  //! Both versions, however, adhere to the standard convention that the
  //! coefficient \f$c\f$ should be a real number.\n\n
  //! For more details on the different flavours of Givens rotation and a
  //! description of the variant of computation proposed in the new standard
  //! of the Basic Linear Algebra Subroutines Technical Forum (BLAST) see:\n\n
  //! [1] David Bindel, James Demmel, William Kahan and Osni Marques,
  //!   <b>On Computing Givens Rotations Reliably and Efficiently</b>,
  //!   ACM Transactions on Mathematical Software, Vol. 28, No. 2, June 2002,
  //!   pages 206 - 238
  class GivensRotation {

  public:

    //! Enumeration data type for distinguishing the different approaches
    //! for computing the Givens rotation. Possible values are:
    //! - LAPACK
    //! - %OLAS
    typedef enum { LAPACK, OLAS } CompStyle;

    //! Constructor
    GivensRotation( GivensRotation::CompStyle interface );

    //! Default destructor
    ~GivensRotation();

    //! Givens rotation for real valued case
    void gRot( Double f, Double g, Double &c, Double &s, Double &r ) const;

    //! Givens rotation for complex valued case
    void gRot( Complex f, Complex g, Double &c, Complex &s, Complex &r ) const;

  private:

    //! Attribute for setting the type of computational interface
    GivensRotation::CompStyle interface_;

    //! Default constructor

    //! The default constructor is dis-allowed, since we want the user to
    //! decide on the computational interface that should be used for the
    //! computation of the Givens rotation already when generating an
    //! instance of this class.
    GivensRotation();

    //! Our implementation of the real Givens rotation

    //! This basically implements "Algorithm 3" of reference [1], however,
    //! without the scaling to ensure numerical stability, yet.
    void CompRot( Double f, Double g, Double &c, Double &s, Double &r ) const;

    //! Our implementation of the complex Givens rotation

    //! This implements the fast "Algorithm 9" for computing a real Givens
    //! rotation of reference [1]. Note that this algorithm does not employ
    //! scaling to ensure numerical stability.
    void CompRot( Complex f, Complex g, Double &c, Complex &s,
		  Complex &r ) const;

  };

}

#endif
