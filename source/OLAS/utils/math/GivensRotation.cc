#include <string>
#include <cmath>
#include <complex>

#include "MatVec/BLASLAPACKInterface.hh"
#include "GivensRotation.hh"

namespace CoupledField {


  // ***********************
  //   Default Constructor
  // ***********************
  GivensRotation::GivensRotation() {
    EXCEPTION("This constructor should never be called!");
  }


  // **************
  //   Destructor
  // **************
  GivensRotation::~GivensRotation() {
  }


  // ***************
  //   Constructor
  // ***************
  GivensRotation::GivensRotation( GivensRotation::CompStyle interface ) {
    interface_ = interface;
  }


  // *******************************
  //   Real Givens Rotation Driver
  // *******************************
  void GivensRotation::gRot( Double f, Double g, Double &c, Double &s,
                             Double &r ) const {


    // Determine which approach to use for the Givens rotation
    switch( interface_ ) {

    case LAPACK:

#ifdef USE_LAPACK

      // single precision case
      if ( sizeof(Double) == sizeof(float) ) {
        EXCEPTION("Interface for Double = float no longer supported!");
      }

      // double precision case
      else if ( sizeof(Double) == sizeof(float) ) {
        Double myF = f;
        Double myG = g;
        Double myC = c;
        Double myS = s;
        Double myR = r;
        dlartg( &myF, &myG, &myC, &myS, &myR );
      }

      // something else
      else {
        EXCEPTION("Double is neither float nor double! Don't know what to do!");
      }

#else

      EXCEPTION( "Compile with USE_LAPACK to enable support for "
               << "LAPACK's Givens rotation" );

#endif

      break;

    case OLAS:
      CompRot( f, g, c, s, r );
      break;

    default:
      EXCEPTION("Somethings broke! Should not have reached this line of code!");
    }
  }


  // **********************************
  //   Complex Givens Rotation Driver
  // **********************************
  void GivensRotation::gRot( Complex f, Complex g, Double &c, Complex &s,
                             Complex &r ) const {


    // Determine which approach to use for the Givens rotation
    switch( interface_ ) {

    case LAPACK:

#ifdef USE_LAPACK

      // single precision case
      if ( sizeof(Complex) == sizeof(std::complex<float>) ) {
        EXCEPTION("Interface for Complex = std::complex<float> no longer "
               "supported!");
      }

      // double precision case
      else if ( sizeof(Complex) == sizeof(std::complex<Double>) ) {
        std::complex<Double> myF = f;
        std::complex<Double> myG = g;
        Double     myC = c;
        std::complex<Double> myS = s;
        std::complex<Double> myR = r;
        zlartg( &myF, &myG, &myC, &myS, &myR );
        f = myF;
        g = myG;
        c = myC;
        s = myS;
        r = myR;
      }

      // something else
      else {
        EXCEPTION("Complex is neither complex<float> nor complex<double>! "
                 << "Don't know what to do!");
      }

#else

      EXCEPTION("Compile with USE_LAPACK to enable support for LAPACK's "
               << "Givens rotation");

#endif

      break;

    case OLAS:
      CompRot( f, g, c, s, r );
      break;

    default:
      EXCEPTION("Something's broke! Should not have reached this line of code!");
    }
  }


  // ************************
  //   Real Givens Rotation
  // ************************

  // NOTE: This is based on "Algorithm 9" from the paper
  //       David Bindel, James Demmel, William Kahan and Osni Marques,
  //       On Computing Givens Rotations Reliably and Efficiently,
  //       ACM Transactions on Mathematical Software, Vol. 28, No. 2,
  //       June 2002, pages 206 - 238
  //       "Algorithm 9" includes no exception handling or scaling to avoid
  //!      numerical instability.
  void GivensRotation::CompRot( Double f, Double g, Double &c, Double &s,
                                Double &r ) const {


    // Declaration of local variables
    Double fg2, rr;

    // Simple cases
    if ( f == 0.0 || g == 0.0 ) {
      if ( g == 0.0 ) {
        c = 1.0;
        s = 0.0;
        r = f;
      }
      else if ( f == 0.0 ) {
        c = 0.0;
        s = fabs(g) >= 0.0 ? +1.0 : -1.0;
        r = fabs(g);
      }
    }

    // Algorithm 9: Fast Real Complex Givens Rotations without scaling
    else {
      fg2 = f * f + g * g;
      r   = sqrt( fg2 );
      rr  = 1.0 / r;
      c   = fabs(f) * rr;
      s   = g * rr;
      if ( f < 0.0 ) {
        s = -s;
        r = -r;
      }
    }
  }


  // ***************************
  //   Complex Givens Rotation
  // ***************************

  // NOTE: This is based on "Algorithm 3" from the paper
  //       David Bindel, James Demmel, William Kahan and Osni Marques,
  //       On Computing Givens Rotations Reliably and Efficiently,
  //       ACM Transactions on Mathematical Software, Vol. 28, No. 2,
  //       June 2002, pages 206 - 238
  //       However, we do not currently implement the exception handling
  //       or the scaling to avoid numerical instability.
  void GivensRotation::CompRot( Complex f, Complex g, Double &c, Complex &s,
                                Complex &r ) const {


    // Declaration of local variables
    Double fn, gn, d1, f2, g2, fg2;

    // Compute "norm" of f and g
    fn = fabs(f.real()) > fabs(f.imag()) ? fabs(f.real())
      : fabs(f.imag());
    gn = fabs(g.real()) > fabs(g.imag()) ? fabs(g.real())
      : fabs(g.imag());

    // Simple cases
    if ( fn == 0.0 || gn == 0.0 ) {
      if ( gn == 0.0 ) {
        c = 1.0;
        s = 0.0;
        r = f;
      }
      else if ( fn == 0.0 ) {
        c = 0.0;
        d1 = sqrt( g.real() * g.real() + g.imag() * g.imag() );
        r = d1;
        d1 = 1.0 / d1;
        s = conj(g) * d1;
      }
    }

    // Algorithm 3: Fast Complex Givens Rotations when f and g
    //              are "well scaled"
    else {
      f2  = f.real() * f.real() + f.imag() * f.imag();
      g2  = g.real() * g.real() + g.imag() * g.imag();
      fg2 = f2 + g2;
      d1  = 1.0 / sqrt( f2 * fg2 );
      c   = f2 * d1;
      fg2 = fg2 * d1;
      r   = f * fg2;
      s   = f * d1;
      s   = conj(g) * s;
    }
  }

}
