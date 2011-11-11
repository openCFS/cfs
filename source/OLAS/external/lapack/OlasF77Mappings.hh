#ifndef OLASF77MAPPING_HH
#define OLASF77MAPPING_HH

#include <complex>

namespace CoupledField {

  //! \file olasf77mapping.hh
  //! This file contains macros and type definitions required for interfacing
  //! with FORTRAN77 routines. Currently it is only required for interfacing
  //! with the LAPACK library.

  //! Type for real values of double precision. This corresponds to REAL*8.

  //! Type for real values of double precision.
  //! This corresponds to REAL*8.
  typedef double F77real8;

  //! Type for complex values of double precision.

  //! Type for complex values of double precision. This corresponds to
  //! COMPLEX*16. Note that it is not clear whether all FORTRAN77 compilers
  //! represent complex in this fashion. It is the case, however, for the GNU
  //! and Intel compilers. Note also that COMPLEX*16 is no standard type in
  //! FORTRAN77.
  class F77complex16 {
  public:
    F77real8 real;
    F77real8 imag;
    F77complex16() {
      real = 0;
      imag = 0;
    }
    F77complex16 &operator= ( F77complex16 v ) {
      this->real = v.real;
      this->imag = v.imag;
      return *this;
    }
    F77complex16 &operator= ( double v ) {
      this->real = v;
      this->imag = 0;
      return *this;
    }
    F77complex16 &operator+= ( F77complex16 v ) {
      this->real += v.real;
      this->imag += v.imag;
      return *this;
    }
  };

  //@{
  //! Conversion from C++ to our FORTRAN77 data type
  void CC2F77( const double               &v, F77real8     &val );
  void CC2F77( const double               &v, F77complex16 &val );
  void CC2F77( const std::complex<double> &v, F77real8     &val );
  void CC2F77( const std::complex<double> &v, F77complex16 &val );
  //@}

  //@{
  //! Conversion from our FORTRAN77 to C++ data type
  void F772CC( const F77real8     &v, double               &val );
  void F772CC( const F77real8     &v, std::complex<double> &val );
  void F772CC( const F77complex16 &v, std::complex<double> &val );
  void F772CC( const F77complex16 &v, double               &val );
  //@}

  // Define the names of the LAPACK routines
  // The naming scheme unfortunately depends on the compiler!

#define LP_DGBTRF dgbtrf_
#define LP_ZGBTRF zgbtrf_

#define LP_DGBEQU dgbequ_
#define LP_ZGBEQU zgbequ_

#define LP_DLAQGB dlaqgb_
#define LP_ZLAQGB zlaqgb_

#define LP_DGBTRS dgbtrs_
#define LP_ZGBTRS zgbtrs_

#define LP_DGBRFS dgbrfs_
#define LP_ZGBRFS zgbrfs_

#define LP_DLARTG dlartg_
#define LP_ZLARTG zlartg_

#define LP_DPBTRF dpbtrf_
#define LP_ZPBTRF zpbtrf_

#define LP_DPBTRS dpbtrs_
#define LP_ZPBTRS zpbtrs_

  // Generate prototypes for LAPACK routines
  extern "C" {
    void LP_DGBTRF( int*, int*, int*, int*, F77real8    *, int*, int*, int* );
    void LP_ZGBTRF( int*, int*, int*, int*, F77complex16 *, int*, int*, int* );

    void LP_DGBEQU( int*, int*, int*, int*, F77real8*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, int* );
    void LP_ZGBEQU( int*, int*, int*, int*, F77complex16*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, int* );

    void LP_DLAQGB( int*, int*, int*, int*, F77real8*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, char* );
    void LP_ZLAQGB( int*, int*, int*, int*, F77complex16*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, char* );

    void LP_DGBTRS( char*, int*, int*, int*, int*, F77real8*, int*, int*,
		    F77real8*, int*, int* );
    void LP_ZGBTRS( char*, int*, int*, int*, int*, F77complex16*, int*, int*,
		    F77complex16*, int*, int* );

    void LP_DGBRFS( char*, int*, int*, int*, int*, const F77real8*, int*,
		    F77real8*, int*, int*, const F77real8*, int*, F77real8*,
		    int*, F77real8*, F77real8*, F77real8*, int*, int* );
    void LP_ZGBRFS( char*, int*, int*, int*, int*, const F77complex16*, int*,
		    F77complex16*, int*, int*, const F77complex16*, int*,
		    F77complex16*, int*, F77real8*, F77real8*, F77complex16*,
		    F77real8*, int* );

    void LP_DLARTG( F77real8*, F77real8*, F77real8*, F77real8*, F77real8* );
    void LP_ZLARTG( F77complex16*, F77complex16*, F77real8*, F77complex16*,
		    F77complex16* );

    void LP_DPBTRF( char*, int*, int*, double*, int*, int* );
    void LP_ZPBTRF( char*, int*, int*, double*, int*, int* );

    void LP_DPBTRS( char*, int*, int*, int*, double*, int*, double*, int*,
                    int*);
    void LP_ZPBTRS( char*, int*, int*, int*, double*, int*, double*, int*,
                    int*);

  }

}

#endif
