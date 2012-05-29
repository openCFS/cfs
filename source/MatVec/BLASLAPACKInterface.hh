#ifndef FILE_BLASLAPACKINTERFACE
#define FILE_BLASLAPACKINTERFACE

#include <def_use_blas.hh>
#include <def_use_lapack.hh>

#include <def_fortran_interface.hh>

#include <complex>

namespace CoupledField 
{

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

#ifdef USE_BLAS

  extern "C"{
    // matrix matrix multiplication
    void FORTRAN_CALL F77NAME(dgemm)(char*, char*, int*, int*, int*, double*, double*, int* , double*, int*, double*, double*, int*);
    
    // complex matrix matrix multiplication
    void FORTRAN_CALL F77NAME(zgemm)(char*, char*, int*, int*, int*, std::complex<double>*, std::complex<double>*,int*, std::complex<double>*,int*,std::complex<double>*, std::complex<double>*,int*);
    
    // matrix vector multiplication
    void FORTRAN_CALL F77NAME(dgemv)(char*, int*, int*, double*, double*, int*, const double*, int*, double*, double*, int*);
    
    // complex matrix vector multiplication
    void FORTRAN_CALL F77NAME(zgemv)(char*, int*, int*, std::complex<double>*, std::complex<double>*,int*, const std::complex<double>*, int*, std::complex<double>*, std::complex<double>*,int*);
    
  }

#endif

#ifdef USE_LAPACK

  //! Prototypes for LAPACK routines
  extern "C" {
    
    //! LU decomoposition of a general matrix
    void FORTRAN_CALL F77NAME(dgetrf)( int* M, int *N, double* A, int* lda, int* IPIV, int* INFO );

    //! generate inverse of a matrix given its LU decomposition
    void FORTRAN_CALL F77NAME(dgetri)( int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO );

    //! Cholesky decomposition of a symmetric matrix
    void FORTRAN_CALL F77NAME(dpotrf)( char*, int*, double*, int*, int* );
    
    //! generate inverse of a matrix given its Cholesky decomposition
    void FORTRAN_CALL F77NAME(dpotri)( char*, int*, double* , int*, int* );
    
    //! solves positive definite system
    void FORTRAN_CALL F77NAME(dposv)( char*, int*, int*, double*, int*, double*, int*, int* );
     
    //! solves symmetric problems with complex coefficients
    void FORTRAN_CALL F77NAME(zsysv)( char*, int*, int*, F77complex16*, int*, int*, F77complex16*, int*, F77complex16*, int*, int* );
    
    //! solves hermitian problems
    void FORTRAN_CALL F77NAME(zhesv)( char*, int*, int*, F77complex16*, int*, int*, F77complex16*, int*, F77complex16*, int*, int* );
    
    //! solves general kind of problems with complex coefficients
    void FORTRAN_CALL F77NAME(zgesv)( int*, int*, F77complex16*, int*, int*, F77complex16*, int*, int*);
    //    double LP_DLAMCH (char *CMACHp); // DLAMCH - determine double precision machine parameters
    
    //! determines eigenvalues of a complex - valued matrix
    void FORTRAN_CALL F77NAME(zheev)( char*, char*, int*, F77complex16*, int *, F77real8*, F77complex16*, int*, F77real8* ,int* ); 

    void FORTRAN_CALL F77NAME(dgbtrf)( int*, int*, int*, int*, F77real8    *, int*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbtrf)( int*, int*, int*, int*, F77complex16 *, int*, int*, int* );

    void FORTRAN_CALL F77NAME(dgbequ)( int*, int*, int*, int*, F77real8*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, int* );
    void FORTRAN_CALL F77NAME(zgbequ)( int*, int*, int*, int*, F77complex16*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, int* );

    void FORTRAN_CALL F77NAME(dlaqgb)( int*, int*, int*, int*, F77real8*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, char* );
    void FORTRAN_CALL F77NAME(zlaqgb)( int*, int*, int*, int*, F77complex16*, int*, F77real8*,
		    F77real8*, F77real8*, F77real8*, F77real8*, char* );

    void FORTRAN_CALL F77NAME(dgbtrs)( char*, int*, int*, int*, int*, F77real8*, int*, int*,
		    F77real8*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbtrs)( char*, int*, int*, int*, int*, F77complex16*, int*, int*,
		    F77complex16*, int*, int* );

    void FORTRAN_CALL F77NAME(dgbrfs)( char*, int*, int*, int*, int*, const F77real8*, int*,
		    F77real8*, int*, int*, const F77real8*, int*, F77real8*,
		    int*, F77real8*, F77real8*, F77real8*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbrfs)( char*, int*, int*, int*, int*, const F77complex16*, int*,
		    F77complex16*, int*, int*, const F77complex16*, int*,
		    F77complex16*, int*, F77real8*, F77real8*, F77complex16*,
		    F77real8*, int* );

    void FORTRAN_CALL F77NAME(dlartg)( F77real8*, F77real8*, F77real8*, F77real8*, F77real8* );
    void FORTRAN_CALL F77NAME(zlartg)( F77complex16*, F77complex16*, F77real8*, F77complex16*,
		    F77complex16* );

    void FORTRAN_CALL F77NAME(dpbtrf)( char*, int*, int*, double*, int*, int* );
    void FORTRAN_CALL F77NAME(zpbtrf)( char*, int*, int*, double*, int*, int* );

    void FORTRAN_CALL F77NAME(dpbtrs)( char*, int*, int*, int*, double*, int*, double*, int*,
                    int*);
    void FORTRAN_CALL F77NAME(zpbtrs)( char*, int*, int*, int*, double*, int*, double*, int*,
                    int*);

    // ! To be continued ...

  }
#endif

}
#endif
