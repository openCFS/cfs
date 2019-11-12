#ifndef FILE_BLASLAPACKINTERFACE
#define FILE_BLASLAPACKINTERFACE

#include <def_use_blas.hh>

#include <def_cfs_fortran_interface.hh>

#include <complex>

namespace CoupledField 
{

  //! \file BLASLAPACKInterface.hh
  //! This file contains macros and type definitions required for interfacing
  //! with FORTRAN77 routines. Currently it is required for interfacing with
  //! the BLAS and LAPACK libraries.
  //! Consider using LAPACKE in the future: http://www.netlib.org/lapack/lapacke.html

   //! Make sure to add the function also to cmake_modules/CheckFortranRuntime.cmake
   //! to have proper name mangling also for netlib and openblas
  extern "C"{
    // matrix matrix multiplication
    void dgemm(char*, char*, int*, int*, int*, double*, double*, int* ,
               double*, int*, double*, double*, int*);
    
    // complex matrix matrix multiplication
    void zgemm(char*, char*, int*, int*, int*, std::complex<double>*,
               std::complex<double>*,int*, std::complex<double>*,int*,
               std::complex<double>*, std::complex<double>*,int*);
    
    // matrix vector multiplication
    void dgemv(char*, int*, int*, double*, double*, int*, const double*,
               int*, double*, double*, int*);
    
    // complex matrix vector multiplication
    void zgemv(char*, int*, int*, std::complex<double>*, std::complex<double>*,
               int*, const std::complex<double>*, int*, std::complex<double>*,
               std::complex<double>*,int*);
    
    void cblas_zdscal(const int __N, const double __alpha, void *__X, const int __incX);
    double cblas_dznrm2(const int __N, const void *__X, const int __incX);

  }


  //! Prototypes for LAPACK routines
  extern "C" {
    
    //! LU decomoposition of a general matrix
    void dgetrf( int* M, int *N, double* A, int* lda, int* IPIV, int* INFO );

    void zgetrf( const int* m, const int* n, std::complex<double>* a,
                 const int* lda, int* ipiv, int* info );

    //! generate inverse of a matrix given its LU decomposition
    void dgetri( int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO );

    void zgetri( const int* n, std::complex<double>* a, const int* lda,
                 const int* ipiv, std::complex<double>* work, const int* lwork,
				 int* info );

    //! Cholesky decomposition of a symmetric matrix
    void dpotrf( char*, int*, double*, int*, int* );
    
    //! generate inverse of a matrix given its Cholesky decomposition
    void dpotri( char*, int*, double* , int*, int* );
    
    //! solves positive definite system
    void dposv( char*, int*, int*, double*, int*, double*, int*, int* );
     
    //! solves symmetric problems with complex coefficients
    void zsysv( char*, int*, int*, std::complex<double>*, int*, int*,
                std::complex<double>*, int*, std::complex<double>*, int*, int* );
    
    //! solves hermitian problems
    void zhesv( char*, int*, int*, std::complex<double>*, int*, int*,
                std::complex<double>*, int*, std::complex<double>*, int*, int* );
    
    //! solves general kind of problems with complex coefficients
    void zgesv( int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
    //    double LP_DLAMCH (char *CMACHp); // DLAMCH - determine double precision machine parameters
    
    //! determines eigenvalues of a complex - valued matrix
    void zheev( char*, char*, int*, std::complex<double>*, int *, double*,
                std::complex<double>*, int*, double* ,int* ); 

    void dgbtrf( int*, int*, int*, int*, double    *, int*, int*, int* );
    void zgbtrf( int*, int*, int*, int*, std::complex<double> *, int*, int*, int* );

    void dgbequ( int*, int*, int*, int*, double*, int*, double*,
                 double*, double*, double*, double*, int* );
    void zgbequ( int*, int*, int*, int*, std::complex<double>*, int*, double*,
                 double*, double*, double*, double*, int* );

    void dlaqgb( int*, int*, int*, int*, double*, int*, double*,
                 double*, double*, double*, double*, char* );
    void zlaqgb( int*, int*, int*, int*, std::complex<double>*, int*, double*,
                 double*, double*, double*, double*, char* );

    void dgbtrs( char*, int*, int*, int*, int*, double*, int*, int*,
                 double*, int*, int* );
    void zgbtrs( char*, int*, int*, int*, int*, std::complex<double>*, int*, int*,
                 std::complex<double>*, int*, int* );

    void dgbrfs( char*, int*, int*, int*, int*, const double*, int*,
                 double*, int*, int*, const double*, int*, double*,
                 int*, double*, double*, double*, int*, int* );
    void zgbrfs( char*, int*, int*, int*, int*, const std::complex<double>*, int*,
                 std::complex<double>*, int*, int*, const std::complex<double>*, int*,
                 std::complex<double>*, int*, double*, double*, std::complex<double>*,
                 double*, int* );

    void dlartg( double*, double*, double*, double*, double* );
    void zlartg( std::complex<double>*, std::complex<double>*, double*, std::complex<double>*,
                 std::complex<double>* );

    void dpbtrf( char*, int*, int*, double*, int*, int* );
    void zpbtrf( char*, int*, int*, double*, int*, int* );

    void dpbtrs( char*, int*, int*, int*, double*, int*, double*, int*,
                 int*);
    void zpbtrs( char*, int*, int*, int*, double*, int*, double*, int*,
                 int*);
    void dsyev( char*, char*, int*, double*, int*, double*, double*, int*, int* );

    void dgecon( char*, int*, double*, int*, double*, double*, double*, int*, int*);


    void dgesvd(char* , char* , int *, int *, double*, int *, double*,
                      double*, int *, double*, int *, double*, int *, int *);

    int dscal(int *, double *, double *, int *);

    double dznrm2_( int*, std::complex<double>*, int*);

    void zdscal_(int* N,double* DA,std::complex<double>* ZX, int* INCX);


    void znaupd(int *IDO, const char *BMAT,int *N,const char *WHICH,int *NEV,double *TOL,std::complex<double> *RESID,int *NCV,std::complex<double> *V,
                 int *LDV,int *IPARAM,int *IPNTR,std::complex<double> *WORKD,std::complex<double> *WORKL,int *LWORKL, double *RWORK,
                 int *INFO,int _BMAT, int _WHICH    // The length of the actual WHICH argument
                );

    void zneupd(bool *RVEC,
                 const char *HOWMNY,
                 bool *SELECT,
                 std::complex<double> *D,
				 std::complex<double> *Z,
                 int *LDZ,
				 std::complex<double> *SIGMA,
				 std::complex<double> *WORKEV,
                 const char *BMAT,
                 int *N,
                 const char *WHICH,
                 int *NEV,
                 double *TOL,
				 std::complex<double> *RESID,
                 int *NCV,
				 std::complex<double> *V,
                 int *LDV,
                 int *IPARAM,
                 int *IPNTR,
                 std::complex<double> *WORKD,
				 std::complex<double> *WORKL,
                 int *LWORKL,
                 double *RWORK,
                 int *INFO
                );
    // ! To be continued ...
  }

}
#endif
