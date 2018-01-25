#ifndef FILE_BLASLAPACKINTERFACE
#define FILE_BLASLAPACKINTERFACE

#include <def_use_blas.hh>
#include <def_use_lapack.hh>

#include <def_cfs_fortran_interface.hh>

#include <complex>

#include <def_use_adolc.hh>

namespace CoupledField 
{

  //! \file BLASLAPACKInterface.hh
  //! This file contains macros and type definitions required for interfacing
  //! with FORTRAN77 routines. Currently it is required for interfacing with
  //! the BLAS and LAPACK libraries.
  //! Consider using LAPACKE in the future: http://www.netlib.org/lapack/lapacke.html

#ifdef USE_BLAS

#ifdef USE_ADOLC
// methods are currently not defined but could easily be added in the future. "Fleissarbeit"
    inline void dgemm(char*, char*, int*, int*, int*, Double*, Double*, int* ,
               Double*, int*, Double*, Double*, int*){EXCEPTION("Not implemented.");}
    
    // complex matrix matrix multiplication
    inline void zgemm(char*, char*, int*, int*, int*, std::complex<Double>*,
               std::complex<Double>*,int*, std::complex<Double>*,int*,
               std::complex<Double>*, std::complex<Double>*,int*){EXCEPTION("Not implemented.");}
    
    // matrix vector multiplication
    inline void dgemv(char*, int*, int*, Double*, Double*, int*, const Double*,
               int*, Double*, Double*, int*){EXCEPTION("Not implemented.");}
    
    // complex matrix vector multiplication
    inline void zgemv(char*, int*, int*, std::complex<Double>*, std::complex<Double>*,
               int*, const std::complex<Double>*, int*, std::complex<Double>*,
               std::complex<Double>*,int*){EXCEPTION("Not implemented.");}
    
#else
    
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
    
  }
#endif
#endif

#ifdef USE_LAPACK

#ifdef USE_ADOLC


    //! LU decomoposition of a general matrix
    void dgetrf( int* M, int *N, Double* A, int* lda, int* IPIV, int* INFO ){EXCEPTION("Not implemented.");}

    void zgetrf( const int* m, const int* n, std::complex<Double>* a,
                 const int* lda, int* ipiv, int* info ){EXCEPTION("Not implemented.");}

    //! generate inverse of a matrix given its LU decomposition
    void dgetri( int* N, Double* A, int* lda, int* IPIV, Double* WORK, int* lwork, int* INFO ){EXCEPTION("Not implemented.");}

    void zgetri( const int* n, std::complex<Double>* a, const int* lda,
                 const int* ipiv, std::complex<Double>* work, const int* lwork,
				 int* info ){EXCEPTION("Not implemented.");}

    //! Cholesky decomposition of a symmetric matrix
    void dpotrf( char*, int*, Double*, int*, int* ){EXCEPTION("Not implemented.");}
    
    //! generate inverse of a matrix given its Cholesky decomposition
    void dpotri( char*, int*, Double* , int*, int* ){EXCEPTION("Not implemented.");}
    
    //! solves positive definite system
    void dposv( char*, int*, int*, Double*, int*, Double*, int*, int* ){EXCEPTION("Not implemented.");}
     
    //! solves symmetric problems with complex coefficients
    void zsysv( char*, int*, int*, std::complex<Double>*, int*, int*,
                std::complex<Double>*, int*, std::complex<Double>*, int*, int* ){EXCEPTION("Not implemented.");}
    
    //! solves hermitian problems
    void zhesv( char*, int*, int*, std::complex<Double>*, int*, int*,
                std::complex<Double>*, int*, std::complex<Double>*, int*, int* ){EXCEPTION("Not implemented.");}
    
    //! solves general kind of problems with complex coefficients
    void zgesv( int*, int*, std::complex<Double>*, int*, int*, std::complex<Double>*, int*, int*){EXCEPTION("Not implemented.");}
    //    Double LP_DLAMCH (char *CMACHp); // DLAMCH - determine Double precision machine parameters
    
    //! determines eigenvalues of a complex - valued matrix
    void zheev( char*, char*, int*, std::complex<Double>*, int *, Double*,
                std::complex<Double>*, int*, Double* ,int* ){EXCEPTION("Not implemented.");}

    void dgbtrf( int*, int*, int*, int*, Double    *, int*, int*, int* );
    void zgbtrf( int*, int*, int*, int*, std::complex<Double> *, int*, int*, int* ){EXCEPTION("Not implemented.");}

    void dgbequ( int*, int*, int*, int*, Double*, int*, Double*,
                 Double*, Double*, Double*, Double*, int* ){EXCEPTION("Not implemented.");}
    void zgbequ( int*, int*, int*, int*, std::complex<Double>*, int*, Double*,
                 Double*, Double*, Double*, Double*, int* ){EXCEPTION("Not implemented.");}

    void dlaqgb( int*, int*, int*, int*, Double*, int*, Double*,
                 Double*, Double*, Double*, Double*, char* ){EXCEPTION("Not implemented.");}
    void zlaqgb( int*, int*, int*, int*, std::complex<Double>*, int*, Double*,
                 Double*, Double*, Double*, Double*, char* ){EXCEPTION("Not implemented.");}

    void dgbtrs( char*, int*, int*, int*, int*, Double*, int*, int*,
                 Double*, int*, int* ){EXCEPTION("Not implemented.");}
    void zgbtrs( char*, int*, int*, int*, int*, std::complex<Double>*, int*, int*,
                 std::complex<Double>*, int*, int* ){EXCEPTION("Not implemented.");}

    void dgbrfs( char*, int*, int*, int*, int*, const Double*, int*,
                 Double*, int*, int*, const Double*, int*, Double*,
                 int*, Double*, Double*, Double*, int*, int* ){EXCEPTION("Not implemented.");}
    void zgbrfs( char*, int*, int*, int*, int*, const std::complex<Double>*, int*,
                 std::complex<Double>*, int*, int*, const std::complex<Double>*, int*,
                 std::complex<Double>*, int*, Double*, Double*, std::complex<Double>*,
                 Double*, int* ){EXCEPTION("Not implemented.");}

    void dlartg( Double*, Double*, Double*, Double*, Double* ){EXCEPTION("Not implemented.");}
    void zlartg( std::complex<Double>*, std::complex<Double>*, Double*, std::complex<Double>*,
                 std::complex<Double>* ){EXCEPTION("Not implemented.");}

    void dpbtrf( char*, int*, int*, Double*, int*, int* ){EXCEPTION("Not implemented.");}
    void zpbtrf( char*, int*, int*, Double*, int*, int* ){EXCEPTION("Not implemented.");}

    void dpbtrs( char*, int*, int*, int*, Double*, int*, Double*, int*,
                 int*){EXCEPTION("Not implemented.");}
    void zpbtrs( char*, int*, int*, int*, Double*, int*, Double*, int*,
                 int*){EXCEPTION("Not implemented.");}
    void dsyev( char*, char*, int*, Double*, int*, Double*, Double*, int*, int* ){EXCEPTION("Not implemented.");}

    void dgecon( char*, int*, Double*, int*, Double*, Double*, Double*, int*, int*){EXCEPTION("Not implemented.");}

#else


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

    // ! To be continued ...

  }
#endif
#endif
}
#endif
