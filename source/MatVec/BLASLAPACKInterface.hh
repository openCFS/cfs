#ifndef FILE_BLASLAPACKINTERFACE
#define FILE_BLASLAPACKINTERFACE

#include <def_use_blas.hh>
#include <def_use_lapack.hh>

#include <def_fortran_interface.hh>

#include <complex>

namespace CoupledField 
{

  //! \file BLASLAPACKInterface.hh
  //! This file contains macros and type definitions required for interfacing
  //! with FORTRAN77 routines. Currently it is required for interfacing with
  //! the BLAS and LAPACK libraries.

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
    void FORTRAN_CALL F77NAME(zsysv)( char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int* );
    
    //! solves hermitian problems
    void FORTRAN_CALL F77NAME(zhesv)( char*, int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, std::complex<double>*, int*, int* );
    
    //! solves general kind of problems with complex coefficients
    void FORTRAN_CALL F77NAME(zgesv)( int*, int*, std::complex<double>*, int*, int*, std::complex<double>*, int*, int*);
    //    double LP_DLAMCH (char *CMACHp); // DLAMCH - determine double precision machine parameters
    
    //! determines eigenvalues of a complex - valued matrix
    void FORTRAN_CALL F77NAME(zheev)( char*, char*, int*, std::complex<double>*, int *, double*, std::complex<double>*, int*, double* ,int* ); 

    void FORTRAN_CALL F77NAME(dgbtrf)( int*, int*, int*, int*, double    *, int*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbtrf)( int*, int*, int*, int*, std::complex<double> *, int*, int*, int* );

    void FORTRAN_CALL F77NAME(dgbequ)( int*, int*, int*, int*, double*, int*, double*,
		    double*, double*, double*, double*, int* );
    void FORTRAN_CALL F77NAME(zgbequ)( int*, int*, int*, int*, std::complex<double>*, int*, double*,
		    double*, double*, double*, double*, int* );

    void FORTRAN_CALL F77NAME(dlaqgb)( int*, int*, int*, int*, double*, int*, double*,
		    double*, double*, double*, double*, char* );
    void FORTRAN_CALL F77NAME(zlaqgb)( int*, int*, int*, int*, std::complex<double>*, int*, double*,
		    double*, double*, double*, double*, char* );

    void FORTRAN_CALL F77NAME(dgbtrs)( char*, int*, int*, int*, int*, double*, int*, int*,
		    double*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbtrs)( char*, int*, int*, int*, int*, std::complex<double>*, int*, int*,
		    std::complex<double>*, int*, int* );

    void FORTRAN_CALL F77NAME(dgbrfs)( char*, int*, int*, int*, int*, const double*, int*,
		    double*, int*, int*, const double*, int*, double*,
		    int*, double*, double*, double*, int*, int* );
    void FORTRAN_CALL F77NAME(zgbrfs)( char*, int*, int*, int*, int*, const std::complex<double>*, int*,
		    std::complex<double>*, int*, int*, const std::complex<double>*, int*,
		    std::complex<double>*, int*, double*, double*, std::complex<double>*,
		    double*, int* );

    void FORTRAN_CALL F77NAME(dlartg)( double*, double*, double*, double*, double* );
    void FORTRAN_CALL F77NAME(zlartg)( std::complex<double>*, std::complex<double>*, double*, std::complex<double>*,
		    std::complex<double>* );

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
