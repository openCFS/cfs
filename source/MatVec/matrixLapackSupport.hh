// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MATRIX_LAPACK_SUPPORT
#define FILE_MATRIX_LAPACK_SUPPORT

#include <def_use_lapack.hh>
#include "OLAS/external/lapack/olasf77mapping.hh"

namespace CoupledField 
{

#ifdef USE_LAPACK
 
#define LP_ZSYSV zsysv_ // solves symmetric matrices
#define LP_ZGESV zgesv_ // solves general type of matrices
#define LP_ZHESV zhesv_ // solves  ZHESV computes the solution to a complex system of linear equations
  // A * X = B, where A is an N-by-N Hermitian matrix and X and B are N-by-NRHS matrices
#define LP_DLAMCH dlamch_ // tests data types ...
#define LP_ZHEEV zheev_ //ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
  //  complex Hermitian matrix A.

  //! Prototypes for LAPACK routines
  
  extern "C" {
    //! solves symmetric problems with complex coefficients
    void LP_ZSYSV( char*, int*, int*, F77complex16*, int*, int*, F77complex16*, int*, F77complex16*, int*, int*);
    //! solves hermitian problems
    void LP_ZHESV( char*, int*, int*, F77complex16*, int*, int*, F77complex16*, int*, F77complex16*, int*, int*);
    //! solves general kind of problems with complex coefficients
    void LP_ZGESV( int*, int*, F77complex16*, int*, int*, F77complex16*, int*, int*);
    //    double LP_DLAMCH (char *CMACHp); // DLAMCH - determine double precision machine parameters
    //! determines eigenvalues of a complex - valued matrix
    void LP_ZHEEV( char*, char*, int*, F77complex16*, int *, F77real8*, F77complex16*, int*, F77real8* ,int* ); 
    // ! To be continued ...

  }
}
#endif
#endif
