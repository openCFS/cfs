#ifndef BLAS_LAPACK_HEADER_INCLUDED
#define BLAS_LAPACK_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define BLAS_LAPACK_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define BLAS_LAPACK_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define BLAS_LAPACK_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define BLAS_LAPACK_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define dtrsv BLAS_LAPACK_GLOBAL(dtrsv, DTRSV)
#define dgemv BLAS_LAPACK_GLOBAL(dgemv, DGEMV)
#define dtrsm BLAS_LAPACK_GLOBAL(dtrsm, DTRSM)
#define dgemm BLAS_LAPACK_GLOBAL(dgemm, DGEMM)
#define dsyrk BLAS_LAPACK_GLOBAL(dsyrk, DSYRK)
#define dger BLAS_LAPACK_GLOBAL(dger, DGER)
#define dscal BLAS_LAPACK_GLOBAL(dscal, DSCAL)
#define dpotrf BLAS_LAPACK_GLOBAL(dpotrf, DPOTRF)
#define ztrsv BLAS_LAPACK_GLOBAL(ztrsv, ZTRSV)
#define zgemv BLAS_LAPACK_GLOBAL(zgemv, ZGEMV)
#define ztrsm BLAS_LAPACK_GLOBAL(ztrsm, ZTRSM)
#define zgemm BLAS_LAPACK_GLOBAL(zgemm, ZGEMM)
#define zherk BLAS_LAPACK_GLOBAL(zherk, ZHERK)
#define zgeru BLAS_LAPACK_GLOBAL(zgeru, ZGERU)
#define zscal BLAS_LAPACK_GLOBAL(zscal, ZSCAL)
#define zpotrf BLAS_LAPACK_GLOBAL(zpotrf, ZPOTRF)

#endif
