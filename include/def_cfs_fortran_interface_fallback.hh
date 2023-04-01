#ifndef CFS_FORTRAN_INTERFACE_HEADER_INCLUDED
#define CFS_FORTRAN_INTERFACE_HEADER_INCLUDED

/* Mangling for Fortran global symbols without underscores. */
#define CFS_FORTRAN_INTERFACE_GLOBAL(name,NAME) name##_

/* Mangling for Fortran global symbols with underscores. */
#define CFS_FORTRAN_INTERFACE_GLOBAL_(name,NAME) name##_

/* Mangling for Fortran module symbols without underscores. */
#define CFS_FORTRAN_INTERFACE_MODULE(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/* Mangling for Fortran module symbols with underscores. */
#define CFS_FORTRAN_INTERFACE_MODULE_(mod_name,name, mod_NAME,NAME) __##mod_name##_MOD_##name

/*--------------------------------------------------------------------------*/
/* Mangle some symbols automatically.                                       */
#define dgecon CFS_FORTRAN_INTERFACE_GLOBAL(dgecon, DGECON)
#define dgemm CFS_FORTRAN_INTERFACE_GLOBAL(dgemm, DGEMM)
#define zgemm CFS_FORTRAN_INTERFACE_GLOBAL(zgemm, ZGEMM)
#define dgemv CFS_FORTRAN_INTERFACE_GLOBAL(dgemv, DGEMV)
#define zgemv CFS_FORTRAN_INTERFACE_GLOBAL(zgemv, ZGEMV)
#define dgetrf CFS_FORTRAN_INTERFACE_GLOBAL(dgetrf, DGETRF)
#define dgetri CFS_FORTRAN_INTERFACE_GLOBAL(dgetri, DGETRI)
#define dpotrf CFS_FORTRAN_INTERFACE_GLOBAL(dpotrf, DPOTRF)
#define dpotri CFS_FORTRAN_INTERFACE_GLOBAL(dpotri, DPOTRI)
#define dposv CFS_FORTRAN_INTERFACE_GLOBAL(dposv, DPOSV)
#define zsysv CFS_FORTRAN_INTERFACE_GLOBAL(zsysv, ZSYSV)
#define zhesv CFS_FORTRAN_INTERFACE_GLOBAL(zhesv, ZHESV)
#define zgesv CFS_FORTRAN_INTERFACE_GLOBAL(zgesv, ZGESV)
#define zheev CFS_FORTRAN_INTERFACE_GLOBAL(zheev, ZHEEV)
#define dgbtrf CFS_FORTRAN_INTERFACE_GLOBAL(dgbtrf, DGBTRF)
#define zgbtrf CFS_FORTRAN_INTERFACE_GLOBAL(zgbtrf, ZGBTRF)
#define dgbequ CFS_FORTRAN_INTERFACE_GLOBAL(dgbequ, DGBEQU)
#define zgbequ CFS_FORTRAN_INTERFACE_GLOBAL(zgbequ, ZGBEQU)
#define dlaqgb CFS_FORTRAN_INTERFACE_GLOBAL(dlaqgb, DLAQGB)
#define zlaqgb CFS_FORTRAN_INTERFACE_GLOBAL(zlaqgb, ZLAQGB)
#define dgbtrs CFS_FORTRAN_INTERFACE_GLOBAL(dgbtrs, DGBTRS)
#define zgbtrs CFS_FORTRAN_INTERFACE_GLOBAL(zgbtrs, ZGBTRS)
#define dgbrfs CFS_FORTRAN_INTERFACE_GLOBAL(dgbrfs, DGBRFS)
#define zgbrfs CFS_FORTRAN_INTERFACE_GLOBAL(zgbrfs, ZGBRFS)
#define dlartg CFS_FORTRAN_INTERFACE_GLOBAL(dlartg, DLARTG)
#define zlartg CFS_FORTRAN_INTERFACE_GLOBAL(zlartg, ZLARTG)
#define dpbtrf CFS_FORTRAN_INTERFACE_GLOBAL(dpbtrf, DPBTRF)
#define zpbtrf CFS_FORTRAN_INTERFACE_GLOBAL(zpbtrf, ZPBTRF)
#define dpbtrs CFS_FORTRAN_INTERFACE_GLOBAL(dpbtrs, DPBTRS)
#define zpbtrs CFS_FORTRAN_INTERFACE_GLOBAL(zpbtrs, ZPBTRS)
#define dsyev CFS_FORTRAN_INTERFACE_GLOBAL(dsyev, DSYEV)
#define zgetrf CFS_FORTRAN_INTERFACE_GLOBAL(zgetrf, ZGETRF)
#define zgetri CFS_FORTRAN_INTERFACE_GLOBAL(zgetri, ZGETRI)
#define ilaver CFS_FORTRAN_INTERFACE_GLOBAL(ilaver, ILAVER)
#define cblas_zdscal CFS_FORTRAN_INTERFACE_GLOBAL_(cblas_zdscal, CBLAS_ZDSCAL)
#define cblas_dznrm2 CFS_FORTRAN_INTERFACE_GLOBAL_(cblas_dznrm2, CBLAS_DZNRM2)
#define dsaupd CFS_FORTRAN_INTERFACE_GLOBAL(dsaupd, DSAUPD)
#define dseupd CFS_FORTRAN_INTERFACE_GLOBAL(dseupd, DSEUPD)
#define debug CFS_FORTRAN_INTERFACE_GLOBAL(debug, DEBUG)
#define znaupd CFS_FORTRAN_INTERFACE_GLOBAL(znaupd, ZNAUPD)
#define zneupd CFS_FORTRAN_INTERFACE_GLOBAL(zneupd, ZNEUPD)
#define pardisoinit CFS_FORTRAN_INTERFACE_GLOBAL(pardisoinit, PARDISOINIT)
#define pardiso CFS_FORTRAN_INTERFACE_GLOBAL(pardiso, PARDISO)
#define dscal CFS_FORTRAN_INTERFACE_GLOBAL(dscal, DSCAL)
#define dgesvd CFS_FORTRAN_INTERFACE_GLOBAL(dgesvd, DGESVD)
#define scpip30 CFS_FORTRAN_INTERFACE_GLOBAL(scpip30, SCPIP30)
#define snopta CFS_FORTRAN_INTERFACE_GLOBAL(snopta, SNOPTA)
#define sninit CFS_FORTRAN_INTERFACE_GLOBAL(sninit, SNINIT)
#define sngeti CFS_FORTRAN_INTERFACE_GLOBAL(sngeti, SNGETI)
#define sngetr CFS_FORTRAN_INTERFACE_GLOBAL(sngetr, SNGETR)
#define snset CFS_FORTRAN_INTERFACE_GLOBAL(snset, SNSET)
#define sngetc CFS_FORTRAN_INTERFACE_GLOBAL(sngetc, SNGETC)
#define snseti CFS_FORTRAN_INTERFACE_GLOBAL(snseti, SNSETI)
#define snsetr CFS_FORTRAN_INTERFACE_GLOBAL(snsetr, SNSETR)
#define snspec CFS_FORTRAN_INTERFACE_GLOBAL(snspec, SNSPEC)
#define snmema CFS_FORTRAN_INTERFACE_GLOBAL(snmema, SNMEMA)
#define snjac CFS_FORTRAN_INTERFACE_GLOBAL(snjac, SNJAC)
#define snopenappend CFS_FORTRAN_INTERFACE_GLOBAL(snopenappend, SNOPENAPPEND)
#define snclose CFS_FORTRAN_INTERFACE_GLOBAL(snclose, SNCLOSE)


#endif // CFS_FORTRAN_INTERFACE_HEADER_INCLUDED
