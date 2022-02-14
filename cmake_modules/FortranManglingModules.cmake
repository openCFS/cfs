# this file is included in CheckFortranRuntime.cmake
# whenever you use a new call from a fortran lib, add it here 
# and consider updating include/def_cfs_fortran_interface_fallback.hh

FortranCInterface_HEADER("${CFS_BINARY_DIR}/include/def_cfs_fortran_interface.hh"
  MACRO_NAMESPACE "CFS_FORTRAN_INTERFACE_"
  SYMBOLS
  # BLAS and LAPACK
  dgecon
  dgemm
  zgemm
  dgemv
  zgemv
  dgetrf
  dgetri
  dpotrf
  dpotri
  dposv
  zsysv
  zhesv
  zgesv
  zheev
  dgbtrf
  zgbtrf
  dgbequ
  zgbequ
  dlaqgb
  zlaqgb
  dgbtrs
  zgbtrs
  dgbrfs
  zgbrfs
  dlartg
  zlartg
  dpbtrf
  zpbtrf
  dpbtrs
  zpbtrs
  dsyev
  zgetrf
  zgetri
  ilaver
  cblas_zdscal
  cblas_dznrm2
  # ARPACK
  dsaupd
  dseupd
  debug
  znaupd
  zneupd
  # Pardiso
  pardisoinit
  pardiso
  # additionally for PALM
  dscal
  dgesvd
  # scpip
  scpip30
  # snopt
  snopta
  sninit
  sngeti
  sngetr
  snset
  sngetc
  snseti
  snsetr
  snspec
  snmema
  snjac
  snopenappend
  snclose
)

