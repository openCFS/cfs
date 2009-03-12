SET(GOTOBLAS_FOUND 0)

#-------------------------------------------------------------------------------
# Look for metis header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("GotoBLAS"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgoto.a"
  "${CFS_DEPS_ROOT}/gotoblas/build_gotoblas.pl"
  "build_gotoblas.log")

#-------------------------------------------------------------------------------
# Determine paths of GOTOBLAS libraries.
#-------------------------------------------------------------------------------
SET(GOTOBLAS_SERIAL_LIB
  ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgoto.a;${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a
  CACHE FILEPATH "GotoBLAS serial library.")

SET(GOTOBLAS_OPENMP_LIB
  ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgoto_p.a;${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a
  CACHE FILEPATH "GotoBLAS serial library.")

MARK_AS_ADVANCED(GOTOBLAS_SERIAL_LIB)
MARK_AS_ADVANCED(GOTOBLAS_OPENMP_LIB)

IF(CFS_BLAS_LAPACK STREQUAL "GOTO")
IF(USE_OPENMP)
  IF(GOTOBLAS_OPENMP_LIB)
    SET(GOTOBLAS_FOUND 1)
    SET(BLAS_LIBRARY "${GOTOBLAS_OPENMP_LIB}")
    SET(LAPACK_LIBRARY "${GOTOBLAS_OPENMP_LIB}")
  ENDIF(GOTOBLAS_OPENMP_LIB)
ELSE(USE_OPENMP)
  IF(GOTOBLAS_SERIAL_LIB)
    SET(GOTOBLAS_FOUND 1)
    SET(BLAS_LIBRARY "${GOTOBLAS_SERIAL_LIB}")
    SET(LAPACK_LIBRARY "${GOTOBLAS_SERIAL_LIB}")
  ENDIF(GOTOBLAS_SERIAL_LIB)
ENDIF(USE_OPENMP)
ENDIF(CFS_BLAS_LAPACK STREQUAL "GOTO")
