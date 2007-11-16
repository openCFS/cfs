SET(ACML_FOUND 0)

#-------------------------------------------------------------------------------
# Determine paths of ACML libraries.
#-------------------------------------------------------------------------------
FIND_LIBRARY(BLAS_ACML_SERIAL_LIB
  NAMES acml
  PATHS ${FORTRAN_POSSIBLE_LIB_PATHS}
  )

FIND_LIBRARY(BLAS_ACML_OPENMP_LIB
  NAMES acml_mp
  PATHS ${FORTRAN_POSSIBLE_LIB_PATHS}
  )

MARK_AS_ADVANCED(BLAS_ACML_SERIAL_LIB)
MARK_AS_ADVANCED(BLAS_ACML_OPENMP_LIB)

IF(USE_OPENMP)
  IF(BLAS_ACML_OPENMP_LIB)
    SET(ACML_FOUND 1)
    SET(LAPACK_ACML_OPENMP_LIB "${BLAS_ACML_OPENMP_LIB}" CACHE FILEPATH "AMD Performance Library (ACML).")
    MARK_AS_ADVANCED(LAPACK_ACML_OPENMP_LIB)
  ENDIF(BLAS_ACML_OPENMP_LIB)
ELSE(USE_OPENMP)
  IF(BLAS_ACML_SERIAL_LIB)
    SET(ACML_FOUND 1)
    SET(LAPACK_ACML_SERIAL_LIB "${BLAS_ACML_SERIAL_LIB}" CACHE FILEPATH "AMD Performance Library (ACML).")
    MARK_AS_ADVANCED(LAPACK_ACML_SERIAL_LIB)
  ENDIF(BLAS_ACML_SERIAL_LIB)
ENDIF(USE_OPENMP)

IF(ACML_FOUND)
  IF(USE_OPENMP)
    SET(ACML_HEADER_FILE "${BLAS_ACML_OPENMP_LIB}")
  ELSE(USE_OPENMP)
    SET(ACML_HEADER_FILE "${BLAS_ACML_SERIAL_LIB}")
  ENDIF(USE_OPENMP)

  STRING(REGEX REPLACE "/gfortran64.*" "/../include/acml.h"
    ACML_HEADER_FILE
    "${ACML_HEADER_FILE}")
#  MESSAGE("ACML_HEADER_FILE ${ACML_HEADER_FILE} SERIAL")

  #-----------------------------------------------------------------------------
  # Determine version of ACML by reading its version header
  #-----------------------------------------------------------------------------
  FILE(READ
    "${ACML_HEADER_FILE}"
    ACML_HEADER
    )

#  MESSAGE("ACML_HEADER ${ACML_HEADER}")

  STRING(REGEX MATCH " ACML version [0-9]+\\.[0-9]+\\.[0-9]+"
    CFS_ACML_VERSION
    "${ACML_HEADER}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
    CFS_ACML_VERSION "${CFS_ACML_VERSION}")

#  MESSAGE("CFS_ACML_VERSION ${CFS_ACML_VERSION}")
  
ENDIF(ACML_FOUND)