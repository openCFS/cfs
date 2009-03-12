SET(LAPACK_FOUND 0)

#-------------------------------------------------------------------------------
# Look for lapack lib.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("Netlib BLAS/LAPACK"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a"
  "${CFS_DEPS_ROOT}/lapack/build_lapack.pl"
  "build_lapack.log")

#-------------------------------------------------------------------------------
# Determine paths of LAPACK libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(NETLIB_LAPACK_LIBRARY_DEBUG
    NAMES lapackd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(NETLIB_LAPACK_LIBRARY_RELEASE
    NAMES lapack.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(NETLIB_BLAS_LIBRARY_DEBUG
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libblas.a
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_BLAS_LIBRARY_RELEASE
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libblas.a
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_LAPACK_LIBRARY_DEBUG
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a
    CACHE FILEPATH "Netlib LAPACK library.")
  SET(NETLIB_LAPACK_LIBRARY_RELEASE
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a
    CACHE FILEPATH "Netlib LAPACK library.")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of LAPACK libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(NETLIB_BLAS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(NETLIB_BLAS_LIBRARY_RELEASE)
MARK_AS_ADVANCED(NETLIB_LAPACK_LIBRARY_DEBUG)
MARK_AS_ADVANCED(NETLIB_LAPACK_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set LAPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
IF(DEBUG)
  SET(BLAS_LIBRARY "${NETLIB_BLAS_LIBRARY_DEBUG}")
  SET(LAPACK_LIBRARY "${NETLIB_LAPACK_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(BLAS_LIBRARY "${NETLIB_BLAS_LIBRARY_RELEASE}")
  SET(LAPACK_LIBRARY "${NETLIB_LAPACK_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs
IF(EXISTS "${LAPACK_LIBRARY}")
  SET(LAPACK_FOUND 1)
ELSE(EXISTS "${LAPACK_LIBRARY}")
  MESSAGE(FATAL_ERROR "Could not find ${LAPACK_LIBRARY}")
ENDIF(EXISTS "${LAPACK_LIBRARY}")
ENDIF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
