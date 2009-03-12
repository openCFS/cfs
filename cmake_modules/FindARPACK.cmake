SET(ARPACK_FOUND 0)

#-------------------------------------------------------------------------------
# Look for ARPACK header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("ARPACK"
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libarpack.a"
  "${CFS_DEPS_ROOT}/arpack/build_arpack.pl"
  "build_arpack.log")

SET(ARPACK_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of ARPACK libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(ARPACK_LIBRARY_DEBUG
    NAMES arpackd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(ARPACK_LIBRARY_RELEASE
    NAMES arpack.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(ARPACK_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libarpack.a"
    CACHE FILEPATH "ARPACK library")
  SET(ARPACK_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libarpack.a"
    CACHE FILEPATH "ARPACK library")

ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of ARPACK libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(ARPACK_LIBRARY_DEBUG)
MARK_AS_ADVANCED(ARPACK_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set ARPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(ARPACK_LIBRARY "${ARPACK_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(ARPACK_LIBRARY "${ARPACK_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${ARPACK_LIBRARY}")
  
  SET(ARPACK_FOUND 1)
  
ELSE(EXISTS "${ARPACK_LIBRARY}")
  
  MESSAGE(FATAL_ERROR "Could not find ${ARPACK_LIBRARY}")
  
ENDIF(EXISTS "${ARPACK_LIBRARY}")
