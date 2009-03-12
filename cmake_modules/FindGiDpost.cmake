SET(GIDPOST_FOUND 0)

#-------------------------------------------------------------------------------
# Look for metis header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("GiDpost"
  "${CFS_BINARY_DIR}/include/gidpost.h"
  "${CFS_DEPS_ROOT}/gidpost/build_gidpost.pl"
  "build_gidpost.log")

SET(GIDPOST_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of GIDPOST libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(GIDPOST_LIBRARY_DEBUG
    NAMES gidpostd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(GIDPOST_LIBRARY_RELEASE
    NAMES gidpost.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(GIDPOST_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
    CACHE FILEPATH "GiDpost library")
  SET(GIDPOST_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
    CACHE FILEPATH "GiDpost library")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of GIDPOST libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(GIDPOST_LIBRARY_DEBUG)
MARK_AS_ADVANCED(GIDPOST_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set GIDPOST_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${GIDPOST_LIBRARY}"
    AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
  SET(GIDPOST_FOUND 1)
  
ELSE(EXISTS "${GIDPOST_LIBRARY}"
    AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${GIDPOST_LIBRARY} or ${GIDPOST_INCLUDE_DIR}/gidpost.h")
  
ENDIF(EXISTS "${GIDPOST_LIBRARY}"
  AND EXISTS "${GIDPOST_INCLUDE_DIR}/gidpost.h")

IF(GIDPOST_FOUND)
  #-----------------------------------------------------------------------------
  # Determine version of GiDpost by reading its version header
  #-----------------------------------------------------------------------------
  FILE(READ
    "${GIDPOST_INCLUDE_DIR}/gidpost.h"
    GIDPOST_HEADER
    )

  STRING(REGEX MATCH " gidpost [0-9]+\\.[0-9]+"
    CFS_GIDPOST_VERSION
    "${GIDPOST_HEADER}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+"
        CFS_GIDPOST_VERSION "${CFS_GIDPOST_VERSION}")

#  MESSAGE("CFS_GIDPOST_VERSION ${CFS_GIDPOST_VERSION}")

ENDIF(GIDPOST_FOUND)

IF(NOT GIDPOST_FOUND)
  MESSAGE("Warning: GIDPOST could not be found! Please specify proper paths.")
ENDIF(NOT GIDPOST_FOUND)
