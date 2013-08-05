SET(CCMIO_FOUND 0)

#-------------------------------------------------------------------------------
# Look for metis header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("CCMIO"
  "${CFS_BINARY_DIR}/include/ccmio.h"
  "${CFS_DEPS_ROOT}/ccmio/build_ccmio.pl"
  "build_ccmio.log")

SET(CCMIO_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of CCMIO libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(CCMIO_LIBRARY_DEBUG
    NAMES cgnsd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(CCMIO_LIBRARY_RELEASE
    NAMES cgns.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(CCMIO_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libccmio.a;${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libadf.a"
    CACHE FILEPATH "STARCCM+ I/O Library")
  SET(CCMIO_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libccmio.a;${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libadf.a"
    CACHE FILEPATH "STARCCM+ I/O System Library")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of CCMIO libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(CCMIO_LIBRARY_DEBUG)
MARK_AS_ADVANCED(CCMIO_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set CCMIO_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(CCMIO_LIBRARY "${CCMIO_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(CCMIO_LIBRARY "${CCMIO_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${CCMIO_INCLUDE_DIR}/ccmio.h")
 
  SET(CCMIO_FOUND 1)
  
ELSE(EXISTS "${CCMIO_INCLUDE_DIR}/ccmio.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${CCMIO_LIBRARY} or ${CCMIO_INCLUDE_DIR}/ccmio.h")
  
ENDIF(EXISTS "${CCMIO_INCLUDE_DIR}/ccmio.h")

