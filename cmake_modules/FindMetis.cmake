SET(METIS_FOUND 0)

#-------------------------------------------------------------------------------
# Look for metis header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("METIS"
  "${CFS_BINARY_DIR}/include/metis.h"
  "${CFS_DEPS_ROOT}/metis/build_metis.pl"
  "build_metis.log")

SET(METIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of METIS libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(METIS_LIBRARY_DEBUG
    NAMES metisd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(METIS_LIBRARY_RELEASE
    NAMES metis.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(METIS_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libmetis.a")
  SET(METIS_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libmetis.a")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of METIS libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(METIS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(METIS_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set METIS_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(METIS_LIBRARY "${METIS_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(METIS_LIBRARY "${METIS_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${METIS_LIBRARY}"
    AND EXISTS "${METIS_INCLUDE_DIR}/metis.h")
  
  SET(METIS_FOUND 1)
  
ELSE(EXISTS "${METIS_LIBRARY}"
    AND EXISTS "${METIS_INCLUDE_DIR}/metis.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${METIS_LIBRARY} or ${METIS_INCLUDE_DIR}/metis.h")
  
ENDIF(EXISTS "${METIS_LIBRARY}"
  AND EXISTS "${METIS_INCLUDE_DIR}/metis.h")

#-----------------------------------------------------------------------------
# Determine version of METIS by preprocessing its version header
#-----------------------------------------------------------------------------
SET(METIS_TEST_FILE "${CFS_BINARY_DIR}/include/get_metis_infos.hh")
FILE(WRITE ${METIS_TEST_FILE}
  "#include <${METIS_INCLUDE_DIR}/defs.h>\n\n"
  "<CFS_METIS_VERSION>METISTITLE</CFS_METIS_VERSION>\n")

FILE(TO_NATIVE_PATH ${METIS_TEST_FILE} METIS_TEST_FILE)

IF(CMAKE_VC_COMPILER_TESTS_RUN)
  SET(METIS_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP \"${METIS_TEST_FILE}\"")
ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
  SET(METIS_VERSION_CMD "${CMAKE_CXX_COMPILER} -E \"${METIS_TEST_FILE}\"")
ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)

EXEC_PROGRAM(${METIS_VERSION_CMD}
  ARGS
  OUTPUT_VARIABLE METIS_VERSION_STR
  RETURN_VALUE RETVAL)

#  MESSAGE("METIS_VERSION_STR ${METIS_VERSION_STR}")
FILE(REMOVE ${METIS_TEST_FILE})

STRING(REGEX MATCH
  "<CFS_METIS_VERSION>.*</CFS_METIS_VERSION>"
  CFS_METIS_VERSION "${METIS_VERSION_STR}")
STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
  CFS_METIS_VERSION "${CFS_METIS_VERSION}")

# MESSAGE("CFS_METIS_VERSION ${CFS_METIS_VERSION}")
