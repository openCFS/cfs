SET(CGNS_FOUND 0)

#-------------------------------------------------------------------------------
# Look for metis header.
#-------------------------------------------------------------------------------
BUILD_EXTLIB("CGNS & adfviewer"
  "${CFS_BINARY_DIR}/include/cgnslib.h"
  "${CFS_DEPS_ROOT}/cgns/build_cgns.pl"
  "build_cgns.log")

SET(CGNS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of CGNS libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(CGNS_LIBRARY_DEBUG
    NAMES cgnsd.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(CGNS_LIBRARY_RELEASE
    NAMES cgns.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  SET(CGNS_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libcgns.a"
    CACHE FILEPATH "CFD General Notation System Library")
  SET(CGNS_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libcgns.a"
    CACHE FILEPATH "CFD General Notation System Library")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of CGNS libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(CGNS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(CGNS_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set CGNS_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(CGNS_LIBRARY "${CGNS_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(CGNS_LIBRARY "${CGNS_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${CGNS_LIBRARY}"
    AND EXISTS "${CGNS_INCLUDE_DIR}/cgnslib.h")
 
  SET(CGNS_FOUND 1)
  
ELSE(EXISTS "${CGNS_LIBRARY}"
    AND EXISTS "${CGNS_INCLUDE_DIR}/cgnslib.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${CGNS_LIBRARY} or ${CGNS_INCLUDE_DIR}/cgnslib.h")
  
ENDIF(EXISTS "${CGNS_LIBRARY}"
  AND EXISTS "${CGNS_INCLUDE_DIR}/cgnslib.h")

#-----------------------------------------------------------------------------
# Determine version of CGNS by preprocessing its version header
#-----------------------------------------------------------------------------
#SET(CGNS_TEST_FILE "${CFS_BINARY_DIR}/include/get_metis_infos.hh")
#FILE(WRITE ${CGNS_TEST_FILE}
#  "#include <${CGNS_INCLUDE_DIR}/defs.h>\n\n"
#  "<CFS_CGNS_VERSION>CGNSTITLE</CFS_CGNS_VERSION>\n")

#FILE(TO_NATIVE_PATH ${CGNS_TEST_FILE} CGNS_TEST_FILE

#IF(CMAKE_VC_COMPILER_TESTS_RUN)
#  SET(CGNS_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP \"${CGNS_TEST_FILE}\"")
#ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
#  SET(CGNS_VERSION_CMD "${CMAKE_CXX_COMPILER} -E \"${CGNS_TEST_FILE}\"")
#ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)
#
#EXEC_PROGRAM(${CGNS_VERSION_CMD}
#  ARGS
#  OUTPUT_VARIABLE CGNS_VERSION_STR
#  RETURN_VALUE RETVAL)

#  MESSAGE("CGNS_VERSION_STR ${CGNS_VERSION_STR}")
#FILE(REMOVE ${CGNS_TEST_FILE})

#STRING(REGEX MATCH
#  "<CFS_CGNS_VERSION>.*</CFS_CGNS_VERSION>"
#  CFS_CGNS_VERSION "${CGNS_VERSION_STR}")
#STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
#  CFS_CGNS_VERSION "${CFS_CGNS_VERSION}")

# MESSAGE("CFS_CGNS_VERSION ${CFS_CGNS_VERSION}")
