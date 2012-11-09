SET(FFTW_FOUND 0)

#-------------------------------------------------------------------------------
# Build FFTW 
#-------------------------------------------------------------------------------
BUILD_EXTLIB("fftw"
  "${CFS_BINARY_DIR}/include/fftw3.h"
  "${CFS_DEPS_ROOT}/fftw/build_fftw.pl"
  "build_fftw.log")

SET(FFTW_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of FFTW libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  FIND_FILE(FFTW_LIBRARY_DEBUG
    NAMES fftw.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_FILE(FFTW_LIBRARY_RELEASE
    NAMES fftw.lib
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ELSE(WIN32)
  FIND_LIBRARY(FFTW_LIBRARY_DEBUG
    NAMES fftw3
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
  FIND_LIBRARY(FFTW_LIBRARY_RELEASE
    NAMES fftw3 
    PATHS ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
    )
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of FFTW libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(FFTW_LIBRARY_DEBUG)
MARK_AS_ADVANCED(FFTW_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set FFTW_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(FFTW_LIBRARY "${FFTW_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(FFTW_LIBRARY "${FFTW_LIBRARY_RELEASE}")
ENDIF(DEBUG)

# Safety check for existance of libs and headers
IF(EXISTS "${FFTW_LIBRARY}"
    AND EXISTS "${FFTW_INCLUDE_DIR}/fftw3.h")
 
  SET(FFTW_FOUND 1)
  
ELSE(EXISTS "${FFTW_LIBRARY}"
    AND EXISTS "${FFTW_INCLUDE_DIR}/fftw3.h")
  
  MESSAGE(FATAL_ERROR "Could not find ${FFTW_LIBRARY} or ${FFTW_INCLUDE_DIR}/fftw3.h")
  
ENDIF(EXISTS "${FFTW_LIBRARY}"
  AND EXISTS "${FFTW_INCLUDE_DIR}/fftw3.h")

#-----------------------------------------------------------------------------
# Determine version of FFTW by preprocessing its version header OSBOLETE!
#-----------------------------------------------------------------------------
#SET(FFTW_TEST_FILE "${CFS_BINARY_DIR}/include/get_metis_infos.hh")
#FILE(WRITE ${FFTW_TEST_FILE}
#  "#include <${FFTW_INCLUDE_DIR}/defs.h>\n\n"
#  "<CFS_FFTW_VERSION>FFTWTITLE</CFS_FFTW_VERSION>\n")

#FILE(TO_NATIVE_PATH ${FFTW_TEST_FILE} FFTW_TEST_FILE

#IF(CMAKE_VC_COMPILER_TESTS_RUN)
#  SET(FFTW_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP \"${FFTW_TEST_FILE}\"")
#ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
#  SET(FFTW_VERSION_CMD "${CMAKE_CXX_COMPILER} -E \"${FFTW_TEST_FILE}\"")
#ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)
#
#EXEC_PROGRAM(${FFTW_VERSION_CMD}
#  ARGS
#  OUTPUT_VARIABLE FFTW_VERSION_STR
#  RETURN_VALUE RETVAL)

#  MESSAGE("FFTW_VERSION_STR ${FFTW_VERSION_STR}")
#FILE(REMOVE ${FFTW_TEST_FILE})

#STRING(REGEX MATCH
#  "<CFS_FFTW_VERSION>.*</CFS_FFTW_VERSION>"
#  CFS_FFTW_VERSION "${FFTW_VERSION_STR}")
#STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+"
#  CFS_FFTW_VERSION "${CFS_FFTW_VERSION}")

# MESSAGE("CFS_FFTW_VERSION ${CFS_FFTW_VERSION}")