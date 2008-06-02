# - Find Xerces-C library.
# This module finds if Xerces is installed and determines where the
# library is. It also determines what the name of the library is.
# This code sets the following variables:
#  XERCES_LIBRARY        = path to XERCES library (xerces-c)
#  XERCES_INCLUDE_DIR    = path to XERCES headers
#  XERCES_FOUND          = bool which says if BLAS has been found
#  XERCES_VERSION_MAJOR
#  XERCES_VERSION_MINOR
#  XERCES_VERSION_REVISION
#  XERCES_VERSION

SET(XERCES_FOUND 0)

#-------------------------------------------------------------------------------
# Determine paths of XERCES libraries.
#-------------------------------------------------------------------------------
SET (XERCES_POSSIBLE_LIB_PATHS
  ${CFSDEPS_LIBRARY_DIR}/xercesc_2.8.0
  ${CFSDEPS_LIBRARY_DIR}
  /usr/lib64
  /usr/lib
  /usr/local/lib64
  /usr/local/lib
)

FIND_LIBRARY(XERCES_LIBRARY
  NAMES xerces-c
  PATHS ${XERCES_POSSIBLE_LIB_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

#-------------------------------------------------------------------------------
# Mark paths of XERCES libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(XERCES_LIBRARY)


#-------------------------------------------------------------------------------
# Look for Xerces header.
#-------------------------------------------------------------------------------
SET (XERCES_POSSIBLE_INCLUDE_PATHS
  ${CFSDEPS_INCLUDE_DIR}/xercesc_2.8.0
  ${CFSDEPS_INCLUDE_DIR}
  /usr/include
  /usr/local/include
  )

FIND_PATH(XERCES_INCLUDE_DIR
  NAMES xercesc/util/XercesVersion.hpp 
  PATHS ${XERCES_POSSIBLE_INCLUDE_PATHS}
  NO_DEFAULT_PATH
  NO_CMAKE_ENVIRONMENT_PATH
  NO_CMAKE_PATH
  NO_SYSTEM_ENVIRONMENT_PATH
  NO_CMAKE_SYSTEM_PATH
  )

MARK_AS_ADVANCED(XERCES_INCLUDE_DIR)


IF(XERCES_LIBRARY AND XERCES_INCLUDE_DIR)
  SET(XERCES_FOUND 1)
ENDIF(XERCES_LIBRARY AND XERCES_INCLUDE_DIR)


IF(XERCES_FOUND)

  SET(XERCES_TEST_FILE "${CFS_BINARY_DIR}/include/get_xerces_infos.hh")
  FILE(WRITE ${XERCES_TEST_FILE}
    "#include <${XERCES_INCLUDE_DIR}/xercesc/util/XercesVersion.hpp>\n\n"
    "<CFS_XERCES_VERSION>XERCES_VERSION_MAJOR.XERCES_VERSION_MINOR.XERCES_VERSION_REVISION</CFS_XERCES_VERSION>\n")

  IF(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(XERCES_VERSION_CMD "${CFS_CXX_COMPILER} /nologo /EP ${XERCES_TEST_FILE}")
  ELSE(CMAKE_VC_COMPILER_TESTS_RUN)
    SET(XERCES_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${XERCES_TEST_FILE}")
  ENDIF(CMAKE_VC_COMPILER_TESTS_RUN)

  EXEC_PROGRAM(${XERCES_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE XERCES_VERSION_STR
    RETURN_VALUE RETVAL)

  FILE(REMOVE ${XERCES_TEST_FILE})

  STRING(REGEX MATCH "<CFS_XERCES_VERSION>.*</CFS_XERCES_VERSION>" CFS_XERCES_VERSION "${XERCES_VERSION_STR}")
  STRING(REGEX REPLACE " " "" CFS_XERCES_VERSION "${CFS_XERCES_VERSION}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" CFS_XERCES_VERSION "${CFS_XERCES_VERSION}")

#  MESSAGE("CFS_XERCES_VERSION ${CFS_XERCES_VERSION}")

  STRING(REGEX REPLACE "xerces-c\\.so" "xerces-c.a"
    XERCES_LIBRARY "${XERCES_LIBRARY}")

  SET(XERCES_LIBRARY "${XERCES_LIBRARY}" CACHE PATH "Path to Xerces-C library." FORCE)

ENDIF(XERCES_FOUND)

IF(XERCES_LIBRARY)
  IF(XERCES_INCLUDE_DIR)
    SET(XERCES_FOUND 1)
  ENDIF(XERCES_INCLUDE_DIR)
ENDIF(XERCES_LIBRARY)

MARK_AS_ADVANCED(
  XERCES_LIBRARY
  XERCES_INCLUDE_DIR
  XERCES_FOUND
  XERCES_VERSION_MAJOR
  XERCES_VERSION_MINOR
  XERCES_VERSION_REVISION
  XERCES_VERSION
  )

IF(NOT XERCES_FOUND)
  MESSAGE("Warning: XERCES could not be found! Please specify proper paths.")
ENDIF(NOT XERCES_FOUND)

