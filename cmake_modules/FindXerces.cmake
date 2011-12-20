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
BUILD_EXTLIB("Xerces-C"
  "${CFS_BINARY_DIR}/include/xercesc"
  "${CFS_DEPS_ROOT}/xerces/build_xerces.pl"
  "build_xercesc.log")

#-------------------------------------------------------------------------------
# Determine paths of XERCES libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(XERCES_LIBRARY_RELEASE ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/xerces-c_2.lib)
  SET(XERCES_LIBRARY_DEBUG ${CFS_BINARY_DIR}/lib/${CFS_ARCH_STR}/xerces-c_2d.lib)
ELSE(WIN32)
  SET(XERCES_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libxerces-c.a")

  IF(CFS_DISTRO STREQUAL "MACOSX")
    SET(XERCES_LIBRARY "${XERCES_LIBRARY};-framework CoreFoundation;-framework CoreServices")
  ENDIF(CFS_DISTRO STREQUAL "MACOSX")

  SET(XERCES_LIBRARY_DEBUG
    "${XERCES_LIBRARY}"
    CACHE FILEPATH "Xerces-C debug library.")

  SET(XERCES_LIBRARY_RELEASE
    "${XERCES_LIBRARY}"
    CACHE FILEPATH "Xerces-C release library.")
ENDIF(WIN32)

MARK_AS_ADVANCED(XERCES_LIBRARY_DEBUG)
MARK_AS_ADVANCED(XERCES_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set XERCES_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(XERCES_LIBRARY "${XERCES_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(XERCES_LIBRARY "${XERCES_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(XERCES_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

SET(XERCES_FOUND 1)

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

