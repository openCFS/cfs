SET(CGAL_FOUND 0)

#-------------------------------------------------------------------------------
# Determine path of CGAL library.
#-------------------------------------------------------------------------------
SET(GMP_URL "${LSE17_SOURCES_DIR}/gmp")
SET(GMP_BZ2 "gmp-4.2.4.tar.bz2")
SET(GMP_MD5 "fc1e3b3a2a5038d4d74138d0b9cf8dbe")
INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/gmp/External_gmp.cmake")

SET(MPFR_URL "${LSE17_SOURCES_DIR}/mpfr")
SET(MPFR_BZ2 "mpfr-2.4.0.tar.bz2")
SET(MPFR_MD5 "f5916d785d4f7e7282057f6a3ebff9ce")
INCLUDE("${CFS_SOURCE_DIR}/cfsdeps/mpfr/External_mpfr.cmake")

BUILD_EXTLIB("CGAL"
  "${CFS_BINARY_DIR}/include/CGAL"
  "${CFS_DEPS_ROOT}/cgal/build_cgal.pl"
  "build_cgal.log")

#-------------------------------------------------------------------------------
# Mark path of CGAL library as advanced.
#-------------------------------------------------------------------------------
SET(CGAL_LIBRARY_DEBUG
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libCGAL.a"
  CACHE FILEPATH "CGAL debug library.")

SET(CGAL_LIBRARY_RELEASE
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libCGAL.a"
  CACHE FILEPATH "CGAL release library.")

MARK_AS_ADVANCED(CGAL_LIBRARY_DEBUG)
MARK_AS_ADVANCED(CGAL_LIBRARY_RELEASE)

IF(DEBUG)
  SET(CGAL_LIBRARY "${CGAL_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(CGAL_LIBRARY "${CGAL_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(CGAL_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
SET(CGAL_FOUND 1)

IF(CGAL_FOUND)

  #-----------------------------------------------------------------------------
  # Determine version of CGAL by preprocessing its version header
  #-----------------------------------------------------------------------------
  SET(CGAL_TEST_FILE "${CFS_BINARY_DIR}/include/get_cgal_infos.hh")
  FILE(WRITE ${CGAL_TEST_FILE}
    "#include <${CGAL_INCLUDE_DIR}/CGAL/version.h>\n\n"
    "<CFS_CGAL_VERSION>CGAL_VERSION</CFS_CGAL_VERSION>\n")

  SET(CGAL_VERSION_CMD "${CMAKE_CXX_COMPILER} -E ${CGAL_TEST_FILE}")

  EXEC_PROGRAM(${CGAL_VERSION_CMD}
    ARGS
    OUTPUT_VARIABLE CGAL_VERSION_STR
    RETURN_VALUE RETVAL)

#  MESSAGE("CGAL_VERSION_STR ${CGAL_VERSION_STR}")
  FILE(REMOVE ${CGAL_TEST_FILE})

  STRING(REGEX MATCH
    "<CFS_CGAL_VERSION>.*</CFS_CGAL_VERSION>"
    CFS_CGAL_VERSION "${CGAL_VERSION_STR}")
  STRING(REGEX MATCH "[0-9]+\\.[0-9]+"
    CFS_CGAL_VERSION "${CFS_CGAL_VERSION}")

#  MESSAGE("CFS_CGAL_VERSION ${CFS_CGAL_VERSION}")

ENDIF(CGAL_FOUND)
