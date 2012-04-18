#-------------------------------------------------------------------------------
# Multiple precision floating-point computations with correct rounding
# Needed by CGAL
#
# Project Homepage
# http://mpfrlib.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to mpfr sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(mpfr_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/mpfr")
set(mpfr_source  "${mpfr_prefix}/src/mpfr")
set(mpfr_install  "${CMAKE_CURRENT_BINARY_DIR}")


if(CFS_ARCH STREQUAL "I386")
  SET(ABI "32")
else(CFS_ARCH STREQUAL "I386")
  SET(ABI "64")
endif(CFS_ARCH STREQUAL "I386")

SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/mpfr/mpfr-configure.sh.in")
SET(CONF "${mpfr_prefix}/mpfr-configure.sh")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# The mpfr external project
#-------------------------------------------------------------------------------
ExternalProject_Add(mpfr
  DEPENDS gmp
  PREFIX "${mpfr_prefix}"
  SOURCE_DIR "${mpfr_source}"
  URL ${MPFR_URL}/${MPFR_BZ2}
  URL_MD5 ${MPFR_MD5}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND sh ${CONF}
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile install
)

SET(MPFR_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of MPFR libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(MPFR_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}mpfr${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "MPFR library.")

MARK_AS_ADVANCED(MPFR_LIBRARY)
MARK_AS_ADVANCED(MPFR_INCLUDE_DIR)

