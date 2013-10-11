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

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/mpfr/mpfr-configure.cmake.in")
SET(CONF "${mpfr_prefix}/mpfr-configure.cmake")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since lse17 may not be
# accessible from behind firewalls.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://ftp.uni-erlangen.de/mirrors/GNU/mpfr/${MPFR_BZ2}"
  "http://ftp.gwdg.de/pub/misc/gnu/ftp/gnu/mpfr/${MPFR_BZ2}"
  "http://ftp.gnu.org/gnu/mpfr/${MPFR_BZ2}"
)

#-------------------------------------------------------------------------------
# Try to download sources into CFSDEPS cache directory.
#-------------------------------------------------------------------------------
DOWNLOAD_CFSDEPS(
  "${CFS_DEPS_CACHE_DIR}/sources/mpfr/${MPFR_BZ2}"
  ${MPFR_MD5}
  "${MIRRORS}"
)

#-------------------------------------------------------------------------------
# The mpfr external project
#-------------------------------------------------------------------------------
ExternalProject_Add(mpfr
  DEPENDS gmp
  PREFIX "${mpfr_prefix}"
  SOURCE_DIR "${mpfr_source}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/mpfr
  URL ${MPFR_URL}/${MPFR_BZ2}
  URL_MD5 ${MPFR_MD5}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile install
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  mpfr
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

