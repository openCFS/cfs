#-------------------------------------------------------------------------------
# The GNU MP Bignum Library
# Needed by MPFR and CGAL
#
# Project Homepage
# http://gmplib.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to gmp sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(gmp_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/gmp")
set(gmp_source  "${gmp_prefix}/src/gmp")
set(gmp_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-configure.cmake.in")
SET(CONF "${gmp_prefix}/gmp-configure.cmake")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-patch.cmake.in")
SET(PFN "${gmp_prefix}/gmp-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since lse17 may not be
# accessible from behind firewalls.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://ftp.uni-erlangen.de/mirrors/GNU/gmp/${GMP_BZ2}"
  "http://ftp.gwdg.de/pub/misc/gnu/ftp/gnu/gmp/${GMP_BZ2}"
  "http://ftp.gnu.org/gnu/gmp/${GMP_BZ2}"
)

#-------------------------------------------------------------------------------
# Try to download sources into CFSDEPS cache directory.
#-------------------------------------------------------------------------------
DOWNLOAD_CFSDEPS(
  "${CFS_DEPS_CACHE_DIR}/sources/gmp/${GMP_BZ2}"
  ${GMP_MD5}
  "${MIRRORS}"
)

#-------------------------------------------------------------------------------
# The gmp external project
#-------------------------------------------------------------------------------
ExternalProject_Add(gmp
  PREFIX "${gmp_prefix}"
  SOURCE_DIR "${gmp_source}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/gmp
  URL ${GMP_URL}/${GMP_BZ2}
  URL_MD5 ${GMP_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
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
  gmp
)

SET(GMP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of GMP libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(GMP_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "GMP library.")

SET(GMPXX_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "GMPXX library.")

MARK_AS_ADVANCED(GMP_LIBRARY)
MARK_AS_ADVANCED(GMPXX_LIBRARY)
MARK_AS_ADVANCED(GMP_INCLUDE_DIR)

