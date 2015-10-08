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
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://ftp.uni-erlangen.de/mirrors/GNU/mpfr/${MPFR_BZ2}"
  "http://ftp.gwdg.de/pub/misc/gnu/ftp/gnu/mpfr/${MPFR_BZ2}"
  "http://ftp.gnu.org/gnu/mpfr/${MPFR_BZ2}"
  "${MPFR_URL}/${MPFR_BZ2}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/mpfr/${MPFR_BZ2}")
SET(MD5_SUM ${MPFR_MD5})

SET(DLFN "${mpfr_prefix}/mpfr-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "mpfr" "${MPFR_VER}")  
  
# This should be either PREFIX_DIR/src (install manifest is used for zipping)
# or PREFIX_DIR/install (install directory will be zipped)
SET(TMP_DIR "${mpfr_prefix}/src")

SET(ZIPFROMCACHE "${mpfr_prefix}/mpfr-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${mpfr_prefix}/mpfr-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The mpfr external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(mpfr
    PREFIX "${mpfr_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(mpfr
    DEPENDS gmp
    PREFIX "${mpfr_prefix}"
    SOURCE_DIR "${mpfr_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${MPFR_MD5}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile install
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(mpfr cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${mpfr_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(mpfr cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

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

