#-------------------------------------------------------------------------------
# Fastest Fourier Transform in the West (C subroutine library)
#
# Project Homepage
# http://www.fftw.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to fftw sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(fftw_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/fftw")
set(fftw_source  "${fftw_prefix}/src/fftw")
set(fftw_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/fftw/fftw-configure.cmake.in")
SET(CONF "${fftw_prefix}/fftw-configure.cmake")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.fftw.org/${FFTW_GZ}"
  "ftp://ftp.fftw.org/pub/fftw/${FFTW_GZ}"
  "${FFTW_URL}/${FFTW_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/fftw/${FFTW_GZ}")
SET(MD5_SUM ${FFTW_MD5})

SET(DLFN "${fftw_prefix}/fftw-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "fftw" "${FFTW_VER}") 
  
# This should be either PREFIX_DIR/src (install manifest is used for zipping)
# or PREFIX_DIR/install (install directory will be zipped)
SET(TMP_DIR "${fftw_prefix}/src")

SET(ZIPFROMCACHE "${fftw_prefix}/fftw-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${fftw_prefix}/fftw-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The fftw external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(fftw
    PREFIX "${fftw_prefix}"
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
  ExternalProject_Add(fftw
    PREFIX "${fftw_prefix}"
    SOURCE_DIR "${fftw_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${FFTW_MD5}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
    BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile install
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(fftw cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${mpfr_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(fftw cfsdeps_zipToCache
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
  fftw
)

SET(FFTW_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of FFTW libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
IF(WIN32)
  SET(FFTW_LIBRARY
    "${LD}/${CMAKE_IMPORT_LIBRARY_PREFIX}fftw3${CMAKE_IMPORT_LIBRARY_SUFFIX}"
    CACHE FILEPATH "FFTW library.")
ELSE(WIN32)
  SET(FFTW_LIBRARY
      "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}fftw3${CMAKE_STATIC_LIBRARY_SUFFIX}"
      CACHE FILEPATH "FFTW library.")
ENDIF(WIN32)

MARK_AS_ADVANCED(FFTW_LIBRARY)
MARK_AS_ADVANCED(FFTW_INCLUDE_DIR)
