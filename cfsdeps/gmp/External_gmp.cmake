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
set(gmp_install  "${gmp_prefix}/install")

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
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://ftp.uni-erlangen.de/mirrors/GNU/gmp/${GMP_BZ2}"
  "http://ftp.gwdg.de/pub/misc/gnu/ftp/gnu/gmp/${GMP_BZ2}"
  "http://ftp.gnu.org/gnu/gmp/${GMP_BZ2}"
  "${GMP_URL}/${GMP_BZ2}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/gmp/${GMP_BZ2}")
SET(MD5_SUM ${GMP_MD5})

SET(DLFN "${gmp_prefix}/gmp-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-post_install.cmake.in")
SET(PI "${gmp_prefix}/gmp-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "gmp" "${GMP_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${gmp_install}")

SET(ZIPFROMCACHE "${gmp_prefix}/gmp-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${gmp_prefix}/gmp-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

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

#-------------------------------------------------------------------------------
# The gmp external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(gmp
    PREFIX "${gmp_prefix}"
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
  if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    # GMP does not use CMake but automake: We cannot use the CMake Ninja-generator,
    # we use make instead to build GMP
    find_program(GMP_MAKE_PROGRAM make)
  else("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    set(GMP_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM})
  endif("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  ExternalProject_Add(gmp
    PREFIX "${gmp_prefix}"
    SOURCE_DIR "${gmp_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${GMP_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
    BUILD_COMMAND ${GMP_MAKE_PROGRAM} -f Makefile
    INSTALL_COMMAND ${GMP_MAKE_PROGRAM} -f Makefile install
    BUILD_BYPRODUCTS ${GMP_LIBRARY} ${GMPXX_LIBRARY}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(gmp cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${gmp_prefix}
  )
  
  #-------------------------------------------------------------------------------
  # Execute the stuff from gmp-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(gmp post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(gmp cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
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
  gmp
)

SET(GMP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
