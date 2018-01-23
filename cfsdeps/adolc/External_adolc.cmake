#-------------------------------------------------------------------------------
# The ADOL-C Library
#
# Project Homepage
# https://gitlab.com/adol-c/adol-c
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to gmp sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(adolc_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/adolc")
set(adolc_source  "${adolc_prefix}/src/adolc")
set(adolc_install "${adolc_prefix}/install")

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------

# here we need the necessary patches for cfs

#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gmp/gmp-patch.cmake.in")
#SET(PFN "${gmp_prefix}/gmp-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "https://gitlab.com/adol-c/adol-c/repository/${ADOLC_VER}/archive.tar.gz"
  "${ADOLC_URL}/${ADOLC_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/adolc/${ADOLC_GZ}")
SET(MD5_SUM ${ADOLC_MD5})

SET(DLFN "${adolc_prefix}/adolc-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/adolc/adolc-post_install.cmake.in")
SET(PI "${adolc_prefix}/adolc-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "adolc" "${ADOLC_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${adolc_install}")

SET(ZIPFROMCACHE "${adolc_prefix}/adolc-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${adolc_prefix}/adolc-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of ADOL-C libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
# ToDo: needs to be the correct link line
SET(ADOLC_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "ADOL-C library.")

SET(ADOLC_INCLUDE_DIR "${adolc_install}/include/adolc")

MARK_AS_ADVANCED(ADOLC_LIBRARY)
MARK_AS_ADVANCED(ADOLC_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The ADOL-C external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(adolc
    PREFIX "${adolc_prefix}"
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
    set(ADOLC_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "program to build ADOL-C")
  endif("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  # check for configure dependencies
  
  # add external project 
  MARK_AS_ADVANCED(ADOLC_MAKE_PROGRAM)
  ExternalProject_Add(adolc
#    DEPENDS hdf5-static
    PREFIX "${adolc_prefix}"
    SOURCE_DIR "${adolc_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${GMP_MD5}
    PATCH_COMMAND "" #${CMAKE_COMMAND} -P "${PFN}"
    BINARY_DIR "${adolc_source}"
    CONFIGURE_COMMAND "${adolc_source}/configure" --prefix=${adolc_install} --with-cflags=${CFS_C_FLAGS} --with-cxxflags=${CFS_CXX_FLAGS} # don't know if we need FLAGS from CFS -> they lead to a lot of warnings ...
    BUILD_COMMAND ${ADOLC_MAKE_PROGRAM} -f Makefile
    INSTALL_COMMAND ${ADOLC_MAKE_PROGRAM} -f Makefile install
    BUILD_BYPRODUCTS ${ADOLC_LIBRARY}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(adolc cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY "${adolc_prefix}"
  )
  
  #-------------------------------------------------------------------------------
  # Execute the stuff from gmp-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(adolc post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(adolc cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY "${CFS_BINARY_DIR}"
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  adolc
)
