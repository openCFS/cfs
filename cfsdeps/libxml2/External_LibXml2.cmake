# libxml2 is an alternative for xerces. xerces and libxml2 seem to be the only
# xml parser libs able to validate schema and apply schema defaults.
set(libxml2_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/libxml2")
set(libxml2_source  "${libxml2_prefix}/src/libxml2")
set(libxml2_install  "${libxml2_prefix}/libxml2_target")

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://xmlsoft.org/libxml2/${LIBXML2_GZ}"
  "http://ftp.osuosl.org/pub/blfs/conglomeration/libxml2/${LIBXML2_GZ}"
  "${libxml2_URL}/${LIBXML2_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/libxml2/${LIBXML2_GZ}")
SET(MD5_SUM "${LIBXML2_MD5}")

SET(DLFN "${libxml2_prefix}/libxml2-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)
  
# see below for what we need the post_install (pre-zip)
SET(PI "${libxml2_prefix}/ipopt-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/libxml2/libxml2-post_install.cmake.in" "${PI}" @ONLY) 
  
  
# using no cmake build we make do difference between debug and release build 
PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "libxml2" "${LIBXML2_VER}")

# not using cmake build (see above) we don't use manifuest and copy the   
SET(TMP_DIR "${libxml2_install}")

SET(ZIPFROMCACHE "${libxml2_prefix}/libxml2-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${libxml2_prefix}/libxml2-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The libxml2 external project.
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(libxml2
    PREFIX "${libxml2_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE()
  #-------------------------------------------------------------------------------
  # we make a standard configure build in a TMP_DIR target. ZIPTOCACHE uses this and
  # copies lib to lib64/BUILD.
  #-------------------------------------------------------------------------------
  ExternalProject_Add(libxml2
    PREFIX "${libxml2_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/libxml2
    URL ${LOCAL_FILE}
    URL_MD5 ${libxml2_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND "" # at the moment only make and no cmake version ${CMAKE_COMMAND} -P "${PFN}" 
    CONFIGURE_COMMAND ${libxml2_source}/configure --prefix=${libxml2_install} --disable-shared --without-ftp --without-html --without-http --without-icu --without-iconv --without-python --without-modules --without-lzma 
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(libxml2 cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${libxml2_prefix}
  )
  
  # libxml2 make install creates an include/libxml2/libxml but includes <libxml/**.h>
  # move the stuff to make it easy
  # Note that without precompiled packages this step and the lib copying in cfsdeps_zipToCache
  # will be missing. Fix it when necessary!  
  ExternalProject_Add_Step(libxml2 post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
   IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(libxml2 cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS ${CFSDEPS} libxml2)

#-------------------------------------------------------------------------------
# Determine paths of libxml2 libraries.
# There is no separation between BLAS and LAPACK, both are combined. While BLAS is optimized
# only parts of LAPACK are optimized. libxml2 can be build without LAPACK in request.
#-------------------------------------------------------------------------------
SET(LIBXML2_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libxml2.a" CACHE FILEPATH "libxml2 library.")

MARK_AS_ADVANCED(LIBXML2_LIBRARY)

SET(LIBXML2_INCLUDE_DIR "{CFS_BINARY_DIR}/include/libxml2")




