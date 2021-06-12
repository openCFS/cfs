#-------------------------------------------------------------------------------
# muParser - a fast math parser library
#
# Project Homepage
# http://muparser.sourceforge.net/
#-------------------------------------------------------------------------------
set(muparser_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser")
set(muparser_source  "${muparser_prefix}/src/muparser")
# we need no muparser_install as we build directly to CMAKE_CURRENT_BINARY_DIR and select the files to zip from the manifest

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
# https://github.com/beltoforion/muparser/archive/v2.2.6.1.tar.gz
SET(MIRRORS
  "https://github.com/beltoforion/muparser/archive/${MUPARSER_TGZ}"
  "${MUPARSER_URL}/${MUPARSER_TGZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/muparser/${MUPARSER_TGZ}")
SET(MD5_SUM ${MUPARSER_MD5})

SET(DLFN "${muparser_prefix}/muparser-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/muparser/license/" DESTINATION "${CFS_BINARY_DIR}/license/muparser" )



# do make a difference between debug and release build since we are using cmake now 
PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "muparser" "${MUPARSER_VER}")

# we need to set TMP_DIR for ZIP_TO_CACHE, read in cfsdeps_zipToCache.cmake.in such that ZIP_TO_CACHE finds  cfsdeps/muparser/src/muparser/install_manifest.txt
SET(TMP_DIR "${muparser_prefix}")

SET(ZIPFROMCACHE "${muparser_prefix}/muparser-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${muparser_prefix}/muparser-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The muparser external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(muparser
    PREFIX "${muparser_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE()
  # on mac with clang omp.h is not found. This is no issue for cfs itself
  # a very dirty hack is to do ln -s /usr/local/include/omp.h in cfsdeps/muparser/src/muparser/include
 
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  SET(MUPARSER_SHARED_LIBS OFF)
  IF(WIN32)
    SET(MUPARSER_SHARED_LIBS ON)
    IF(DEBUG)
      SET(CMAKE_ARGS "${CMAKE_ARGS}" "-DCMAKE_BUILD_TYPE=DEBUG")
    ENDIF()
  ENDIF()
  ExternalProject_Add(muparser
    PREFIX "${muparser_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/muparser
    URL ${LOCAL_FILE}
    URL_MD5 ${MUPARSER_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CMAKE_ARGS
      ${CMAKE_ARGS}
      # install to cfs, we collect our stuff for precompiled from the manifest
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
      -DBUILD_SHARED_LIBS:BOOL=${MUPARSER_SHARED_LIBS}
      -DBUILD_TESTING:BOOL=OFF
      -DENABLE_OPENMP:BOOL=OFF # makes problems when precompiled is parallel and cfs is serial
      -DENABLE_SAMPLES:BOOL=OFF
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(muparser cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${muparser_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(muparser cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS
SET(CFSDEPS ${CFSDEPS} muparser)

# Determine paths of MUPARSER libraries.
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
IF(WIN32)
  SET(MUPARSER_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/muparser.lib" CACHE FILEPATH "MUPARSER library.")
ELSE(WIN32)
  SET(MUPARSER_LIBRARY "${LD}/libmuparser.a" CACHE FILEPATH "MUPARSER library.")
ENDIF(WIN32)

SET(MUPARSER_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(MUPARSER_INCLUDE_DIR)
MARK_AS_ADVANCED(MUPARSER_LIBRARY)
