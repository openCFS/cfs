#-------------------------------------------------------------------------------
# A Massively Spiffy Yet Delicately Unobtrusive Compression Library
# Needed by Botan, HDF5 and libxml2.
#
# Project Homepage
# http://www.zlib.net/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to zlib sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(zlib_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zlib")
set(zlib_source  "${zlib_prefix}/src/zlib")
set(zlib_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/zlib/zlib-patch.cmake.in")
SET(PFN "${zlib_prefix}/zlib-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${zlib_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

IF(MINGW)
  LIST(APPEND CMAKE_ARGS
    -DCFS_ARCH:STRING=${CFS_ARCH}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/zlib/zlib-patch.cmake.in")
SET(PFN "${zlib_prefix}/zlib-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://zlib.net/${ZLIB_GZ}"
  "http://fossies.org/linux/misc/${ZLIB_GZ}"
  "${ZLIB_URL}/${ZLIB_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/zlib/${ZLIB_GZ}")
SET(MD5_SUM ${ZLIB_MD5})

SET(DLFN "${zlib_prefix}/zlib-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "zlib" "${ZLIB_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${zlib_prefix}")

SET(ZIPFROMCACHE "${zlib_prefix}/zlib-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${zlib_prefix}/zlib-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The zlib external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(zlib
    PREFIX "${zlib_prefix}"
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
  ExternalProject_Add(zlib
    PREFIX ${zlib_prefix}
    SOURCE_DIR ${zlib_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${ZLIB_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(zlib cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${zlib_prefix}
  )

  #-------------------------------------------------------------------------------
  # No zip to cache here but after minizip is built. See below.
  #-------------------------------------------------------------------------------
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
LIST(APPEND CFSDEPS zlib)

IF(MINGW)
  SET(ZLIB_LIB zlibstatic)
  SET(ZLIB_SHARED_LIB zlib)
ELSE(MINGW)
  IF(UNIX)
    SET(ZLIB_LIB z)
    SET(ZLIB_SHARED_LIB z)
  ELSE(UNIX)
    SET(ZLIB_LIB zlibstatic)
    SET(ZLIB_SHARED_LIB zlib)
    IF(DEBUG)
      SET(ZLIB_LIB "${ZLIB_LIB}d")
      SET(ZLIB_SHARED_LIB "${ZLIB_SHARED_LIB}d")
    ENDIF()
  ENDIF(UNIX)
ENDIF(MINGW)

SET(LD ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR})
SET(ZLIB_LIBRARY
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${ZLIB_LIB}${CMAKE_STATIC_LIBRARY_SUFFIX}
  CACHE FILEPATH "zlib library")
SET(ZLIB_SHARED_LIBRARY
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${ZLIB_SHARED_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
IF(MINGW)
  SET(ZLIB_SHARED_LIBRARY "${ZLIB_SHARED_LIBRARY}.a")
ENDIF(MINGW)
SET(ZLIB_SHARED_LIBRARY ${ZLIB_SHARED_LIBRARY} CACHE FILEPATH "zlib shared library")
SET(ZLIB_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "zlib include directory")

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_SHARED_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The minizip external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # Do nothing. The project is included in zlib, so we already unzipped the files!
  #-------------------------------------------------------------------------------
  ExternalProject_Add(minizip
    PREFIX "${zlib_prefix}"
    DOWNLOAD_COMMAND ""
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
  ExternalProject_Add(minizip
    DEPENDS zlib
    PREFIX ${zlib_prefix}
    SOURCE_DIR ${zlib_source}/contrib/minizip
    DOWNLOAD_COMMAND ""
    PATCH_COMMAND ""
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:PATH=${ZLIB_SHARED_LIBRARY}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(minizip cfsdeps_zipToCache
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
LIST(APPEND CFSDEPS minizip)

SET(MINIZIP_SHARED_LIB minizip)
IF(MINGW OR WIN32)
  SET(MINIZIP_SHARED_LIB minizipdll)
ENDIF()

SET(MINIZIP_LIBRARY 
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}minizip_static${CMAKE_STATIC_LIBRARY_SUFFIX}
  CACHE FILEPATH "minizip library")
SET(MINIZIP_SHARED_LIBRARY
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${MINIZIP_SHARED_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
IF(MINGW)
  SET(MINIZIP_SHARED_LIBRARY "${MINIZIP_SHARED_LIBRARY}.a")
ENDIF(MINGW)
SET(MINIZIP_SHARED_LIBRARY ${MINIZIP_SHARED_LIBRARY}
  CACHE FILEPATH "minizip shared library")
SET(MINIZIP_INCLUDE_DIR ${CFS_BINARY_DIR}/include/minizip
  CACHE PATH "minizip include directory")

MARK_AS_ADVANCED(MINIZIP_LIBRARY)
MARK_AS_ADVANCED(MINIZIP_SHARED_LIBRARY)
MARK_AS_ADVANCED(MINIZIP_INCLUDE_DIR)

# Determine version of Minizip by reading its header.
IF(EXISTS "${MINIZIP_INCLUDE_DIR}/zip.h")
  FILE(STRINGS "${MINIZIP_INCLUDE_DIR}/zip.h" MINIZIP_VERSION REGEX "Version [0-9]")
  STRING(STRIP "${MINIZIP_VERSION}" MINIZIP_VERSION)
ELSE()
  SET(MINIZIP_VERSION "N/A")
ENDIF()
