#-------------------------------------------------------------------------------
# bzip2  is a  freely available,  patent free  (see below),  high-quality data
# compressor. It typically  compresses files to within 10% to  15% of the best
# available  techniques (the  PPM family  of statistical  compressors), whilst
# being  around  twice  as  fast  at  compression  and  six  times  faster  at
# decompression. 
#
# Needed by Boost.
#
# Project Homepage
# http://www.bzip.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to bzip2 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(bzip2_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/bzip2")
set(bzip2_source  "${bzip2_prefix}/src/bzip2")
set(bzip2_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${bzip2_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_C_FLAGS:FILEPATH=${CFLAGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
  )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/bzip2/bzip2-patch.cmake.in")
SET(PFN "${bzip2_prefix}/bzip2-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp3.de.freebsd.org/FreeBSD/ports/distfiles/${BZIP2_GZ}"
  "ftp://ftp.jussieu.fr/pub/haiku/releases/r1alpha4/sources/${BZIP2_GZ}"
  "${BZIP2_URL}/${BZIP2_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/bzip2/${BZIP2_GZ}")
SET(MD5_SUM ${BZIP2_MD5})

SET(DLFN "${bzip2_prefix}/bzip2-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "bzip2_${BZIP2_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}_${CMAKE_BUILD_TYPE}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "bzip2_${BZIP2_VER}_${CFS_ARCH_STR}_${CFS_CXX_COMPILER_NAME}_${CFS_CXX_COMPILER_VER}_${CMAKE_BUILD_TYPE}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${bzip2_prefix}")

SET(ZIPFROMCACHE "${bzip2_prefix}/bzip2-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${bzip2_prefix}/bzip2-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The bzip2 external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(bzip2-static
    PREFIX "${bzip2_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  ExternalProject_Add(bzip2-shared
    PREFIX "${bzip2_prefix}"
    DOWNLOAD_COMMAND ""
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project (static)
  #-------------------------------------------------------------------------------
  ExternalProject_Add(bzip2-static
    PREFIX "${bzip2_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/bzip2
    SOURCE_DIR ${bzip2_source}
    URL ${BZIP2_URL}/${BZIP2_GZ}
    URL_MD5 ${BZIP2_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=OFF
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(bzip2-static cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
     DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${bzip2_prefix}
  )
  
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project (shared)
  #-------------------------------------------------------------------------------
  ExternalProject_Add(bzip2-shared
    DEPENDS bzip2-static
    PREFIX "${bzip2_prefix}"
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${bzip2_source}
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=ON
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(bzip2-shared cfsdeps_zipToCache
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
  bzip2-static
  bzip2-shared
)

SET(BZIP2_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of BZIP2 libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(BZIP2_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "BZIP2 library.")

MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR)
MARK_AS_ADVANCED(BZIP2_LIBRARY)
