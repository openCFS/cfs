#-------------------------------------------------------------------------------
# 
# OpenBLAS is a free and efficient and parallel BLAS/ LAPACK implementation.
# It is a successor of the stalled GOTO BLAS
# http://www.openblas.net/
# 
#-------------------------------------------------------------------------------
set(openblas_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/openblas")
set(openblas_source  "${openblas_prefix}/src/openblas")
# we need no openblas_install as we build directly to CMAKE_CURRENT_BINARY_DIR and select the files to zip from the manifest

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "https://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_GZ}"
  "${CFS_DS_SOURCES_DIR}/openblas/${OPENBLAS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/openblas/${OPENBLAS_GZ}")
SET(MD5_SUM ${OPENBLAS_MD5})

SET(DLFN "${openblas_prefix}/openblas-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/openblas/license/" DESTINATION "${CFS_BINARY_DIR}/license/openblas" )



# do make a difference between debug and release build since we are using cmake now 
PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "openblas" "${OPENBLAS_VER}")

# we need to set TMP_DIR for ZIP_TO_CACHE, read in cfsdeps_zipToCache.cmake.in such that ZIP_TO_CACHE finds  cfsdeps/openblas/src/openblas/install_manifest.txt
SET(TMP_DIR "${openblas_prefix}")

SET(ZIPFROMCACHE "${openblas_prefix}/openblas-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${openblas_prefix}/openblas-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The OpenBLAS external project.
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(openblas
    PREFIX "${openblas_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  # CFS_LIB_TARGET from the patched Makefile.install
  # NO_SHARED=1 is an openblas option to omit shared libraries. The standard should be with shared
  # but when done via this external project the shared libs are not built?! Therefore NO_SHARED
  # in the install command
  #-------------------------------------------------------------------------------
  ExternalProject_Add(openblas
    PREFIX "${openblas_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/openblas
    URL ${LOCAL_FILE}
    URL_MD5 ${OPENBLAS_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}
      -DCMAKE_INSTALL_LIBDIR:PATH=${CMAKE_CURRENT_BINARY_DIR}/lib
   )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(openblas cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${openblas_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(openblas cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
set(CFSDEPS ${CFSDEPS} openblas)

#-------------------------------------------------------------------------------
# Determine paths of OpenBLAS libraries.
# There is no separation between BLAS and LAPACK, both are combined. While BLAS is optimized
# only parts of LAPACK are optimized. OpenBLAS can be build without LAPACK in request but then 
# cfs does not link
#-------------------------------------------------------------------------------
SET(BLAS_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/libopenblas.a;-lpthread;-lm")
SET(LAPACK_LIB "${BLAS_LIB}")

set(OPENBLAS_LIBRARY ${BLAS_LIB} CACHE FILEPATH "OpenBLAS library.")
MARK_AS_ADVANCED(OPENBLAS_LIBRARY)

# for OPENBLAS LAPACK_LIBRARY is the same as BLAS_LIBRARY. See also External_LAPACK for netlib and FindIntelMKL
# e.g. lib/libopenblas.a;-lpthread
if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  set(BLAS_LIBRARY "${OPENBLAS_LIBRARY}")
endif()
