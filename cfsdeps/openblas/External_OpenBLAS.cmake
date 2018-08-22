#-------------------------------------------------------------------------------
# 
# OpenBLAS is a free and efficient and parallel BLAS/ LAPACK implementation.
# It is a successor of the stalled GOTO BLAS
# http://www.openblas.net/
# 
#-------------------------------------------------------------------------------
set(openblas_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/openblas")
set(openblas_source  "${openblas_prefix}/src/openblas")
set(openblas_install  "${CFS_BINARY_DIR}/tmp/openblas_target")

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://github.com/xianyi/OpenBLAS/archive/${OPENBLAS_GZ}"
  "${OPENBLAS_URL}/${OPENBLAS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/openblas/${OPENBLAS_GZ}")
SET(MD5_SUM ${OPENBLAS_MD5})

SET(DLFN "${openblas_prefix}/openblas-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)
  
#-------------------------------------------------------------------------------
# currently we don't use the openblas cmake interface, therefore the patch is not used
# Basically add the option CFS_TARGET_LIB to Makefile.install
#-------------------------------------------------------------------------------
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/openblas/openblas-patch.cmake.in")
#SET(PFN "${openblas_prefix}/openblas-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

# using no cmake build we make do difference between debug and release build 
PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "openblas" "${OPENBLAS_VER}")

# not using cmake build (see above) we don't use manifuest and copy the   
SET(TMP_DIR "${openblas_install}")

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
    PATCH_COMMAND "" # at the moment only make and no cmake version ${CMAKE_COMMAND} -P "${PFN}" 
    CONFIGURE_COMMAND ""
    BUILD_COMMAND  ${CMAKE_MAKE_PROGRAM} DYNAMIC_ARCH=1 libs netlib # netlib seems to add lapack which we also need, DYNAMIC_ARCH=1 creates binaries portable to diffrent CPU architectures 
    INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install NO_SHARED=1 "PREFIX=${openblas_install}" 
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
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS ${CFSDEPS} openblas)

#-------------------------------------------------------------------------------
# Determine paths of OpenBLAS libraries.
# There is no separation between BLAS and LAPACK, both are combined. While BLAS is optimized
# only parts of LAPACK are optimized. OpenBLAS can be build without LAPACK in request but then 
# cfs does not link
#-------------------------------------------------------------------------------
SET(BLAS_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libopenblas.a;-lpthread;-lm")
SET(LAPACK_LIB "${BLAS_LIB}")


SET(OPENBLAS_LIBRARY_DEBUG ${BLAS_LIB} CACHE FILEPATH "OpenBLAS library.")
SET(OPENBLAS_LIBRARY_RELEASE  ${BLAS_LIB} CACHE FILEPATH "OpenBLAS library.")

#-------------------------------------------------------------------------------
# Mark paths of OPENBLAS libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(OPENBLAS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(OPENBLAS_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set BLAS/LAPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------

# for OPENBLAS LAPACK_LIBRARY is the same as BLAS_LIBRARY. See also External_LAPACK for netlib and FindIntelMKL
# e.g. lib64/OPENSUSE_TUMBLEWEED_X86_64/libopenblas.a;-lpthread
IF(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")
  IF(DEBUG)
    SET(BLAS_LIBRARY "${OPENBLAS_LIBRARY_DEBUG}")
    SET(LAPACK_LIBRARY "${OPENBLAS_LIBRARY_DEBUG}")
  ELSE(DEBUG)
    SET(BLAS_LIBRARY "${OPENBLAS_LIBRARY_RELEASE}")
    SET(LAPACK_LIBRARY "${OPENBLAS_LIBRARY_RELEASE}")
  ENDIF(DEBUG)
ENDIF(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")
