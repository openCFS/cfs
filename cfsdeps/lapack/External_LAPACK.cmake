#-------------------------------------------------------------------------------
# Set paths to lapack sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(lapack_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/lapack")
set(lapack_source  "${lapack_prefix}/src/lapack")
set(lapack_install  "${CFS_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${lapack_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DBUILD_TESTING:BOOL=OFF
#  -DBUILD_SINGLE:BOOL=OFF
#  -DBUILD_COMPLEX:BOOL=OFF
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

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lapack/lapack-patch.cmake.in")
SET(PFN "${lapack_prefix}/lapack-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp.mirrorservice.org/sites/distfiles.macports.org/atlas/${LAPACK_GZ}"
  "http://www.netlib.org/lapack/${LAPACK_GZ}"
  "${LAPACK_URL}/${LAPACK_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/lapack/${LAPACK_GZ}")
SET(MD5_SUM ${LAPACK_MD5})

SET(DLFN "${lapack_prefix}/lapack-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "lapack_${LAPACK_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}_${CMAKE_BUILD_TYPE}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "lapack_${LAPACK_VER}_${CFS_ARCH_STR}_${FC_ID}_${CMAKE_BUILD_TYPE}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${lapack_prefix}")

SET(ZIPFROMCACHE "${lapack_prefix}/lapack-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${lapack_prefix}/lapack-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The lapack external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(lapack
    PREFIX "${lapack_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(lapack
    PREFIX "${lapack_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/lapack
    URL ${LOCAL_FILE}
    URL_MD5 ${LAPACK_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(lapack cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${lapack_prefix}
  )
  
  IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(lapack cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  lapack
)

#-------------------------------------------------------------------------------
# Determine paths of LAPACK libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(NETLIB_BLAS_LIBRARY_DEBUG
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_BLAS_LIBRARY_RELEASE
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_LAPACK_LIBRARY_DEBUG
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib LAPACK library.")
  SET(NETLIB_LAPACK_LIBRARY_RELEASE
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib LAPACK library.")
ELSE(WIN32)
  SET(BLAS_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libblas.a")
  SET(LAPACK_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/liblapack.a")

  IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
    LIST(APPEND BLAS_LIB ${GFORTRAN_LIBRARY})
  ENDIF()

  SET(NETLIB_BLAS_LIBRARY_DEBUG
    ${BLAS_LIB}
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_BLAS_LIBRARY_RELEASE
    ${BLAS_LIB}
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_LAPACK_LIBRARY_DEBUG
    ${LAPACK_LIB}
    CACHE FILEPATH "Netlib LAPACK library.")
  SET(NETLIB_LAPACK_LIBRARY_RELEASE
    ${LAPACK_LIB}
    CACHE FILEPATH "Netlib LAPACK library.")
ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Mark paths of LAPACK libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(NETLIB_BLAS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(NETLIB_BLAS_LIBRARY_RELEASE)
MARK_AS_ADVANCED(NETLIB_LAPACK_LIBRARY_DEBUG)
MARK_AS_ADVANCED(NETLIB_LAPACK_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set LAPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
  IF(DEBUG)
    SET(BLAS_LIBRARY "${NETLIB_BLAS_LIBRARY_DEBUG}")
    SET(LAPACK_LIBRARY "${NETLIB_LAPACK_LIBRARY_DEBUG}")
  ELSE(DEBUG)
    SET(BLAS_LIBRARY "${NETLIB_BLAS_LIBRARY_RELEASE}")
    SET(LAPACK_LIBRARY "${NETLIB_LAPACK_LIBRARY_RELEASE}")
  ENDIF(DEBUG)

  IF(MINGW)
    LIST(APPEND BLAS_LIBRARY ${GFORTRAN_LIBRARY})
  ENDIF(MINGW)
ENDIF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
