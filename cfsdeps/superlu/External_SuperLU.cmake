#-------------------------------------------------------------------------------
# SuperLU is a general purpose library for the direct solution of large,
# sparse, nonsymmetric systems of linear equations on high performance
# machines.
#
# Project Homepage
# http://crd-legacy.lbl.gov/~xiaoye/SuperLU/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to SuperLU sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(superlu_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/superlu")
set(superlu_source  "${superlu_prefix}/src/superlu")
set(superlu_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(SUPERLU_BLAS_LIBRARY "${BLAS_LIBRARY}")

IF(USE_BLAS_LAPACK STREQUAL "NETLIB")
  LIST(APPEND SUPERLU_BLAS_LIBRARY "${CFS_FORTRAN_LIBS}")
ELSEIF(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  LIST(APPEND SUPERLU_BLAS_LIBRARY "-lm")
ENDIF()

STRING(REPLACE ";" "^" SUPERLU_BLAS_LIBRARY "${SUPERLU_BLAS_LIBRARY}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${superlu_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DSUPERLU_USE_BUNDLED_BLAS=OFF
  -DBLAS_LIBRARIES:FILEPATH=${SUPERLU_BLAS_LIBRARY}
  -DSUPERLU_ENABLE_TESTING=ON
  )

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_CXX_FLAGS:FILEPATH=${CFLAGS}
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/superlu/superlu-patch.cmake.in")
SET(PFN "${superlu_prefix}/superlu-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://crd-legacy.lbl.gov/~xiaoye/SuperLU/${SUPERLU_GZ}"
  "${CFS_DS_SOURCES_DIR}/superlu/${SUPERLU_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/superlu/${SUPERLU_GZ}")
SET(MD5_SUM ${SUPERLU_MD5})

SET(DLFN "${superlu_prefix}/superlu-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

# Since CMake may have problems with symlinks inside archives on Windows, we
# have to unpack the archive using the GNU tar and gunzip utilities.
SET(DLFN_UNTAR "${superlu_prefix}/superlu-download-untar.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cfsdeps/superlu/superlu-download-untar.cmake.in"
  "${DLFN_UNTAR}"
  @ONLY
)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/superlu/license/" DESTINATION "${CFS_BINARY_DIR}/license/superlu" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "superlu" "${SUPERLU_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${superlu_prefix}")

SET(ZIPFROMCACHE "${superlu_prefix}/superlu-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${superlu_prefix}/superlu-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of SuperLU libraries.
#-----------------------------------------------------------------------------
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(SUPERLU_LIBRARY "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}")
ELSE()
  SET(SUPERLU_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE FILEPATH "SuperLU library.")
ENDIF()

MARK_AS_ADVANCED(SUPERLU_LIBRARY)

#-------------------------------------------------------------------------------
# The superlu external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(superlu
    PREFIX "${superlu_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${SUPERLU_LIBRARY}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(superlu
    PREFIX "${superlu_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${SUPERLU_MD5}
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${DLFN_UNTAR}"
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR "^"
    CMAKE_ARGS
      ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${SUPERLU_LIBRARY}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(superlu cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${superlu_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(superlu cfsdeps_zipToCache
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
  superlu
)

SET(SUPERLU_INCLUDE_DIR "${CFS_BINARY_DIR}/include/superlu")
MARK_AS_ADVANCED(SUPERLU_INCLUDE_DIR)

