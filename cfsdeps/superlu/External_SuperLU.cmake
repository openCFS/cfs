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

IF(CFS_BLAS_LAPACK STREQUAL "NETLIB")
  LIST(APPEND SUPERLU_BLAS_LIBRARY "${CFS_FORTRAN_LIBS}")
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
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
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
  "${SUPERLU_URL}/${SUPERLU_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/superlu/${SUPERLU_GZ}")
SET(MD5_SUM ${SUPERLU_MD5})

SET(DLFN "${superlu_prefix}/superlu-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The superlu external project
#-------------------------------------------------------------------------------
ExternalProject_Add(superlu
  PREFIX "${superlu_prefix}"
  URL ${LOCAL_FILE}
  URL_MD5 ${SUPERLU_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  LIST_SEPARATOR "^"
  CMAKE_ARGS
    ${CMAKE_ARGS}
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

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  superlu
)

SET(SUPERLU_INCLUDE_DIR "${CFS_BINARY_DIR}/include/superlu")

#-----------------------------------------------------------------------------
# Determine paths of SuperLU libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")


IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(SUPERLU_LIBRARY "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}")
ELSEIF(MINGW)
  SET(SUPERLU_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE FILEPATH "SuperLU library.")
ELSE()
  SET(SUPERLU_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}superlu${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE FILEPATH "SuperLU library.")
ENDIF()

MARK_AS_ADVANCED(SUPERLU_LIBRARY)
MARK_AS_ADVANCED(SUPERLU_INCLUDE_DIR)

