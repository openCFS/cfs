#-------------------------------------------------------------------------------
# 
# OpenBLAS is a free and efficient and parallel BLAS/ LAPACK implementation.
# It is a successor of the stalled GOTO BLAS
# http://www.openblas.net/
# 
#-------------------------------------------------------------------------------
set(openblas_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/openblas")
set(openblas_source  "${openblas_prefix}/src/openblas")
set(openblas_install  "${CFS_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${openblas_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DBUILD_TESTING:BOOL=OFF
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()


#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://github.com/xianyi/OpenBLAS/tarball/v0.2.10/${OPENBLAS_GZ}"
  "${OPENBLAS_URL}/${OPENBLAS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/openblas/${OPENBLAS_GZ}")
SET(MD5_SUM ${OPENBLAS_MD5})

SET(DLFN "${openblas_prefix}/openblas-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )
  
#-------------------------------------------------------------------------------
# Set names of patch file and template file.
# Basically add the option CFS_TARGET_LIB to Makefile.install
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/openblas/openblas-patch.cmake.in")
SET(PFN "${openblas_prefix}/openblas-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 
  
#-------------------------------------------------------------------------------
# The OpenBLAS external project, CFS_LIB_TARGET from the patched Makefile.install
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
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CONFIGURE_COMMAND ""
  BUILD_COMMAND  ${CMAKE_MAKE_PROGRAM} libs netlib
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install NO_SHARED=1 "PREFIX=${openblas_install}" "CFS_LIB_TARGET=${LIBRARY_OUTPUT_PATH}"
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

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  openblas
)

#-------------------------------------------------------------------------------
# Determine paths of OpenBLAS libraries.
# There is no separation between BLAS and LAPACK, both are combined. While BLAS is optimized
# only parts of LAPACK are optimized. OpenBLAS can be build without LAPACK in request.
#-------------------------------------------------------------------------------
SET(BLAS_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libopenblas.a;-lpthread")
SET(LAPACK_LIB "${BLAS_LIB}")

#seems not to be necessary?!
#IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
#  LIST(APPEND BLAS_LIB ${GFORTRAN_LIBRARY})
#ENDIF()

SET(OPENBLAS_LIBRARY_DEBUG
  ${BLAS_LIB}
  CACHE FILEPATH "OpenBLAS library.")
SET(OPENBLAS_LIBRARY_RELEASE
  ${BLAS_LIB}
  CACHE FILEPATH "OpenBLAS library.")

#-------------------------------------------------------------------------------
# Mark paths of OPENBLAS libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(OPENBLAS_LIBRARY_DEBUG)
MARK_AS_ADVANCED(OPENBLAS_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set BLAS/LAPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")
  IF(DEBUG)
    SET(BLAS_LIBRARY "${OPENBLAS_LIBRARY_DEBUG}")
    SET(LAPACK_LIBRARY "${OPENBLAS_LIBRARY_DEBUG}")
  ELSE(DEBUG)
    SET(BLAS_LIBRARY "${OPENBLAS_LIBRARY_RELEASE}")
    SET(LAPACK_LIBRARY "${OPENBLAS_LIBRARY_RELEASE}")
  ENDIF(DEBUG)

ENDIF(CFS_BLAS_LAPACK STREQUAL "OPENBLAS")
