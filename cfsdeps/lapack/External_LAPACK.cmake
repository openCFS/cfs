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
# The Lapack external project
#-------------------------------------------------------------------------------
ExternalProject_Add(lapack
  PREFIX "${lapack_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/lapack
  URL ${LAPACK_URL}/${LAPACK_GZ}
  URL_MD5 ${LAPACK_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lapack/lapack-patch.cmake.in")
SET(PFN "${lapack_prefix}/lapack-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(lapack custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${lapack_source}
)

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
