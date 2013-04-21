#-------------------------------------------------------------------------------
# Library of Iterative Solvers
#
# Project Homepage
# http://www.ssisc.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to lis sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(lis_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/lis")
set(lis_source  "${lis_prefix}/src/lis")
set(lis_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${lis_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
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

IF(OPENMP_FOUND)
  LIST(APPEND CMAKE_ARGS
    -DLIS_ENABLE_OPENMP:BOOL=ON
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lis/lis-patch.cmake.in")
SET(PFN "${lis_prefix}/lis-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The lis external project
#-------------------------------------------------------------------------------
ExternalProject_Add(lis
  PREFIX "${lis_prefix}"
  SOURCE_DIR "${lis_source}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/lis
  URL ${LIS_URL}/${LIS_GZ}
  URL_MD5 ${LIS_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
#    -DLIS_BUILD_TEST:BOOL=ON
    -DLIS_ENABLE_FORTRAN:BOOL=ON
    -DLIS_ENABLE_SAAMG:BOOL=ON
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  lis
)

SET(LIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of LIS libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(LIS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}lis${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "LIS library.")

MARK_AS_ADVANCED(LIS_LIBRARY)
MARK_AS_ADVANCED(LIS_INCLUDE_DIR)

