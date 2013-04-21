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

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/zlib/zlib-patch.cmake.in")
SET(PFN "${zlib_prefix}/zlib-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The zlib external project
#-------------------------------------------------------------------------------
ExternalProject_Add(zlib
  PREFIX ${zlib_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/zlib
  SOURCE_DIR ${zlib_source}
  URL ${ZLIB_URL}/${ZLIB_GZ}
  URL_MD5 ${ZLIB_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  zlib
)

IF(MINGW)
  SET(ZLIB_LIB zlibstatic)
ELSE(MINGW)
  IF(UNIX)
    SET(ZLIB_LIB z)
  ELSE(UNIX)
    SET(ZLIB_LIB zlibstatic)
    IF(DEBUG)
      SET(ZLIB_LIB "${ZLIB_LIB}d")
    ENDIF()
  ENDIF(UNIX)
ENDIF(MINGW)

SET(ZLIB_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}${ZLIB_LIB}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "zlib library")
SET(ZLIB_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "zlib include directory")

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)
