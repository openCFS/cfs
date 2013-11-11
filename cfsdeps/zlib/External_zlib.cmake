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

IF(MINGW)
  LIST(APPEND CMAKE_ARGS
    -DCFS_ARCH:STRING=${CFS_ARCH}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/zlib/zlib-patch.cmake.in")
SET(PFN "${zlib_prefix}/zlib-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp.pl.pgpi.org/vol/rzm1/GraphicsMagick/delegates/zlib-1.2.7.tar.gz"
  "ftp://ftp.uwsg.indiana.edu/linux/gentoo/distfiles/zlib-1.2.7.tar.gz"
  "${ZLIB_URL}/${ZLIB_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/zlib/${ZLIB_GZ}")
SET(MD5_SUM ${ZLIB_MD5})

SET(DLFN "${zlib_prefix}/zlib-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The zlib external project
#-------------------------------------------------------------------------------
ExternalProject_Add(zlib
  PREFIX ${zlib_prefix}
  SOURCE_DIR ${zlib_source}
  URL ${LOCAL_FILE}
  URL_MD5 ${ZLIB_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(zlib cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${zlib_prefix}
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
