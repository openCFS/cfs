#-------------------------------------------------------------------------------
# Xerces-C++ is a validating XML parser written in a portable subset of C++.
#
# Project Homepage
# http://xerces.apache.org/xerces-c/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to Xerces-C sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(xerces_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/xerces")
set(xerces_source  "${xerces_prefix}/src/xerces")
set(xerces_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${xerces_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
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

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-patch.cmake.in")
SET(PFN "${xerces_prefix}/xerces-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since lse17 may not be
# accessible from behind firewalls.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp.de.cw.net/pub/FreeBSD/ports/distfiles/xerces-c-3.1.1.tar.gz"
  "http://xml.apache.org/dist/xerces-c/3/sources/xerces-c-3.1.1.tar.gz"
)

#-------------------------------------------------------------------------------
# Try to download sources into CFSDEPS cache directory.
#-------------------------------------------------------------------------------
DOWNLOAD_CFSDEPS(
  "${CFS_DEPS_CACHE_DIR}/sources/xerces/${XERCES_GZ}"
  ${XERCES_MD5}
  "${MIRRORS}"
)

#-------------------------------------------------------------------------------
# The xerces external project
#-------------------------------------------------------------------------------
ExternalProject_Add(xerces
  PREFIX "${xerces_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/xerces
  URL ${XERCES_URL}/${XERCES_GZ}
  URL_MD5 ${XERCES_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  xerces
)

SET(XERCES_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of XERCES libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")


IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(XERCES_LIBRARY "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX}")
  LIST(APPEND XERCES_LIBRARY "-framework CoreFoundation")
  LIST(APPEND XERCES_LIBRARY "-framework CoreServices")
ELSEIF(MINGW)
  SET(XERCES_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE FILEPATH "XERCES library.")
ELSE()
  SET(XERCES_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX};-lpthread"
    CACHE FILEPATH "XERCES library.")
ENDIF()

MARK_AS_ADVANCED(XERCES_LIBRARY)
MARK_AS_ADVANCED(XERCES_INCLUDE_DIR)

