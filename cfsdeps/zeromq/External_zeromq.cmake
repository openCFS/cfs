#-------------------------------------------------------------------------------
# ZeroMQ
# Distributed Messaging Library
#
# Project Homepage
# http://zeromq.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to ZeroMQ sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(zeromq_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zeromq")
set(zeromq_source  "${zeromq_prefix}/src/zeromq")
set(zeromq_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${zeromq_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
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
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/zeromq/zeromq-patch.cmake.in")
#SET(PFN "${zeromq_prefix}/zeromq-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://download.zeromq.org/${ZEROMQ_GZ}"
  "${ZEROMQ_URL}/${ZEROMQ_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/zeromq/${ZEROMQ_GZ}")
SET(MD5_SUM ${ZEROMQ_MD5})

SET(DLFN "${zeromq_prefix}/zeromq-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The zeromq external project
#-------------------------------------------------------------------------------
ExternalProject_Add(zeromq
  PREFIX ${zeromq_prefix}
  URL ${LOCAL_FILE}
  URL_MD5 ${ZEROMQ_MD5}
  SOURCE_DIR ${zeromq_source}
  CONFIGURE_COMMAND sh configure --without-libsodium --prefix=${zeromq_prefix}
  CMAKE_ARGS ${CMAKE_ARGS}
  BUILD_IN_SOURCE 1
)

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(zeromq cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${zeromq_prefix}
)

#-------------------------------------------------------------------------------
# Add custom step for downloading and installing C++ binding (header only)
#-------------------------------------------------------------------------------
SET(cppzmq_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cppzmq")
ExternalProject_Add(cppzmq
  DEPENDS zeromq
  PREFIX "${cppzmq_prefix}"
  #  PREFIX "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cppzmq"
  GIT_REPOSITORY "https://github.com/zeromq/cppzmq.git"
  GIT_TAG "master"
  PATCH_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${cppzmq_prefix}/src/cppzmq/zmq.hpp" "${CFS_BINARY_DIR}/include"
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
LIST(APPEND CFSDEPS zeromq cppzmq)

IF(MINGW)
  SET(ZEROMQ_LIB zmqstatic)
  SET(ZEROMQ_SHARED_LIB zmq)
ELSE(MINGW)
  IF(UNIX)
    SET(ZEROMQ_LIB zmq)
    SET(ZEROMQ_SHARED_LIB zmq)
  ELSE(UNIX)
    SET(ZEROMQ_LIB zmqstatic)
    SET(ZEROMQ_SHARED_LIB zmq)
    IF(DEBUG)
      SET(ZEROMQ_LIB "${ZEROMQ_LIB}d")
      SET(ZEROMQ_SHARED_LIB "${ZEROMQ_SHARED_LIB}d")
    ENDIF()
  ENDIF(UNIX)
ENDIF(MINGW)

SET(LD ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR})
SET(ZEROMQ_LIBRARY
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${ZEROMQ_LIB}${CMAKE_STATIC_LIBRARY_SUFFIX}
  CACHE FILEPATH "ZeroMQ library")
SET(ZEROMQ_SHARED_LIBRARY
  ${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${ZEROMQ_SHARED_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
IF(MINGW)
  SET(ZEROMQ_SHARED_LIBRARY "${ZEROMQ_SHARED_LIBRARY}.a")
ENDIF(MINGW)
SET(ZEROMQ_SHARED_LIBRARY ${ZEROMQ_SHARED_LIBRARY} CACHE FILEPATH "ZeroMQ shared library")
SET(ZEROMQ_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "ZeroMQ include directory")

MARK_AS_ADVANCED(ZEROMQ_LIBRARY)
MARK_AS_ADVANCED(ZEROMQ_SHARED_LIBRARY)
MARK_AS_ADVANCED(ZEROMQ_INCLUDE_DIR)

