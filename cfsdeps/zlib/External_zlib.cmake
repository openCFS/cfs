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
)

#-------------------------------------------------------------------------------
# The zlib-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(zlib-shared
  PREFIX ${zlib_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/zlib
  SOURCE_DIR ${zlib_source}
  URL ${ZLIB_URL}/${ZLIB_GZ}
  URL_MD5 ${ZLIB_MD5}
  CMAKE_ARGS
     ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=ON
    )

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
ExternalProject_Add_Step(zlib-shared custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${zlib_source}
)
  
#-------------------------------------------------------------------------------
# The zlib-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(zlib-static
  DEPENDS zlib-shared
  PREFIX ${zlib_prefix}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${zlib_source}  
  CMAKE_ARGS
     ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    )

SET(ZLIB_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libz.a CACHE FILEPATH "zlib library")
SET(ZLIB_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "zlib include directory")

MARK_AS_ADVANCED(ZLIB_LIBRARY)
MARK_AS_ADVANCED(ZLIB_INCLUDE_DIR)
