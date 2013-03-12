#-------------------------------------------------------------------------------
# Set paths to bzip2 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(bzip2_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/bzip2")
set(bzip2_source  "${bzip2_prefix}/src/bzip2")
set(bzip2_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${bzip2_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  )

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_C_FLAGS:FILEPATH=${CFLAGS}
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
# The bzip2 external project
#-------------------------------------------------------------------------------
ExternalProject_Add(bzip2-static
  PREFIX "${bzip2_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/bzip2
  SOURCE_DIR ${bzip2_source}
  URL ${BZIP2_URL}/${BZIP2_GZ}
  URL_MD5 ${BZIP2_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=OFF
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/bzip2/bzip2-patch.cmake.in")
SET(PFN "${bzip2_prefix}/bzip2-patch.cmake")
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
ExternalProject_Add_Step(bzip2-static custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${bzip2_source}
)

#-------------------------------------------------------------------------------
# The bzip2 external project
#-------------------------------------------------------------------------------
ExternalProject_Add(bzip2-shared
  DEPENDS bzip2-static
  PREFIX "${bzip2_prefix}"
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${bzip2_source}
  CMAKE_ARGS
    ${CMAKE_ARGS}
    -DBUILD_SHARED_LIBS:BOOL=ON
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  bzip2-static
  bzip2-shared
)

SET(BZIP2_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of BZIP2 libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(BZIP2_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "BZIP2 library.")

MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR)
MARK_AS_ADVANCED(BZIP2_LIBRARY)
