#-------------------------------------------------------------------------------
# METIS - Serial Graph Partitioning and Fill-reducing Matrix Ordering
# Needed by Ilupack
#
# Project Homepage
# http://glaros.dtc.umn.edu/gkhome/views/metis
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to metis sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(metis_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/metis")
set(metis_source  "${metis_prefix}/src/metis")
set(metis_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${metis_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/metis/metis-patch.cmake.in")
SET(PFN "${metis_prefix}/metis-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The metis external project
#-------------------------------------------------------------------------------
ExternalProject_Add(metis
  PREFIX "${metis_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/metis
  URL ${METIS_URL}/${METIS_GZ}
  URL_MD5 ${METIS_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  metis
)

SET(METIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of METIS libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(METIS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "METIS library.")

MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
MARK_AS_ADVANCED(METIS_LIBRARY)
