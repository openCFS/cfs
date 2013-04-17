#-------------------------------------------------------------------------------
# GiDpost: a C / C++ / Fortran library to create postprocess files for GiD
#
# Project Homepage
# http://gid.cimne.upc.es/gid-plus/tools/gidpost
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to gidpost sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(gidpost_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/gidpost")
set(gidpost_source  "${gidpost_prefix}/src/gidpost")
set(gidpost_install  "${CFS_BINARY_DIR}")

IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
  SET(GIDPOST_C_FLAGS "-DINTEL_COMPILER ${CFSDEPS_C_FLAGS}")
ELSE()
  SET(GIDPOST_C_FLAGS "-DgFortran ${CFSDEPS_C_FLAGS}")
ENDIF()

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${gidpost_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${GIDPOST_C_FLAGS}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DZLIB_INCLUDE_DIRS:FILEPATH=${CFS_BINARY_DIR}/include
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  -DZLIB_LIBRARIES:FILEPATH=${ZLIB_LIBRARY}
  -DHDF5_DIR:FILEPATH=${CFS_BINARY_DIR}
  -DHDF5_INCLUDE_DIRS:FILEPATH=${CFS_BINARY_DIR}/include
  -DHDF5:BOOL=ON
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DBUILD_FORTRAN_EXAMPLES:BOOL=ON
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gidpost/gidpost-patch.cmake.in")
SET(PFN "${gidpost_prefix}/gidpost-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# The GiDpost external project
#-------------------------------------------------------------------------------
ExternalProject_Add(gidpost
  DEPENDS hdf5-static zlib
  PREFIX "${gidpost_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/gidpost
  URL ${GIDPOST_URL}/${GIDPOST_ZIP}
  URL_MD5 ${GIDPOST_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CMAKE_ARGS
     ${CMAKE_ARGS}
  )

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  gidpost
)

#-------------------------------------------------------------------------------
# These variables are used to find Gidpost by other projects
#-------------------------------------------------------------------------------
SET(GIDPOST_INCLUDE_DIR
  "${CFS_BINARY_DIR}/include"
  CACHE
  FILEPATH
  "GiDpost include dir"
  FORCE)

#-------------------------------------------------------------------------------
# Determine paths of GIDPOST libraries.
#-------------------------------------------------------------------------------
SET(GIDPOST_LIBRARY_DEBUG
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)
SET(GIDPOST_LIBRARY_RELEASE
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)

#-------------------------------------------------------------------------------
# Mark paths of GIDPOST libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(GIDPOST_INCLUDE_DIR)
MARK_AS_ADVANCED(GIDPOST_LIBRARY_DEBUG)
MARK_AS_ADVANCED(GIDPOST_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set GIDPOST_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_RELEASE}")
ENDIF(DEBUG)
