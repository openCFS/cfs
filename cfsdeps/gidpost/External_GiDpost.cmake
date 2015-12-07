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

STRING(REPLACE ";" "," GID_FORTRAN_LIBS "${CFS_FORTRAN_LIBS}")
STRING(REPLACE ";" "," GID_HDF5_LIBRARY "${HDF5_LIBRARY}")

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
  -DZLIB_INCLUDE_DIR:FILEPATH=${CFS_BINARY_DIR}/include
  -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
  -DZLIB_LIBRARIES:FILEPATH=${ZLIB_LIBRARY}
  -DHDF5_DIR:FILEPATH=${CFS_BINARY_DIR}/cmake/hdf5
  -DHDF5_C_LIBRARY:PATH=${GID_HDF5_LIBRARY}
  -DHDF5_INCLUDE_DIR:FILEPATH=${CFS_BINARY_DIR}/include
  -DHDF5_LT_LIBRARY:FILEPATH=${HDF5_LT_LIBRARY}
  -DHDF5:BOOL=ON
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DBUILD_FORTRAN_EXAMPLES:BOOL=ON
  -DGID_FORTRAN_LIBS=${GID_FORTRAN_LIBS}
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
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://www.gidhome.com/pub/Tools/gidpost/old/${GIDPOST_ZIP}"
  "${CFS_FAU_MIRROR}/sources/gidpost/${GIDPOST_ZIP}"
  "${GIDPOST_URL}/${GIDPOST_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/gidpost/${GIDPOST_ZIP}")
SET(MD5_SUM ${GIDPOST_MD5})

SET(DLFN "${gidpost_prefix}/gidpost-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "gidpost" "${GIDPOST_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${gidpost_prefix}")

SET(ZIPFROMCACHE "${gidpost_prefix}/gidpost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${gidpost_prefix}/gidpost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The GiDpost external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(gidpost
    PREFIX "${gidpost_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(gidpost
    DEPENDS hdf5-static zlib
    PREFIX "${gidpost_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${GIDPOST_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR ,
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(gidpost cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${gidpost_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(gidpost cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

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
