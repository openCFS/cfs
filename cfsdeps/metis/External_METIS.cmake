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

set(metis_c_flags ${CFS_SUPPRESSIONS})
if(CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(metis_c_flags "-Wno-implicit-function-declaration ") # Apple clang version 12.0.0 
endif()

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${metis_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
   -DCMAKE_C_FLAGS:STRING=${metis_c_flags}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX})

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT})
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE})
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/metis/metis-patch.cmake.in")
SET(PFN "${metis_prefix}/metis-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "ftp://ftp1.rrzn.uni-hannover.de/pub/mirror/bsd/FreeBSD/ports/distfiles/${METIS_GZ}"
  "http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD/${METIS_GZ}"
  "${CFS_DS_SOURCES_DIR}/metis/${METIS_GZ}")
  
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/metis/${METIS_GZ}")
SET(MD5_SUM ${METIS_MD5})

SET(DLFN "${metis_prefix}/metis-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/metis/license/" DESTINATION "${CFS_BINARY_DIR}/license/metis" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "metis" "${METIS_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${metis_prefix}")

SET(ZIPFROMCACHE "${metis_prefix}/metis-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${metis_prefix}/metis-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of METIS libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
IF(NOT USE_ILUPACK_PARALLEL)
SET(METIS_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}metis${CMAKE_STATIC_LIBRARY_SUFFIX};"
    CACHE FILEPATH "METIS library.")
ENDIF()
MARK_AS_ADVANCED(METIS_LIBRARY)

#-------------------------------------------------------------------------------
# The metis external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(metis
    PREFIX "${metis_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${METIS_LIBRARY}
  )
  # This is still required for the build to work properly since we need to replace the metis related files from ilupack with old metis
  IF(USE_ILUPACK_PARALLEL)
    add_dependencies(metis ilupack)
  ENDIF()

ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(metis
    PREFIX "${metis_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${METIS_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${METIS_LIBRARY}    
  )

  # The ilupack has a metis.h from a newer version of metis which is not compatible with CFS, so building it in this way 
  # replaces the ilupack metis.h with the older metis library's metis.h.
  IF(USE_ILUPACK_PARALLEL)
    add_dependencies(metis ilupack)
  ENDIF()

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(metis cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${metis_prefix}
  )
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(metis cfsdeps_zipToCache
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
SET(CFSDEPS ${CFSDEPS} metis)

SET(METIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
