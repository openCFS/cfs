#-------------------------------------------------------------------------------
# Library of Iterative Solvers
#
# Project Homepage
# http://www.ssisc.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to lis sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(lis_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/lis")
set(lis_source  "${lis_prefix}/src/lis")
set(lis_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${lis_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
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

IF(OPENMP_FOUND)
  LIST(APPEND CMAKE_ARGS
    -DLIS_ENABLE_OPENMP:BOOL=ON
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lis/lis-patch.cmake.in")
SET(PFN "${lis_prefix}/lis-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.ssisc.org/lis/dl/${LIS_GZ}"
  "${LIS_URL}/${LIS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/lis/${LIS_GZ}")
SET(MD5_SUM ${LIS_MD5})

SET(DLFN "${lis_prefix}/lis-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "lis" "${LIS_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${lis_prefix}")

SET(ZIPFROMCACHE "${lis_prefix}/lis-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${lis_prefix}/lis-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of LIS libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(LIS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}lis${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "LIS library.")
MARK_AS_ADVANCED(LIS_LIBRARY)

#-------------------------------------------------------------------------------
# The lis external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(lis
    PREFIX "${lis_prefix}"
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
  ExternalProject_Add(lis
    PREFIX "${lis_prefix}"
    SOURCE_DIR "${lis_source}"
    URL ${LOCAL_FILE}
    URL_MD5 ${LIS_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DLIS_BUILD_TEST:BOOL=ON
      -DLIS_ENABLE_FORTRAN:BOOL=ON
      -DLIS_ENABLE_SAAMG:BOOL=ON
    BUILD_BYPRODUCTS ${LIS_LIBRARY}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(lis cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${lis_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(lis cfsdeps_zipToCache
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
  lis
)

SET(LIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(LIS_INCLUDE_DIR)

