#-------------------------------------------------------------------------------
# STAR-CCM+ I/O library
# Needed for reading .ccm files.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to CCMIO sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ccmio_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ccmio")
set(ccmio_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(ccmio_source  "${ccmio_prefix}/src/ccmio")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${ccmio_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFLAGS}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  # Explicitly set build architectures and  system SDK root dir to match
  # those given in CMake.
  
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
  )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ccmio/ccmio-patch.cmake.in")
SET(PFN "${ccmio_prefix}/ccmio-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://portal.nersc.gov/svn/visit/trunk/third_party/${CCMIO_GZ}"
# The following files have mostly the same contents as the first one, but have
# a slightly different size and MD5 sum.
#  "http://priede.bf.lu.lv/ftp/pub/GIS/bibliotekas/dazadas/${CCMIO_GZ}"
#  "ftp://www.bf.lu.lv/pub/TIS/bibliotekas/dazadas/${CCMIO_GZ}"
#  "ftp://196.203.130.15/pub/logiciels/visit/${CCMIO_GZ}"
  "${CCMIO_URL}/${CCMIO_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ccmio/${CCMIO_GZ}")
SET(MD5_SUM ${CCMIO_MD5})

SET(DLFN "${ccmio_prefix}/ccmio-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "ccmio" "${CCMIO_VER}")  
  
# This should be either PREFIX_DIR/src (install manifest is used for zipping)
# or PREFIX_DIR/install (install directory will be zipped)
SET(TMP_DIR "${ccmio_prefix}/src")

SET(ZIPFROMCACHE "${ccmio_prefix}/ccmio-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ccmio_prefix}/ccmio-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The CCMIO-static external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ccmio-static
    PREFIX "${ccmio_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  ExternalProject_Add(ccmio-static
    DEPENDS cgns
    PREFIX ${ccmio_prefix}
    SOURCE_DIR ${ccmio_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${CCMIO_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR ,
    CMAKE_ARGS
       ${CMAKE_ARGS}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ccmio-static cfsdeps_download
     COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
     DEPENDERS download
     DEPENDS "${DLFN}"
     WORKING_DIRECTORY ${ccmio_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ccmio-static cfsdeps_zipToCache
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
  ccmio-static
)

#-------------------------------------------------------------------------------
# Determine paths of CCMIO libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CCMIO_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}ccmio${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "CCMIO library.")

SET(CCMIO_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "CCMIO include directory")

MARK_AS_ADVANCED(CCMIO_INCLUDE_DIR)
MARK_AS_ADVANCED(CCMIO_LIBRARY)
