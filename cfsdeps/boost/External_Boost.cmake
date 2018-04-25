#-------------------------------------------------------------------------------
# Boost provides free peer-reviewed portable C++ source libraries.
# Needed by Botan and CGAL
#
# Project Homepage
# http://www.boost.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to Boost sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(BOOST_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/boost")
set(BOOST_source  "${BOOST_prefix}/src/boost")
set(BOOST_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${BOOST_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
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
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/boost/boost-patch.cmake.in")
#SET(PFN "${ARPACK_prefix}/arpack-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/boost/${BOOST_GZ}")

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "boost" "${BOOST_VER}") 

# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${BOOST_prefix}")

SET(ZIPFROMCACHE "${BOOST_prefix}/boost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ARPACK_prefix}/boost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of BOOST libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

#-------------------------------------------------------------------------------
# The BOOST external project
#-------------------------------------------------------------------------------

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(boost
    PREFIX "${BOOST_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(boost
    PREFIX "${BOOST_prefix}"
    URL "https://dl.bintray.com/boostorg/release/${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}.0/source/${BOOST_GZ}"
    URL_MD5 "${BOOST_MD5}"
    CONFIGURE_COMMAND ""
    BINARY_DIR ${BOOST_source}
    BUILD_COMMAND ./bootstrap.sh --prefix=${BOOST_install}
    INSTALL_COMMAND ./b2 install
  )

  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(boost cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS ${CFSDEPS} boost)

