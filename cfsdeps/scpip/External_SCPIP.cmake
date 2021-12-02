#-------------------------------------------------------------------------------
# SCPIP - an efficient software tool for the solution of
# structural optimization problems. 
# Implements the Method of Moving Asymptotes (MMA) from Krister Svanberg
#
# Project Homepage
# http://www.mathematik.uni-wuerzburg.de/~zillober
# http://www.mathematik.uni-wuerzburg.de/~zillober/pubs/manual30.pdf
#-------------------------------------------------------------------------------
#
# SCPIP is not open source! Users need to obtain a permission to use the
# code for academic purupose!  
# 
# The code 
#-------------------------------------------------------------------------------
# Set paths to scpip sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(scpip_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/scpip")
set(scpip_source  "${scpip_prefix}/src/scpip")
set(scpip_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
set(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${scpip_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:FILEPATH=-w
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX} )

if(CFS_DISTRO STREQUAL "MACOSX")
  set(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT})
endif()

if(CMAKE_TOOLCHAIN_FILE)
  list(APPEND CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE})
endif()

# the patch file decrypts and patches the content
set_from_env(CFS_KEY_SCPIP)
set(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/scpip/scpip-patch.cmake.in")
set(PFN "${scpip_prefix}/scpip-patch.cmake")
configure_file("${PFN_TEMPL}" "${PFN}" @ONLY) 

# SCPIP is closed source, the code and binary may only be used after 
# signing an agreement with Ch. Zillober
# CFS expects scpip.tar.bz2 in CFS_DEPS_CACHE_DIR/source/scpip

set_from_env(CFS_DOWNLOAD_SCPIP)
set(MIRRORS "${CFS_DS_SOURCES_DIR}/${CFS_DOWNLOAD_SCPIP}/${SCPIP_ZIP}")
set(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/scpip/${SCPIP_ZIP}")
set(MD5_SUM ${SCPIP_MD5})

set(DLFN "${scpip_prefix}/scpip-download.cmake")
# uses the variables MIRRORS, LOCAL_FILE and MD5_SUM
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "scpip" "${SCPIP_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
set(TMP_DIR "${scpip_prefix}")

set(ZIPFROMCACHE "${scpip_prefix}/scpip-zipFromCache.cmake")
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${scpip_prefix}/scpip-zipToCache.cmake")
configure_file("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of SCPIP libraries.
#-----------------------------------------------------------------------------
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
if(UNIX) # also mac
  set(SCPIP_LIBRARY "${LD}/libscpip.a" CACHE FILEPATH "SCPIP library.")
else()
  set(SCPIP_LIBRARY "${LD}/scpip.lib" CACHE FILEPATH "SCPIP library.")
endif()
mark_as_advanced(SCPIP_LIBRARY)

#-------------------------------------------------------------------------------
# The scpip external project
#-------------------------------------------------------------------------------
if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(scpip
    PREFIX "${scpip_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${SCPIP_LIBRARY} )
else()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------

  # variables can be set in .cfs_platform_defaults.cmake, with cmake -D or as environment variables.
  if(NOT EXISTS LOCAL_FILE AND NOT DEFINED CFS_DOWNLOAD_SCPIP)
    message(WARNING "no local scpip sources exists and CFS_DOWNLOAD_SCPIP not set.")
  endif()

  if(NOT DEFINED CFS_KEY_SCPIP)
    message(WARNING "CFS_KEY_SCPIP not set.")
  endif()  
  
  ExternalProject_Add(scpip
    PREFIX "${scpip_prefix}"
    SOURCE_DIR "${scpip_source}"
    URL "${LOCAL_FILE}"
    URL_MD5 "${SCPIP_MD5}"
    # the patch needs to upack the encrypted file and apply patches to enable compiling - both in the patch command
    # the patch does the unpacking with password
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${SCPIP_LIBRARY} )
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(scpip cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${scpip_prefix} )

  if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(scpip cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR} )
  endif()
endif()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
set(CFSDEPS ${CFSDEPS} scpip )

set(SCPIP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
mark_as_advanced(SCPIP_INCLUDE_DIR)
