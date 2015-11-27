#-------------------------------------------------------------------------------
# SCPIP - - an efficient software tool for the solution of
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
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${scpip_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:FILEPATH=-w
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
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/scpip/scpip-patch.cmake.in")
SET(PFN "${scpip_prefix}/scpip-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

# SCPIP is closed source, the code and binary may only be used after 
# signing an agreement with Ch. Zillober
# CFS expects scpip.tar.bz2 in CFS_DEPS_CACHE_DIR/source/scpip

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/scpip/${SCPIP_BZ2}")
SET(MD5_SUM ${SCPIP_MD5})

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "scpip" "${SCPIP_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${scpip_prefix}")

SET(ZIPFROMCACHE "${scpip_prefix}/scpip-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${scpip_prefix}/scpip-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The scpip external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
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
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(scpip
    PREFIX "${scpip_prefix}"
    SOURCE_DIR "${scpip_source}"
    URL ${CFS_DEPS_CACHE_DIR}/sources/scpip/${SCPIP_BZ2}
    URL_MD5 ${SCPIP_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )
    
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(scpip cfsdeps_zipToCache
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
  scpip
)

SET(SCPIP_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of SCPIP libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(SCPIP_LIBRARY "${LD}/libscpip.a" CACHE FILEPATH "SCPIP library.")

MARK_AS_ADVANCED(SCPIP_LIBRARY)
MARK_AS_ADVANCED(SCPIP_INCLUDE_DIR)

