#-------------------------------------------------------------------------------
# FEAST Eigenvalue Solver
#
# Project Homepage
# http://www.feast-solver.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to FEAST sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(feast_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/feast")
set(feast_source  "${feast_prefix}/src/feast")

#-------------------------------------------------------------------------------
# Configure FEAST by copying the config file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/feast/make.inc.in")
SET(CONF "${feast_prefix}/make.inc"	)
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY)

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.ecs.umass.edu/~polizzi/feast/m3-0/feast_3.0.tgz"
  "${FEAST_URL}/${FEAST_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/feast/${FEAST_GZ}")
SET(MD5_SUM ${FEAST_MD5})

SET(DLFN "${feast_prefix}/feast-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "feast" "${FEAST_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${feast_prefix}")

SET(ZIPFROMCACHE "${feast_prefix}/feast-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${feast_prefix}/feast-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of FEAST libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(FEAST_LIBS "feast" "feast_sparse")
foreach(LIB IN LISTS FEAST_LIBS)
  LIST(APPEND LIBS "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${LIB}${CMAKE_STATIC_LIBRARY_SUFFIX}")
endforeach()
SET(FEAST_LIBRARY ${LIBS} CACHE FILEPATH "FEAST library.")
MARK_AS_ADVANCED(FEAST_LIBRARY)

#-------------------------------------------------------------------------------
# The FEAST external project
#-------------------------------------------------------------------------------
# determine which feast library to copy over in the install step
if(USE_FEAST_COMMUNITY_PRECOMPILED)
  set(FEAST_LIB_DIR "${feast_source}/${FEAST_VER}/lib/x64")
else(USE_FEAST_COMMUNITY_PRECOMPILED)
  set(FEAST_LIB_DIR "${feast_source}/${FEAST_VER}/lib/${CFS_ARCH_STR}")
endif(USE_FEAST_COMMUNITY_PRECOMPILED)

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(feast
    PREFIX "${feast_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${FEAST_LIB_DIR} ${LD}
    BUILD_BYPRODUCTS ${FEAST_LIBRARY}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    # FEAST does not use CMake but make only: We cannot use the CMake Ninja-generator,
    # we use make instead to build FEAST
    find_program(FEAST_MAKE_PROGRAM make)
  else("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    set(FEAST_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "program to build FEAST")
  endif("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  MARK_AS_ADVANCED(FEAST_MAKE_PROGRAM)
  # here we should decide if we need to build feast, or only copy the pre-compiled libs a below
  if(USE_FEAST_COMMUNITY_PRECOMPILED)
  # only copy over
  ExternalProject_Add(feast
    PREFIX "${feast_prefix}"
    SOURCE_DIR "${feast_source}"
    URL ${LOCAL_FILE}
    URL_MD5 "${FEAST_MD5}"
    BINARY_DIR "${feast_source}/${FEAST_VER}/src"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory ${FEAST_LIB_DIR} ${LD}
    BUILD_BYPRODUCTS ${FEAST_LIBRARY}
  )
  else(USE_FEAST_COMMUNITY_PRECOMPILED)
  # compile feast
  ExternalProject_Add(feast
    PREFIX "${feast_prefix}"
    SOURCE_DIR "${feast_source}"
    URL ${LOCAL_FILE}
    URL_MD5 "${FEAST_MD5}"
    BINARY_DIR "${feast_source}/${FEAST_VER}/src"
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy "${CONF}" "${feast_source}/${FEAST_VER}/src" # copy over file
    BUILD_COMMAND ${FEAST_MAKE_PROGRAM} "ARCH=${CFS_ARCH_STR}" "LIB=feast" "all"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory "${FEAST_LIB_DIR}" ${LD}
    BUILD_BYPRODUCTS ${FEAST_LIBRARY}
  )
  endif(USE_FEAST_COMMUNITY_PRECOMPILED)

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(feast cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${feast_prefix}
  )

  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(feast cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
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
  feast
)

#-----------------------------------------------------------------------------
# Determine paths of FEAST includes.
#-----------------------------------------------------------------------------
SET(FEAST_INCLUDE_DIR "${feast_source}/${FEAST_VER}/include" CACHE PATH "include path for FEAST")
MARK_AS_ADVANCED(FEAST_INCLUDE_DIR)
