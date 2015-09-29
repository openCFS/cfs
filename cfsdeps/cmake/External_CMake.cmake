#-------------------------------------------------------------------------------
# CMake, the cross-platform, open-source build system.
# Needed by HDF5 and ParaView.
#
# Project Homepage
# http://www.cmake.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to cmake sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cmake_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cmake")
set(cmake_source  "${cmake_prefix}/src/cmake")
set(cmake_install  "${CFS_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cmake/cmake-configure.cmake.in")
SET(CONF "${cmake_prefix}/cmake-configure.cmake")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(INST_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cmake/cmake-install.cmake.in")
SET(INST "${cmake_prefix}/cmake-install.cmake")
CONFIGURE_FILE("${INST_TEMPL}" "${INST}" @ONLY) 

IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "cmake_${CMAKE_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}_${CMAKE_BUILD_TYPE}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "cmake_${CMAKE_VER}_${CFS_ARCH_STR}_${FC_ID}_${CMAKE_BUILD_TYPE}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${cmake_prefix}")

SET(ZIPFROMCACHE "${cmake_prefix}/cmake-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${cmake_prefix}/cmake-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The cmake external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cmake
    PREFIX "${cmake_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
ELSE("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cmake
    PREFIX "${cmake_prefix}"
    SOURCE_DIR "${cmake_source}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/cmake
    URL ${CMAKE_URL}/${CMAKE_GZ}
    URL_MD5 ${CMAKE_MD5}
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
    BUILD_COMMAND  ${CMAKE_MAKE_PROGRAM}
    INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INST}
  )
  
  IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(cmake cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

LIST(APPEND CFS_PV_DEPENDENCIES cmake)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  cmake
)
