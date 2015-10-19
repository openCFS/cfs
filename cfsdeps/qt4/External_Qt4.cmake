#-------------------------------------------------------------------------------
# Qt Cross-Platform Application Framework
#
# Project Homepage
# http://qt.digia.com/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to qt4 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(qt4_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/qt4")
set(qt4_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(qt4_source  "${qt4_prefix}/src/qt4")

configure_file("${CFS_SOURCE_DIR}/cfsdeps/qt4/qt4-configure-step.cmake.in"
  "${qt4_prefix}/qt4-configure-step.cmake"
  @ONLY)

set(Qt4_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${qt4_prefix}/qt4-configure-step.cmake)

configure_file("${CFS_SOURCE_DIR}/cfsdeps/qt4/qt4-build-step.cmake.in"
  "${qt4_prefix}/qt4-build-step.cmake"
  @ONLY)

set(Qt4_BUILD_COMMAND ${CMAKE_COMMAND} -P ${qt4_prefix}/qt4-build-step.cmake)

configure_file("${CFS_SOURCE_DIR}/cfsdeps/qt4/qt4-install-step.cmake.in"
  "${qt4_prefix}/qt4-install-step.cmake"
  @ONLY)

set(Qt4_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${qt4_prefix}/qt4-install-step.cmake)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/qt4/qt4-patch.cmake.in")
SET(PFN "${qt4_prefix}/qt4-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://distfiles.macports.org/qt4-mac/${QT4_GZ}"
  "http://pkgs.fedoraproject.org/repo/pkgs/qt/${QT4_GZ}/${QT4_MD5}/${QT4_GZ}"
  "http://download.qt-project.org/archive/qt/4.8/4.8.2/${QT4_GZ}"
  "${QT4_URL}/${QT4_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/qt4/${QT4_GZ}")
SET(MD5_SUM ${QT4_MD5})

SET(DLFN "${qt4_prefix}/qt4-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "qt4" "${QT4_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${qt4_prefix}")

SET(ZIPFROMCACHE "${qt4_prefix}/qt4-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${qt4_prefix}/qt4-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The qt4 external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(qt4
    PREFIX "${qt4_prefix}"
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
  ExternalProject_Add(qt4
    PREFIX ${qt4_prefix}
    SOURCE_DIR ${qt4_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${QT4_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CONFIGURE_COMMAND ${Qt4_CONFIGURE_COMMAND}
    BUILD_COMMAND ${Qt4_BUILD_COMMAND}
    INSTALL_COMMAND ${Qt4_INSTALL_COMMAND}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(qt4 cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${qt4_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(qt4 cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

SET(QT_QMAKE_EXECUTABLE "${CFS_BINARY_DIR}/qt4/bin/qmake")
IF(WIN32 AND NOT MINGW)
  SET($ENV{PATH} "${CFS_BINARY_DIR}/qt4/bin;$ENV{PATH}")
ELSE()
  SET($ENV{PATH} "${CFS_BINARY_DIR}/qt4/bin:$ENV{PATH}")
ENDIF()

LIST(APPEND CFS_PV_DEPENDENCIES qt4)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  qt4
)


