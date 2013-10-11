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
# Set up a list of publicly available mirrors, since lse17 may not be
# accessible from behind firewalls.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://distfiles.macports.org/qt4-mac/qt-everywhere-opensource-src-4.8.2.tar.gz"
  "http://pkgs.fedoraproject.org/repo/pkgs/qt/qt-everywhere-opensource-src-4.8.2.tar.gz/3c1146ddf56247e16782f96910a8423b/qt-everywhere-opensource-src-4.8.2.tar.gz"
  "http://download.qt-project.org/archive/qt/4.8/4.8.2/qt-everywhere-opensource-src-4.8.2.tar.gz"
)

#-------------------------------------------------------------------------------
# Try to download sources into CFSDEPS cache directory.
#-------------------------------------------------------------------------------
DOWNLOAD_CFSDEPS(
  "${CFS_DEPS_CACHE_DIR}/sources/qt4/${QT4_GZ}"
  ${QT4_MD5}
  "${MIRRORS}"
)

#-------------------------------------------------------------------------------
# The qt4 external project
#-------------------------------------------------------------------------------
ExternalProject_Add(qt4
  PREFIX ${qt4_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4
  SOURCE_DIR ${qt4_source}
  URL ${QT4_URL}/${QT4_GZ}
  URL_MD5 ${QT4_MD5}
  PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  CONFIGURE_COMMAND ${Qt4_CONFIGURE_COMMAND}
  BUILD_COMMAND ${Qt4_BUILD_COMMAND}
  INSTALL_COMMAND ${Qt4_INSTALL_COMMAND}
  )


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


