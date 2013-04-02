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
# The qt4 external project
#-------------------------------------------------------------------------------
ExternalProject_Add(qt4
  PREFIX ${qt4_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/qt4
  SOURCE_DIR ${qt4_source}
  URL ${QT4_URL}/${QT4_GZ}
  URL_MD5 ${QT4_MD5}
  CONFIGURE_COMMAND ${Qt4_CONFIGURE_COMMAND}
  BUILD_COMMAND ${Qt4_BUILD_COMMAND}
  INSTALL_COMMAND ${Qt4_INSTALL_COMMAND}
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/qt4/qt4-patch.cmake.in")
SET(PFN "${qt4_prefix}/qt4-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# We do not use the PATCH_COMMAND  of ExternalProject_Add since we do not only
# want to apply the patch script  during configuration time but also if it has
# changed.  Therefore,   we  need  a   dependency  on  the   configured  patch
# script. This can be achieved by  adding an additional build step between the
# download and configure steps.
#
# NOTE: The  patch script should  be designed  in such a  way, that it  can be
# applied to  an already patched  source tree. This  is due to the  fact, that
# ExternalProject_Add only extracts the source if the MD5 sum has has changed.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(qt4 custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${qt4_source}
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


