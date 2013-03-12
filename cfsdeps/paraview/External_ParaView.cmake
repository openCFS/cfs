#-------------------------------------------------------------------------------
# Set prefix path and path to paraview sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(paraview_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/paraview")
set(paraview_install  "${CMAKE_CURRENT_BINARY_DIR}/paraview")
set(paraview_source  "${paraview_prefix}/src/paraview")

#-------------------------------------------------------------------------------
# The paraview-download external project
#-------------------------------------------------------------------------------
ExternalProject_Add(paraview-download
  PREFIX ${paraview_prefix}
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/paraview
  SOURCE_DIR ${paraview_source}
  URL ${PARAVIEW_URL}/${PARAVIEW_GZ}
  URL_MD5 ${PARAVIEW_MD5}
  CMAKE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

configure_file("${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-configure-step.cmake.in"
  "${paraview_prefix}/paraview-configure-step.cmake"
  @ONLY)

set(PV_CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${paraview_prefix}/paraview-configure-step.cmake)

configure_file("${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-build-step.cmake.in"
  "${paraview_prefix}/paraview-build-step.cmake"
  @ONLY)

set(PV_BUILD_COMMAND ${CMAKE_COMMAND} -P ${paraview_prefix}/paraview-build-step.cmake)

configure_file("${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-install-step.cmake.in"
  "${paraview_prefix}/paraview-install-step.cmake"
  @ONLY)

set(PV_INSTALL_COMMAND ${CMAKE_COMMAND} -P ${paraview_prefix}/paraview-install-step.cmake)

# MESSAGE(FATAL_ERROR "CFS_PV_CMAKE_COMMAND ${CFS_PV_CMAKE_COMMAND}")
ExternalProject_Add(paraview-superbuild
  DEPENDS ${CFS_PV_DEPENDS} paraview-download
  PREFIX ${paraview_prefix}
  DOWNLOAD_COMMAND ""
  CONFIGURE_COMMAND ${PV_CONFIGURE_COMMAND}
  BUILD_COMMAND ${PV_BUILD_COMMAND}
  INSTALL_COMMAND ""
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/paraview/paraview-patch.cmake.in")
SET(PFN "${paraview_prefix}/paraview-patch.cmake")
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
ExternalProject_Add_Step(paraview-superbuild custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${paraview_source}
)

# ln -s libhdf5_debug.so.1.8.8 libhdf5.so

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  paraview-superbuild
)


