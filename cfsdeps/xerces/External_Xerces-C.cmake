#-------------------------------------------------------------------------------
# Xerces-C++ is a validating XML parser written in a portable subset of C++.
#
# Project Homepage
# http://xerces.apache.org/xerces-c/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to xerces sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(xerces_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/xerces")
set(xerces_source  "${xerces_prefix}/src/xerces")
set(xerces_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-configure.cmake.in")
SET(CONF "${xerces_prefix}/xerces-configure.cmake")
CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of build file and template file.
#-------------------------------------------------------------------------------
SET(BUILD_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-build.cmake.in")
SET(BUILD "${xerces_prefix}/xerces-build.cmake")
CONFIGURE_FILE("${BUILD_TEMPL}" "${BUILD}" @ONLY) 

#-------------------------------------------------------------------------------
# Set names of configure file and template file.
#-------------------------------------------------------------------------------
SET(INST_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-install.cmake.in")
SET(INST "${xerces_prefix}/xerces-install.cmake")
CONFIGURE_FILE("${INST_TEMPL}" "${INST}" @ONLY) 

#-------------------------------------------------------------------------------
# The xerces external project
#-------------------------------------------------------------------------------
ExternalProject_Add(xerces
  PREFIX "${xerces_prefix}"
  SOURCE_DIR "${xerces_source}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/xerces
  URL ${XERCES_URL}/${XERCES_GZ}
  URL_MD5 ${XERCES_MD5}
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ${CMAKE_COMMAND} -P ${CONF}
  BUILD_COMMAND  ${CMAKE_COMMAND} -P ${BUILD}
  INSTALL_COMMAND  ${CMAKE_COMMAND} -P ${INST}
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-patch.cmake.in")
#SET(PFN "${xerces_prefix}/xerces-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

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
#ExternalProject_Add_Step(xerces custom_patch
#   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
#   DEPENDEES download
#   DEPENDERS configure
#   DEPENDS "${PFN}"
#   WORKING_DIRECTORY ${xerces_source}
#)


SET(XERCES_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-----------------------------------------------------------------------------
# Determine paths of XERCES libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(XERCES_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "XERCES library.")

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(XERCES_LIBRARY "${XERCES_LIBRARY};-framework CoreFoundation;-framework CoreServices")
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

MARK_AS_ADVANCED(XERCES_LIBRARY)
MARK_AS_ADVANCED(XERCES_INCLUDE_DIR)

