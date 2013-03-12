#-------------------------------------------------------------------------------
# Set prefix path and path to bzip2 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(bzip2_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/bzip2")
set(bzip2_source  "${bzip2_prefix}/src/bzip2-static")

#-------------------------------------------------------------------------------
# The bzip2-static external project
#-------------------------------------------------------------------------------
ExternalProject_Add(bzip2-static
  PREFIX "${bzip2_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/bzip2
  URL ${BZIP2_URL}/${BZIP2_GZ}
  URL_MD5 ${BZIP2_MD5}
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} bzip2 bzip2recover
  INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/bzip2/bzip2-static-patch.cmake.in")
SET(PFN "${bzip2_prefix}/bzip2-static-patch.cmake")
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
ExternalProject_Add_Step(bzip2-static custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${bzip2_source}
)

#-------------------------------------------------------------------------------
# The bzip2-static external project
#-------------------------------------------------------------------------------
set(bzip2_source  "${bzip2_prefix}/src/bzip2-shared")

ExternalProject_Add(bzip2-shared
  PREFIX "${bzip2_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/bzip2
  URL ${BZIP2_URL}/${BZIP2_GZ}
  URL_MD5 ${BZIP2_MD5}
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -f Makefile-libbz2_so
  INSTALL_COMMAND ""
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/bzip2/bzip2-shared-patch.cmake.in")
SET(PFN "${bzip2_prefix}/bzip2-shared-patch.cmake")
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
ExternalProject_Add_Step(bzip2-shared custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${bzip2_source}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  bzip2-static
  bzip2-shared
)

SET(BZIP2_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}bz2${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "bzip2 library")
SET(BZIP2_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "bzip2 include directory")

MARK_AS_ADVANCED(BZIP2_LIBRARY)
MARK_AS_ADVANCED(BZIP2_INCLUDE_DIR)
