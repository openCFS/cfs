#-------------------------------------------------------------------------------
# Set paths to gidpost sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(gidpost_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/gidpost")
set(gidpost_source  "${gidpost_prefix}/src/gidpost")
set(gidpost_install  "${CFS_BINARY_DIR}")

#-------------------------------------------------------------------------------
# The GiDpost external project
#-------------------------------------------------------------------------------
ExternalProject_Add(gidpost
  PREFIX "${gidpost_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/gidpost
  URL ${GIDPOST_URL}/${GIDPOST_ZIP}
  URL_MD5 ${GIDPOST_MD5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${gidpost_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
    -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  )

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/gidpost/gidpost-patch.cmake.in")
SET(PFN "${gidpost_prefix}/gidpost-patch.cmake")
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
ExternalProject_Add_Step(gidpost custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${gidpost_source}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  gidpost
)

#-------------------------------------------------------------------------------
# These variables are used to find Gidpost by other projects
#-------------------------------------------------------------------------------
SET(GIDPOST_INCLUDE_DIR "${CFS_BINARY_DIR}/include" CACHE FILEPATH "GiDpost include dir" FORCE)

#-------------------------------------------------------------------------------
# Determine paths of GIDPOST libraries.
#-------------------------------------------------------------------------------
SET(GIDPOST_LIBRARY_DEBUG
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)
SET(GIDPOST_LIBRARY_RELEASE
  "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/libgidpost.a"
  CACHE FILEPATH "GiDpost library" FORCE)

#-------------------------------------------------------------------------------
# Mark paths of GIDPOST libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(GIDPOST_LIBRARY_DEBUG)
MARK_AS_ADVANCED(GIDPOST_LIBRARY_RELEASE)

#-------------------------------------------------------------------------------
# Set GIDPOST_LIBRARY according to configuration
#-------------------------------------------------------------------------------
IF(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_DEBUG}")
ELSE(DEBUG)
  SET(GIDPOST_LIBRARY "${GIDPOST_LIBRARY_RELEASE}")
ENDIF(DEBUG)

SET(CFS_GIDPOST_VERSION "1.71")
