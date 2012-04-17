#-------------------------------------------------------------------------------
# Set paths to metis sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(metis_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/metis")
set(metis_source  "${metis_prefix}/src/metis")
set(metis_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# The metis external project
#-------------------------------------------------------------------------------
ExternalProject_Add(metis
  PREFIX "${metis_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/metis
  URL ${METIS_URL}/${METIS_GZ}
  URL_MD5 ${METIS_MD5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${metis_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/metis/metis-patch.cmake.in")
SET(PFN "${metis_prefix}/metis-patch.cmake")
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
ExternalProject_Add_Step(metis custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${metis_source}
)

SET(METIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of METIS libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(METIS_LIBRARY
  "${LD}/libmetis.a"
  CACHE FILEPATH "METIS library.")

MARK_AS_ADVANCED(METIS_INCLUDE_DIR)
MARK_AS_ADVANCED(METIS_LIBRARY)
