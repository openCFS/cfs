#-------------------------------------------------------------------------------
# Set paths to muparser sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(muparser_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser")
set(muparser_source  "${muparser_prefix}/src/muparser")
set(muparser_install  "${CMAKE_CURRENT_BINARY_DIR}")

#-------------------------------------------------------------------------------
# The muparser external project
#-------------------------------------------------------------------------------
ExternalProject_Add(muparser
  PREFIX "${muparser_prefix}"
  DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/muparser
  URL ${MUPARSER_URL}/${MUPARSER_ZIP}
  URL_MD5 ${MUPARSER_MD5}
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${muparser_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/muparser/muparser-patch.cmake.in")
SET(PFN "${muparser_prefix}/muparser-patch.cmake")
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
ExternalProject_Add_Step(muparser custom_patch
   COMMAND ${CMAKE_COMMAND} -P "${PFN}"
   DEPENDEES download
   DEPENDERS configure
   DEPENDS "${PFN}"
   WORKING_DIRECTORY ${muparser_source}
)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  muparser
)

SET(MUPARSER_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of MUPARSER libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(MUPARSER_LIBRARY
  "${LD}/libmuparser.a"
  CACHE FILEPATH "MUPARSER library.")

MARK_AS_ADVANCED(MUPARSER_INCLUDE_DIR)
MARK_AS_ADVANCED(MUPARSER_LIBRARY)
