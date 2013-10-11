#-------------------------------------------------------------------------------
# muParser - a fast math parser library
#
# Project Homepage
# http://muparser.sourceforge.net/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to muparser sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(muparser_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/muparser")
set(muparser_source  "${muparser_prefix}/src/muparser")
set(muparser_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${muparser_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
)

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://pkgs.fedoraproject.org/repo/pkgs/muParser/muparser_v2_2_2.zip/6d77b5cb8096fe2c50afe36ad41bc14a/muparser_v2_2_2.zip"
  "ftp://ftp.jp.netbsd.org/pub/pkgsrc/distfiles/muparser_v2_2_2.zip"
  "${MUPARSER_URL}/${MUPARSER_ZIP}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/muparser/${MUPARSER_ZIP}")
SET(MD5_SUM ${MUPARSER_MD5})

SET(DLFN "${muparser_prefix}/muparser-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
  )

#-------------------------------------------------------------------------------
# The muparser external project
#-------------------------------------------------------------------------------
ExternalProject_Add(muparser
  PREFIX "${muparser_prefix}"
  URL ${LOCAL_FILE}
  URL_MD5 ${MUPARSER_MD5}
  CMAKE_ARGS
    ${CMAKE_ARGS}
)

#-------------------------------------------------------------------------------
# Add custom download step to be able to download from a list of mirrors
# instead of just a single URL.
#-------------------------------------------------------------------------------
ExternalProject_Add_Step(muparser cfsdeps_download
   COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
   DEPENDERS download
   DEPENDS "${DLFN}"
   WORKING_DIRECTORY ${muparser_prefix}
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
