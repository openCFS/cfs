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
  "http://pkgs.fedoraproject.org/repo/pkgs/muParser/${MUPARSER_ZIP}/${MUPARSER_MD5}/${MUPARSER_ZIP}"
  "ftp://ftp.jp.netbsd.org/pub/pkgsrc/distfiles/${MUPARSER_ZIP}"
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

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "muparser" "${muparser_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${muparser_prefix}")

SET(ZIPFROMCACHE "${muparser_prefix}/muparser-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${muparser_prefix}/muparser-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The muparser external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(muparser
    PREFIX "${muparser_prefix}"
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
  ExternalProject_Add(muparser
    DEPENDS boost
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
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(muparser cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

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
