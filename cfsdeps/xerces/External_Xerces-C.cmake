#-------------------------------------------------------------------------------
# Xerces-C++ is a validating XML parser written in a portable subset of C++.
#
# Project Homepage
# http://xerces.apache.org/xerces-c/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to Xerces-C sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(xerces_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/xerces")
set(xerces_source  "${xerces_prefix}/src/xerces")
set(xerces_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${xerces_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  )

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_CXX_FLAGS:FILEPATH=${CFLAGS}
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
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/xerces/xerces-patch.cmake.in")
SET(PFN "${xerces_prefix}/xerces-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://archive.apache.org/dist/xerces/c/3/sources/${XERCES_GZ}"
  "http://xml.apache.org/dist/xerces-c/3/sources/${XERCES_GZ}"
  "http://archive.apache.org/dist/xerces/c/3/sources/${XERCES_GZ}"
  "${CFS_DS_SOURCES_DIR}/xerces/${XERCES_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/xerces/${XERCES_GZ}")
SET(MD5_SUM ${XERCES_MD5})

SET(DLFN "${xerces_prefix}/xerces-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/xerces/license/" DESTINATION "${CFS_BINARY_DIR}/license/xerces" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "xerces" "${XERCES_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${xerces_prefix}")

SET(ZIPFROMCACHE "${xerces_prefix}/xerces-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${xerces_prefix}/xerces-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of XERCES libraries.
#-----------------------------------------------------------------------------
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")

IF(CFS_DISTRO STREQUAL "MACOSX")
  SET(XERCES_LIBRARY "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX}")
  LIST(APPEND XERCES_LIBRARY "-framework CoreFoundation")
  LIST(APPEND XERCES_LIBRARY "-framework CoreServices")
ELSEIF(UNIX)
  SET(XERCES_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX};-lpthread"
    CACHE FILEPATH "XERCES library.")
ELSE()
  SET(XERCES_LIBRARY
    "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}xerces-c${CMAKE_STATIC_LIBRARY_SUFFIX}"
    CACHE FILEPATH "XERCES library.")
ENDIF()

#-------------------------------------------------------------------------------
# The xerces external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(xerces
    PREFIX "${xerces_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${XERCES_LIBRARY}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(xerces
    PREFIX "${xerces_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${XERCES_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${XERCES_LIBRARY}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(xerces cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${xerces_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(xerces cfsdeps_zipToCache
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
SET(CFSDEPS ${CFSDEPS} xerces)

SET(XERCES_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(XERCES_LIBRARY)
MARK_AS_ADVANCED(XERCES_INCLUDE_DIR)

