#-------------------------------------------------------------------------------
# CFD General Notation System (CGNS)
# Needed for ADF routines by STARCCM+ reader. 
# To build cgnsview which can be used to view .cgns files (HDF5 and ADF) and .ccm files (ADF)
# install it manually or reenable BUILD_CGNSTOOL.
#
# Project Homepage
# http://www.cgns.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to CGNS sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cgns_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgns")
set(cgns_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(cgns_source  "${cgns_prefix}/src/cgns")

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${cgns_install}
  # CMAKE_INSTALL_LIBDIR is ignored (4.3.0), hence we need to patch
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_BUILD_TYPE:STRING=Release
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCGNS_BUILD_CGNSTOOLS:BOOL=OFF # handles only cgnstools, we remove tools by patch
  -DCGNS_BUILD_SHARED:BOOL=OFF
  -DCGNS_BUILD_TESTING:BOOL=OFF
  -DCGNS_ENABLE_64BIT:BOOL=OFF
  -DCGNS_ENABLE_BASE_SCOPE:BOOL=OFF
  -DCGNS_ENABLE_FORTRAN:BOOL=OFF
  -DCGNS_ENABLE_HDF5:BOOL=ON
  -DCGNS_ENABLE_TESTS:BOOL=OFF
  -DCGNS_USE_SHARED:BOOL=OFF
  # we have no hdf5 lib and include dir but need hdf5-config.cmake
  -DHDF5_DIR:FILEPATH=${CMAKE_CURRENT_BINARY_DIR}/share/cmake # for hdf5-config.cmake
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFLAGS}
)

#intel compiler has problems with emty rpath commands
IF(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  SET(CMAKE_ARGS ${CMAKE_ARGS} -DCMAKE_SKIP_RPATH:BOOL=TRUE)
ENDIF()

IF(CFS_DISTRO STREQUAL "MACOSX")
  # Explicitly set build architectures and  system SDK root dir to match
  # those given in CMake.
  
  SET(CMAKE_ARGS
    ${CMAKE_ARGS}
    -DCMAKE_OSX_SYSROOT:PATH=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_ARCHITECTURES:STRING=${CMAKE_OSX_ARCHITECTURES}
  )
ENDIF(CFS_DISTRO STREQUAL "MACOSX")

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgns/cgns-patch.cmake.in")
SET(PFN "${cgns_prefix}/cgns-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "https://github.com/CGNS/CGNS/archive/refs/tags/${CGNS_GZ}"
  "${CFS_DS_SOURCES_DIR}/cgns/${CGNS_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgns/${CGNS_GZ}")
SET(MD5_SUM ${CGNS_MD5})

SET(DLFN "${cgns_prefix}/cgns-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/cgns/license/" DESTINATION "${CFS_BINARY_DIR}/license/cgns" )

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "cgns" "${CGNS_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${cgns_prefix}")

SET(ZIPFROMCACHE "${cgns_prefix}/cgns-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${cgns_prefix}/cgns-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of CGNS libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
SET(CGNS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}cgns${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "CGNS library.")

SET(CGNS_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "CGNS include directory")

MARK_AS_ADVANCED(CGNS_INCLUDE_DIR)
MARK_AS_ADVANCED(CGNS_LIBRARY)

#-------------------------------------------------------------------------------
# The CGNS external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cgns
    PREFIX "${cgns_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${CGNS_LIBRARY}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cgns
    DEPENDS hdf5 zlib
    PREFIX ${cgns_prefix}
    SOURCE_DIR ${cgns_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${CGNS_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR ,
    CMAKE_ARGS
      ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${CGNS_LIBRARY}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(cgns cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
     DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${cgns_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(cgns cfsdeps_zipToCache
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
set(CFSDEPS ${CFSDEPS} cgns)
