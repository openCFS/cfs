#-------------------------------------------------------------------------------
# Set prefix path and path to hdf5 sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(hdf5_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/hdf5")
set(hdf5_install  "${CMAKE_CURRENT_BINARY_DIR}")
set(hdf5_source  "${hdf5_prefix}/src/hdf5")

# intel compiler 2018 are not compatible with glibc 2.27, see
# https://software.intel.com/en-us/forums/intel-c-compiler/topic/777003
set(hdf5_c_flags "${CFSDEPS_C_FLAGS}")
set(hdf5_cxx_flags ${CFS_SUPPRESSIONS})
if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel") 
  IF(NOT WIN32)
    set(hdf5_c_flags "-D_Float32=float -D_Float64='long double' -D_Float32x=double -D_Float64x='long double' ${hdf5_c_flags}")
  ENDIF()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # Apple clang version 12.0.0 needs -Wno-implicit-function-declaration
  # Linux clang version 15.0.4 needs -Wno-int-conversion
  set(hdf5_c_flags " -Wno-implicit-function-declaration -Wno-int-conversion")
  set(hdf5_cxx_flags " -Wno-implicit-function-declaration -Wno-int-conversion")
endif()

#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${hdf5_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DHDF5_INSTALL_BIN_DIR:PATH=bin
    -DHDF5_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}
    # We do not want to see warning messages from external projects
    -DCMAKE_C_FLAGS:STRING=${hdf5_c_flags}
    -DCMAKE_CXX_FLAGS:STRING=${hdf5_cxx_flags}
    -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  )
ELSE(WIN32)
  SET(CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${hdf5_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  #  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${hdf5_install}/bin
  #  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${hdf5_install}/${LIB_SUFFIX}
  #  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${hdf5_install}/${LIB_SUFFIX}
    -DHDF5_INSTALL_BIN_DIR:PATH=bin
    -DHDF5_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}
    # We do not want to see warning messages from external projects
    -DCMAKE_C_FLAGS:STRING=${hdf5_c_flags}
    -DCMAKE_CXX_FLAGS:STRING=${hdf5_cxx_flags}
    -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  )
ENDIF(WIN32)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

SET(HDF5_BUILD_TYPE static)

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/hdf5/hdf5-patch.cmake.in")
SET(PFN "${hdf5_prefix}/hdf5-${HDF5_BUILD_TYPE}-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the openCFS development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------

SET(MIRRORS
  "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-${HDF5_VER}/src/${HDF5_BZ2}"
  "https://src.fedoraproject.org/repo/pkgs/hdf5/${HDF5_BZ2}"
  "${CFS_DS_SOURCES_DIR}/hdf5/${HDF5_BZ2}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/hdf5/${HDF5_BZ2}")
SET(MD5_SUM ${HDF5_MD5})

SET(DLFN "${hdf5_prefix}/hdf5-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/hdf5/license/" DESTINATION "${CFS_BINARY_DIR}/license/hdf5" )

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "hdf5" "${HDF5_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${hdf5_prefix}")

SET(ZIPFROMCACHE "${hdf5_prefix}/hdf5-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${hdf5_prefix}/hdf5-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of HDF5 libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(HDF5_LIBRARY_RELEASE ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5)
  SET(HDF5_CPP_LIBRARY_RELEASE ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_cpp)
  SET(HDF5_LIBRARY_DEBUG ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_D)
  SET(HDF5_CPP_LIBRARY_DEBUG ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_cpp_D)
  SET(HDF5_LT_LIBRARY_RELEASE ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_hl)
  SET(HDF5_LT_CPP_LIBRARY_RELEASE ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_hl_cpp)
  SET(HDF5_LT_LIBRARY_DEBUG ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_hl_D)
  SET(HDF5_LT_CPP_LIBRARY_DEBUG ${CFS_BINARY_DIR}/${LIB_SUFFIX}/hdf5_hl_cpp_D)
ELSE(WIN32)
  SET(HDF5_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_debug")
  SET(HDF5_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5")

  SET(HDF5_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp_debug")
  SET(HDF5_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_cpp")

  SET(HDF5_LT_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_debug")
  SET(HDF5_LT_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl")

  SET(HDF5_LT_CPP_LIBRARY_DEBUG
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp_debug")
  SET(HDF5_LT_CPP_LIBRARY_RELEASE
    "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}hdf5_hl_cpp")

ENDIF(WIN32)

#-------------------------------------------------------------------------------
# Set HDF5_LIBRARY* according to configuration
#-------------------------------------------------------------------------------
SET(HDF5_SHARED_LIBRARY_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})

SET(HDF5_SHARED_LIBRARY_RELEASE ${HDF5_LIBRARY_RELEASE})
SET(HDF5_CPP_SHARED_LIBRARY_RELEASE ${HDF5_CPP_LIBRARY_RELEASE})
SET(HDF5_LT_SHARED_LIBRARY_RELEASE ${HDF5_LT_LIBRARY_RELEASE})
SET(HDF5_LT_CPP_SHARED_LIBRARY_RELEASE ${HDF5_LT_CPP_LIBRARY_RELEASE})

SET(HDF5_SHARED_LIBRARY_DEBUG ${HDF5_LIBRARY_DEBUG})
SET(HDF5_CPP_SHARED_LIBRARY_DEBUG ${HDF5_CPP_LIBRARY_DEBUG})
SET(HDF5_LT_SHARED_LIBRARY_DEBUG ${HDF5_LT_LIBRARY_DEBUG})
SET(HDF5_LT_CPP_SHARED_LIBRARY_DEBUG ${HDF5_LT_CPP_LIBRARY_DEBUG})

IF(WIN32)
  SET(HDF5_SHARED_LIBRARY_SUFFIX ${CMAKE_IMPORT_LIBRARY_SUFFIX})
ENDIF()

IF(DEBUG)
  SET(HDF5_SHARED_LIBRARY "${HDF5_SHARED_LIBRARY_DEBUG}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_CPP_SHARED_LIBRARY "${HDF5_CPP_SHARED_LIBRARY_DEBUG}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_LT_SHARED_LIBRARY "${HDF5_LT_SHARED_LIBRARY_DEBUG}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_LT_CPP_SHARED_LIBRARY "${HDF5_LT_CPP_SHARED_LIBRARY_DEBUG}${HDF5_SHARED_LIBRARY_SUFFIX}")

  SET(HDF5_LIBRARY "${HDF5_LIBRARY_DEBUG}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_DEBUG}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_LIBRARY "${HDF5_LT_LIBRARY_DEBUG}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_CPP_LIBRARY "${HDF5_LT_CPP_LIBRARY_DEBUG}${CMAKE_STATIC_LIBRARY_SUFFIX}")
ELSE(DEBUG)
  SET(HDF5_LIBRARY "${HDF5_LIBRARY_RELEASE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_CPP_LIBRARY "${HDF5_CPP_LIBRARY_RELEASE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_LIBRARY "${HDF5_LT_LIBRARY_RELEASE}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(HDF5_LT_CPP_LIBRARY "${HDF5_LT_CPP_LIBRARY_RELEASE}${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(HDF5_SHARED_LIBRARY "${HDF5_SHARED_LIBRARY_RELEASE}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_CPP_SHARED_LIBRARY "${HDF5_CPP_SHARED_LIBRARY_RELEASE}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_LT_SHARED_LIBRARY "${HDF5_LT_SHARED_LIBRARY_RELEASE}${HDF5_SHARED_LIBRARY_SUFFIX}")
  SET(HDF5_LT_SHARED_CPP_LIBRARY "${HDF5_LT_CPP_SHARED_LIBRARY_RELEASE}${HDF5_SHARED_LIBRARY_SUFFIX}")
ENDIF(DEBUG)

set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5_CPP_LIBRARY} ${HDF5_LT_LIBRARY} ${HDF5_LT_CPP_LIBRARY} ${HDF5_SHARED_LIBRARY} ${HDF5_CPP_SHARED_LIBRARY} ${HDF5_LT_SHARED_LIBRARY} ${HDF5_LT_SHARED_CPP_LIBRARY})

#-------------------------------------------------------------------------------
# The hdf5-static external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(hdf5-static
    PREFIX "${hdf5_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  ExternalProject_Add(hdf5-shared
    PREFIX "${hdf5_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${HDF5_LIBRARIES}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(hdf5-static
    DEPENDS zlib
    PREFIX ${hdf5_prefix}
    SOURCE_DIR ${hdf5_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${HDF5_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR ^
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DHDF5_EXTERNAL_LIB_PREFIX:STRING=
      -DHDF5_BUILD_CPP_LIB:BOOL=ON
      -DHDF5_BUILD_HL_LIB:BOOL=ON
      -DHDF5_BUILD_FORTRAN:BOOL=OFF
      -DHDF5_BUILD_EXAMPLES:BOOL=OFF
      -DHDF5_BUILD_TOOLS:BOOL=OFF
      -DHDF5_ENABLE_HSIZET:BOOL=OFF
      -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_LIBRARY}
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DHDF5_BUILD_TOOLS:BOOL=OFF # isn't it sufficient to use system hdf5 tools?!
      # On Mac OS X we can get problems with the system strdup function.
      -DH5_HAVE_STRDUP:BOOL=OFF
    BUILD_BYPRODUCTS
      ${HDF5_LIBRARIES}
    )
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(hdf5-static cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${hdf5_prefix}
  )
  
  SET(HDF5_BUILD_TYPE shared)
  
  #-------------------------------------------------------------------------------
  # Set names of patch file and template file.
  #-------------------------------------------------------------------------------
  SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/hdf5/hdf5-patch.cmake.in")
  SET(PFN "${hdf5_prefix}/hdf5-${HDF5_BUILD_TYPE}-patch.cmake")
  CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 
  
  #-------------------------------------------------------------------------------
  # The hdf5-shared external project
  # We need it to build a shared version of CGNS library 3.1.
  #-------------------------------------------------------------------------------
  ExternalProject_Add(hdf5-shared
    DEPENDS zlib hdf5-static
    PREFIX ${hdf5_prefix}
    DOWNLOAD_COMMAND ""
    SOURCE_DIR ${hdf5_source}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=ON
      -DBUILD_TESTING:BOOL=OFF
      -DHDF5_BUILD_CPP_LIB:BOOL=ON
      -DHDF5_BUILD_HL_LIB:BOOL=ON
      -DHDF5_BUILD_FORTRAN:BOOL=OFF
      -DHDF5_BUILD_EXAMPLES:BOOL=OFF
      -DHDF5_BUILD_TOOLS:BOOL=OFF
      -DHDF5_ENABLE_HSIZET:BOOL=OFF
      -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:FILEPATH=${ZLIB_SHARED_LIBRARY}
  )
    
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(hdf5-shared cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

ExternalProject_Get_Property(hdf5-shared source_dir)

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  hdf5-static
  hdf5-shared
)



IF(UNIX)
  LIST(APPEND HDF5_LIBRARY -ldl)
ENDIF()

SET(HDF5_INCLUDE_DIR ${CFS_BINARY_DIR}/include CACHE PATH "hdf5 include directory")

INCLUDE_DIRECTORIES("${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(HDF5_INCLUDE_DIR)
