#-------------------------------------------------------------------------------
# FLANN - Fast Library for Approximate Nearest Neighbors
# 
# Needed as a BSD-licensed alternative to CGAL for finding nearest neighbors
# in CoefFunctionScatteredData.
#
# Project Homepage
# http://www.cs.ubc.ca/research/flann/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to flann sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(flann_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/flann")
set(flann_source  "${flann_prefix}/src/flann")
set(flann_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/flann/flann-patch.cmake.in")
SET(PFN "${flann_prefix}/flann-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

STRING(REPLACE ";" "," FLANN_HDF5_LIBRARY "${HDF5_LIBRARY};${ZLIB_LIBRARY}")

# set flann suffix
string(REPLACE "lib" "" FLANN_LIB_SUFFIX "${LIB_SUFFIX}")
#-------------------------------------------------------------------------------
# Set common CMake arguments
#-------------------------------------------------------------------------------
SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${flann_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCFS_ARCH_STR:STRING=${CFS_ARCH_STR}
  -DLIB_SUFFIX:STRING=${FLANN_LIB_SUFFIX}/${CFS_ARCH_STR}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DUSE_OPENMP:BOOL=${USE_OPENMP}
  -DBUILD_DOC:BOOL=OFF
  -DHDF5_DIR:FILEPATH=${CFS_BINARY_DIR}/cmake/hdf5
  -DHDF5_C_LIBRARY:PATH=${FLANN_HDF5_LIBRARY}
  -DHDF5_C_HL_LIBRARY:PATH=${HDF5_LT_LIBRARY}
  -DHDF5_CXX_LIBRARY:PATH=${HDF5_CPP_LIBRARY}
  -DHDF5_CXX_HL_LIBRARY:PATH=${HDF5_LT_CPP_LIBRARY}
  -DHDF5_INCLUDE_DIR:FILEPATH=${CFS_BINARY_DIR}/include
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
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/flann/flann-patch.cmake.in")
SET(PFN "${flann_prefix}/flann-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "https://github.com/mariusmuja/flann/archive/${FLANN_VER}.zip"
  #"http://distfiles.macports.org/flann/${FLANN_ZIP}"
  #"http://pkgs.fedoraproject.org/repo/pkgs/flann/${FLANN_ZIP}/${FLANN_MD5}/${FLANN_ZIP}"
  #"http://www.cs.ubc.ca/research/flann/uploads/FLANN/${FLANN_ZIP}"
  "${FLANN_URL}/${FLANN_ZIP}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/flann/${FLANN_ZIP}")
SET(MD5_SUM ${FLANN_MD5})

SET(DLFN "${flann_prefix}/flann-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/flann/license/" DESTINATION "${CFS_BINARY_DIR}/license/flann" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "flann" "${FLANN_VER}") 
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${flann_prefix}")

SET(ZIPFROMCACHE "${flann_prefix}/flann-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${flann_prefix}/flann-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of FLANN libraries.
#-------------------------------------------------------------------------------
IF(UNIX)
  SET(FLANN_LIB flann_cpp_s)
  SET(FLANN_SHARED_LIB flann_cpp_s)
ELSE()
  SET(FLANN_LIB flann_cpp_s)
  SET(FLANN_SHARED_LIB flann)
  IF(DEBUG)
    SET(FLANN_LIB "${FLANN_LIB}d")
    SET(FLANN_SHARED_LIB "${FLANN_SHARED_LIB}d")
  ENDIF()
ENDIF()

SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
IF(WIN32)
  SET(FLANN_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${FLANN_LIB}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "flann library")
  SET(FLANN_SHARED_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${FLANN_SHARED_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
ELSE(WIN32)
  SET(FLANN_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}${FLANN_LIB}${CMAKE_STATIC_LIBRARY_SUFFIX} CACHE FILEPATH "flann library")
  SET(FLANN_SHARED_LIBRARY ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}/${CMAKE_STATIC_LIBRARY_PREFIX}${FLANN_SHARED_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
ENDIF(WIN32)
SET(FLANN_SHARED_LIBRARY ${FLANN_SHARED_LIBRARY} CACHE FILEPATH "flann shared library")
SET(FLANN_INCLUDE_DIR "${CFS_BINARY_DIR}/include/flann" CACHE PATH "flann include directory")
MARK_AS_ADVANCED(FLANN_LIBRARY)
MARK_AS_ADVANCED(FLANN_SHARED_LIBRARY)
MARK_AS_ADVANCED(FLANN_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The flann external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(flann
    PREFIX "${flann_prefix}"
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
  ExternalProject_Add(flann
    DEPENDS hdf5-static
    PREFIX ${flann_prefix}
    SOURCE_DIR ${flann_source}
    URL ${LOCAL_FILE}
    URL_MD5 ${FLANN_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    LIST_SEPARATOR ,
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DBUILD_CUDA_LIB:BOOL=OFF
      -DBUILD_C_BINDINGS=OFF
      -DBUILD_MATLAB_BINDINGS=OFF
      -DBUILD_PYTHON_BINDINGS=OFF
      -DPYTHON_EXECUTABLE=PYTHON_EXECUTABLE_NOTFOUND
    BUILD_BYPRODUCTS ${FLANN_LIBRARY} ${FLANN_SHARED_LIBRARY}
  )
  
  ADD_DEPENDENCIES(flann zlib)
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(flann cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${flann_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(flann cfsdeps_zipToCache
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
  flann
)

