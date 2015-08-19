# libfbi - http://mkirchner.github.com/libfbi/
#
# libfbi is a header-only C++ template library that enables the efficient solution of box intersection problems in an arbitrary number of dimensions. The implementation makes heavy use of C++ metaprogramming and variadic template programming techniques. Despite this complexity, the library provides a straightforward and simple interface that allows easy integration into new and existing projects. 

# kdtree - http://code.google.com/p/kdtree/
# kdtree is a simple, easy to use C library for working with kd-trees.

# Kd-trees are an extension of binary search trees to k-dimensional data. They facilitate very fast searching, and nearest-neighbor queries.

# This particular implementation is designed to be efficient and very easy to use. It is completely written in ANSI/ISO C, and thus completely cross-platform. 

# Octree C++ Class Template - http://nomis80.org/code/octree.html
# An octree is a tree structure for storing sparse 3-D data. What you will find here is source code for the easiest and fastest C++ octree class template around. And it's free. You can find more information on octrees in general from Wikipedia.

# http://www.savarese.com/software/libssrckdtree/

#-------------------------------------------------------------------------------
# Set paths to spacepart sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(spacepart_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/spacepart")
set(libfbi_source  "${spacepart_prefix}/src/libfbi")
set(spacepart_install  "${CMAKE_CURRENT_BINARY_DIR}")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
IF(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Debug")
ELSE(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(CMAKE_BUILD_TYPE STREQUAL "debug")

SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${spacepart_install}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

IF(USE_LIBFBI)
  #-----------------------------------------------------------------------------
  # Set names of patch file and template file.
  #-----------------------------------------------------------------------------
  SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/spacepart/libfbi-patch.cmake.in")
  SET(PFN "${boost_prefix}/libfbi-patch.cmake")
  CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

  IF(WIN32)
    SET(PRECOMPILED_PCKG_NAME "libfbi_${LIBFBI_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}.zip")
  ELSE(WIN32)
    SET(PRECOMPILED_PCKG_NAME "libfbi_${LIBFBI_VER}_${CFS_ARCH_STR}_${FC_ID}.zip")
  ENDIF(WIN32)
  SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
    
  SET(PREFIX_DIR "${spacepart_prefix}")
  
  SET(ZIPFROMCACHE "${spacepart_prefix}/libfbi-zipFromCache.cmake")
  CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)
  
  SET(ZIPTOCACHE "${spacepart_prefix}/libfbi-zipToCache.cmake")
  CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)
  
  #-----------------------------------------------------------------------------
  # The fast box intersection library external project
  #-----------------------------------------------------------------------------
  IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
    #-------------------------------------------------------------------------------
    # If precompiled package exists copy files from cache
    #-------------------------------------------------------------------------------
    ExternalProject_Add(libfbi
      PREFIX "${spacepart_prefix}"
      DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
      PATCH_COMMAND ""
      UPDATE_COMMAND ""
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
    )
  ELSE("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
    #-------------------------------------------------------------------------------
    # If precompiled package does not exist build external project
    #-------------------------------------------------------------------------------
    ExternalProject_Add(libfbi
      DEPENDS boost
      PREFIX "${spacepart_prefix}"
      DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/spacepart
      URL ${LIBFBI_URL}/${LIBFBI_GZ}
      URL_MD5 ${LIBFBI_MD5}
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
      CMAKE_ARGS
        ${CMAKE_ARGS}
        -DBoost_DIR:PATH=${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}
        -DBoost_INCLUDE_DIR:PATH=${CFS_BINARY_DIR}/include
        -DENABLE_BENCHMARK:BOOL=OFF
        -DENABLE_EXAMPLES:BOOL=ON
        -DENABLE_TESTING:BOOL=ON
        -DENABLE_MULTITHREADING:BOOL=${USE_OPENMP}
    )
    
    IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
      #-------------------------------------------------------------------------------
      # Add custom step to zip a precompiled package to the cache.
      #-------------------------------------------------------------------------------
      ExternalProject_Add_Step(libfbi cfsdeps_zipToCache
        COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
        DEPENDEES install
        DEPENDS "${ZIPTOCACHE}"
        WORKING_DIRECTORY ${CFS_BINARY_DIR}
      )
    ENDIF()
  ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
ENDIF(USE_LIBFBI)



