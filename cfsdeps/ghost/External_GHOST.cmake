#-------------------------------------------------------------------------------
# ghost is a hpc kernel, similar to BLAS but on a higher level which
# can use openmp, mpi and cuda.
# ghost is part of the sppexa and developed from FAU.
# ghost is a base for the phist eigenvalue solver but could also be used
# for linear algebra.
#
# We use defined revisions.
#
# ghost requires hwloc. We use cfs-hwloc as the current ghost has a wrong lib location for
# self downloaded hwloc.
#-------------------------------------------------------------------------------

set(GHOST_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ghost")
set(GHOST_SOURCE  "${GHOST_PREFIX}/src/ghost")
set(GHOST_INSTALL  "${GHOST_PREFIX}/install")

# as ghost_REVISION.zip in cfsdeps/source
if(CFS_DEPS_PRECOMPILED)
  set(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ghost/${GHOST_ZIP}")
else()
  set(LOCAL_FILE "${GHOST_SOURCE}/${GHOST_ZIP}")
endif()  


SET(MIRRORS
  "https://bitbucket.org/${GHOST_BB_USER}/${GHOST_BB_PROJECT}/get/${GHOST_ZIP}"
  "${CFS_DS_SOURCES_DIR}/ghost/${GHOST_ZIP}")

SET(MD5_SUM ${GHOST_MD5})

# ghost automatically builds for cuda if it finds it. Add this via CUDA_FOUND to the package name
find_package(CUDA QUIET)

if(CUDA_FOUND)
  set(GHOST_CUDA "cuda")
  set(CUDA_LIBS  "${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so;${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcublas.so;${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcurand.so;${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse.so;" CACHE FILEPATH "cuda libs.")
  message(STATUS "ghost: cuda found on the system" )

else()
  set(GHOST_CUDA "no-cuda")
  message(STATUS "ghost: cuda not found on the system")
endif()

SET(DLFN "${GHOST_PREFIX}/ghost-download.cmake") 
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
set(PI "${GHOST_PREFIX}/ghost-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/ghost/ghost-post_install.cmake.in" "${PI}" @ONLY) 

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "ghost_${GHOST_CUDA}" "${GHOST_REV}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
set(TMP_DIR "${GHOST_INSTALL}")

set(ZIPFROMCACHE "${GHOST_PREFIX}/ghost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${GHOST_PREFIX}/ghost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of ghost libraries.
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
set(GHOST_LIBRARY  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}ghost${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "ghost library.")
set(GHOST_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(GHOST_LIBRARY)
MARK_AS_ADVANCED(GHOST_INCLUDE_DIR)

if(DEBUG)
  set(GHOST_VERBOSITY "1") # print warnings, errors and further important information 
else()
  set(GHOST_VERBOSITY "0") # be totally silent
endif()

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${GHOST_INSTALL}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_BUILD_TYPE:STRING=Release
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  # otherwise the linking of minimal, simple, would faild due to missing -lpciaccess
  -DGHOST_BUILD_TEST:BOOL=OFF
  -DGHOST_USE_MPI:BOOL=OFF
  -DGHOST_USE_OPENMP:BOOL=${USE_OPENMP}
  -DGHOST_VERBOSITY=${GHOST_VERBOSITY}
  # we use or own built stuff, the implicit download via ghost as issues
  -DHWLOC_INCLUDE_DIR=${HWLOC_INCLUDE_DIR}
  -DHWLOC_LIBRARIES=${HWLOC_LIBRARY}
)

#-------------------------------------------------------------------------------
# The ghost external project
#-------------------------------------------------------------------------------
if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ghost
    PREFIX "${GHOST_PREFIX}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
else()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  #message("build new as it doesn't exist ${PRECOMPILED_PCKG_FILE}")
  ExternalProject_Add(ghost
    PREFIX "${GHOST_PREFIX}"
    SOURCE_DIR "${GHOST_SOURCE}"
    URL ${LOCAL_FILE}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    CMAKE_ARGS ${CMAKE_ARGS}
    INSTALL_COMMAND make -f Makefile install
    BUILD_BYPRODUCTS ${GHOST_LIBRARY})

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ghost cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${GHOST_PREFIX})

  #-------------------------------------------------------------------------------
  # Execute the stuff from ghost-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ghost post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install )
  
  if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ghost cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR})
  endif()
endif()

# Add project to global list of CFSDEPS, this allows "make ghost"
set(CFSDEPS ${CFSDEPS} ghost)
