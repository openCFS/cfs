#-------------------------------------------------------------------------------
# phist by Jonas This (DLR) a HPC eigenvalue solver from the SPPEXA project.
# The FAU is part of the SPPEXA and provides the kernel ghost for phist.
# We use phist here without MPI and are fixed for ghost. 
# 
# phist needs a ghost-config.cmake. As there are full path names this interferes
# with the precompiled packages. Therefore we construct a own ghost-config.cmake here!
#-------------------------------------------------------------------------------

set(PHIST_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/phist")
set(PHIST_SOURCE  "${PHIST_PREFIX}/src/phist")
set(PHIST_INSTALL  "${PHIST_PREFIX}/install")

# as phist_REVISION.zip in cfsdeps/source
if(CFS_DEPS_PRECOMPILED)
  set(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/phist/${PHIST_ZIP}")
else()
  set(LOCAL_FILE "${PHIST_SOURCE/${PHIST_ZIP}")
endif()  

SET(MIRRORS
  "https://bitbucket.org/${PHIST_BB_USER}/${PHIST_BB_PROJECT}/get/${PHIST_ZIP}"
  "${CFS_DS_SOURCES_DIR}/phist/${PHIST_ZIP}")

SET(MD5_SUM ${PHIST_MD5})

SET(DLFN "${PHIST_PREFIX}/phist-download.cmake") 
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# configure the ghost-config.make file and copy it in phist-post_download.cmake
set(GC "${PHIST_PREFIX}/ghost-config.cmake")
configure_file("${CFS_SOURCE_DIR}/cfsdeps/phist/ghost-config.cmake.in" "${GC}" @ONLY)

# After download we fake ghost-config.cmake.
set(PD "${PHIST_PREFIX}/phist-post_download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/phist/phist-post_download.cmake.in" "${PD}" @ONLY) 

# After the installation we copy to cfs
set(PI "${PHIST_PREFIX}/phist-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/phist/phist-post_install.cmake.in" "${PI}" @ONLY) 

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "phist" "${PHIST_REV}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
set(TMP_DIR "${PHIST_INSTALL}")

set(ZIPFROMCACHE "${PHIST_PREFIX}/phist-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${PHIST_PREFIX}/phist-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of phist libraries.
set(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
set(PHIST_LIBRARY  "${LD}/libphist_kernels_ghost.a;${LD}/libphist_solvers.a;${LD}/libphist_core.a;${LD}/libphist_tools.a;${LD}/libphist_precon.a" CACHE FILEPATH "phist libraries.")
set(PHIST_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(PHIST_LIBRARY)
MARK_AS_ADVANCED(PHIST_INCLUDE_DIR)

# see phist_config.h.in
if(DEBUG)
  set(PHIST_OUTLEV "3") # verbose info
else()
  set(PHIST_OUTLEV "0") # only print errors
endif()

if(NOT USE_OPENMP)
  message(FATAL_ERROR "phist eigenvalue solver requires OpenMP")
endif()


SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${PHIST_INSTALL}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_MAKE_PROGRAM:FILEPATH=${CMAKE_MAKE_PROGRAM}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  #-DPYTHON_EXECUTABLE:FILEPATH=python2 # not yet python3 compatible
  # We do not want to see warning messages from external projects
  -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
  -DCMAKE_CXX_FLAGS:STRING=${CFSDEPS_CXX_FLAGS}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  # phist settings 
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DPHIST_OUTLEV=${PHIST_OUTLEV}
  -DPHIST_ENABLE_MPI:BOOL=OFF
  # will not compile w/o openmp and concurrently w/o mpi
  -DPHIST_ENABLE_OPENMP:BOOL=ON
  -DPHIST_ENABLE_COMPLEX:BOOL=ON
  -DPHIST_KERNEL_LIB:STRING=ghost
  -DGHOST_DIR:STRING=${CMAKE_CURRENT_BINARY_DIR}/share/ghost
)

#-------------------------------------------------------------------------------
# The phist external project
#-------------------------------------------------------------------------------
if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(phist
    PREFIX "${PHIST_PREFIX}"
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
  ExternalProject_Add(phist
    PREFIX "${PHIST_PREFIX}"
    SOURCE_DIR "${PHIST_SOURCE}"
    URL ${LOCAL_FILE}
    BUILD_IN_SOURCE 1
    # somehow I cannot make the post_downlowad step work, hence we do it here.
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PD}"
    CMAKE_ARGS ${CMAKE_ARGS}
    BUILD_COMMAND make libs
    INSTALL_COMMAND make install
    BUILD_BYPRODUCTS ${PHIST_LIBRARY})

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(phist cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${PHIST_PREFIX})

  ExternalProject_Add_Step(phist post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install )
  
  if("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(phist cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR})
  endif()
endif()

# Add project to global list of CFSDEPS, this allows "make phist"
set(CFSDEPS ${CFSDEPS} phist)
