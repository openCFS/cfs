#-------------------------------------------------------------------------------
# ghost is a hpc kernel, similar to BLAS but on a higher level which
# can use openmp, mpi and cuda.
# ghost is part of the sppexa and developed from FAU.
# ghost is a base for the phist eigenvalue solver but could also be used
# for linear algebra.
# https://github.com/RRZE-HPC/GHOST which is a mirror of https://bitbucket.org/essex/ghost
#
# We use defined revisions.
#
# ghost requires hwloc. We use cfs-hwloc as the current ghost has a wrong lib location for
# self downloaded hwloc.
#-------------------------------------------------------------------------------

set(GHOST_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ghost")
set(GHOST_SOURCE  "${GHOST_PREFIX}/src/ghost")
set(GHOST_INSTALL  "${GHOST_PREFIX}/install")

# we download via svn co https://github.com/RRZE-HPC/GHOST/trunk@r<REVSION> and pack it
# as ghost_REVISION.zip in cfsdeps/source
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ghost/ghost_r${GHOST_REV}.zip")
SET(MD5_SUM ${GHOST_MD5})

# if the cachhed source file does not exist, do svn co ...
SET(DLFN "${GHOST_PREFIX}/ghost-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/ghost/ghost-download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
SET(PI "${GHOST_PREFIX}/ghost-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/ghost/ghost-post_install.cmake.in" "${PI}" @ONLY) 

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "ghost" "${GHOST_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${GHOST_INSTALL}")

SET(ZIPFROMCACHE "${GHOST_PREFIX}/ghost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${GHOST_PREFIX}/ghost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of ghost libraries.
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(GHOST_LIBRARY  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}ghost${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "ghost library.")
SET(GHOST_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(GHOST_LIBRARY)
MARK_AS_ADVANCED(GHOST_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The ghost external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
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
    INSTALL_COMMAND ""
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ghost
    PREFIX "${GHOST_PREFIX}"
    SOURCE_DIR "${GHOST_SOURCE}"
    URL ${LOCAL_FILE}
    URL_MD5 ${GHOST_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    # disable some stuff we probably don't need  
    CONFIGURE_COMMAND ${GHOST_SOURCE}/configure --prefix=${GHOST_INSTALL} --libdir=${GHOST_INSTALL}/${LIB_SUFFIX}/${CFS_ARCH_STR} --enable-static --disable-libxml2 --disable-cairo
    INSTALL_COMMAND ${CONFIGURE_MAKE_PROGRAM} -f Makefile install
    BUILD_BYPRODUCTS ${GHOST_LIBRARY}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ghost cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${GHOST_PREFIX}
  )

  #-------------------------------------------------------------------------------
  # Execute the stuff from ghost-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(ghost post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ghost cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make ghost"
SET(CFSDEPS ${CFSDEPS} ghost)
