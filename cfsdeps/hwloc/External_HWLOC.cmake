#-------------------------------------------------------------------------------
# hwloc gives information about the machine strucuture like cpu, memory 
# controllers, ... . 
# PETSc and ghost make use of hwloc, therefore we offer it here for a common base 
#-------------------------------------------------------------------------------

set(HWLOC_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/hwloc")
set(HWLOC_SOURCE  "${HWLOC_PREFIX}/src/hwloc")
set(HWLOC_INSTALL  "${HWLOC_PREFIX}/install")

# hwloc is build by configure, therefore no cmake args needed
# no patches are required

SET(MIRRORS
  "https://www.open-mpi.org/software/hwloc/v1.11/downloads/${HWLOC_TGZ}"
  "https://fossies.org/linux/misc/${HWLOC_TGZ}"
  "${CFS_DS_SOURCES_DIR}/hwloc/${HWLOC_TGZ}")

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/hwloc/${HWLOC_TGZ}")
SET(MD5_SUM ${HWLOC_MD5})

SET(DLFN "${HWLOC_PREFIX}/hwloc-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/hwloc/hwloc-post_install.cmake.in")
SET(PI "${HWLOC_PREFIX}/hwloc-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "hwloc" "${HWLOC_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${HWLOC_INSTALL}")

SET(ZIPFROMCACHE "${HWLOC_PREFIX}/hwloc-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${HWLOC_PREFIX}/hwloc-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of hwloc libraries.

# we need runtime libs, which are quite different for the systems:
# openSUSE tumbleweed: -lpciaccess;-ludev or strangely /usr/lib64/libpciaccess.so.0;/usr/lib64/libudev.so.1
# leo (ubuntu 14.04): -lnuma;-lOpenCL
# woody (ubuntu 16.04): -lpciaccess;-lnuma
if(HWLOC_SYSTEM_LIBS_DEFAULT)
  set(HWLOC_TMP ${HWLOC_SYSTEM_LIBS_DEFAULT})
else()
  set(HWLOC_TMP "-lpciaccess;-ludev")
endif()  

SET(HWLOC_SYSTEM_LIBS ${HWLOC_TMP} CACHE STRING "Set system libs for hwloc, e.g. '-lpciaccess;-ludev','-lpciaccess;-lnuma','-lnuma;-lOpenCL'")

SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
SET(HWLOC_LIBRARY  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}hwloc.a;${HWLOC_SYSTEM_LIBS}" CACHE FILEPATH "hwloc library, see HWLOC_SYSTEM_LIBS" FORCE)

SET(HWLOC_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

MARK_AS_ADVANCED(HWLOC_LIBRARY)
MARK_AS_ADVANCED(HWLOC_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# The hwloc external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(hwloc
    PREFIX "${HWLOC_PREFIX}"
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
  ExternalProject_Add(hwloc
    PREFIX "${HWLOC_PREFIX}"
    SOURCE_DIR "${HWLOC_SOURCE}"
    URL ${LOCAL_FILE}
    URL_MD5 ${HWLOC_MD5}
    BUILD_IN_SOURCE 1
    PATCH_COMMAND ""
    # disable some stuff we probably don't need  
    CONFIGURE_COMMAND ${HWLOC_SOURCE}/configure --prefix=${HWLOC_INSTALL} --libdir=${HWLOC_INSTALL}/${LIB_SUFFIX} --enable-static --disable-libxml2 --disable-cairo
    INSTALL_COMMAND make -f Makefile install
    BUILD_BYPRODUCTS ${HWLOC_LIBRARY}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(hwloc cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${HWLOC_PREFIX}
  )

  #-------------------------------------------------------------------------------
  # Execute the stuff from hwloc-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(hwloc post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(hwloc cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make hwloc"
SET(CFSDEPS ${CFSDEPS} hwloc)
