#-------------------------------------------------------------------------------
# ARnoldi PACKage -  ARPACK is a collection of  Fortran77 subroutines designed
# to solve large scale eigenvalue problems. 
# We use the maintained arpack-ng variant
#
# Project Homepage
# https://github.com/opencollab/arpack-ng
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to ARPACK sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(ARPACK_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/arpack")
set(ARPACK_source  "${ARPACK_prefix}/src/arpack")
set(ARPACK_install  "${CMAKE_CURRENT_BINARY_DIR}")

SET(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${ARPACK_install}
  -DCMAKE_INSTALL_LIBDIR:PATH=${ARPACK_install}/${LIB_SUFFIX}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  # for our 3.7.0 patch only
  -DSYSTEM_BLAS:BOOL=OFF
  -DSYSTEM_LAPACK:BOOL=OFF
  -DTESTS:BOOL=OFF
)


IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/arpack/arpack-patch.cmake.in")
SET(PFN "${ARPACK_prefix}/arpack-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 


#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "https://github.com/opencollab/arpack-ng/archive/${ARPACK_GZ}"
  "${ARPACK_URL}/${ARPACK_GZ}"
)

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/arpack/${ARPACK_GZ}")
SET(MD5_SUM ${ARPACK_MD5})

SET(DLFN "${ARPACK_prefix}/arpack-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/arpack/arpack-post_install.cmake.in")
SET(PI "${ARPACK_prefix}/arpack-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/arpack/license/" DESTINATION "${CFS_BINARY_DIR}/license/arpack" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "arpack" "${ARPACK_VER}") 

# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${ARPACK_prefix}")

SET(ZIPFROMCACHE "${ARPACK_prefix}/arpack-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${ARPACK_prefix}/arpack-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of ARPACK libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
IF(WIN32)
  SET(ARPACK_LIBRARY "${CFS_BINARY_DIR}/${LIB_SUFFIX}/arpack.lib" CACHE FILEPATH "ARPACK library.")
ELSE(WIN32)
  SET(ARPACK_LIBRARY "${LD}/libarpack.a" CACHE FILEPATH "ARPACK library.")
ENDIF(WIN32)
MARK_AS_ADVANCED(ARPACK_LIBRARY)

#-------------------------------------------------------------------------------
# The ARPACK external project
#-------------------------------------------------------------------------------

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(arpack
    PREFIX "${ARPACK_prefix}"
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
  ExternalProject_Add(arpack
    PREFIX "${ARPACK_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${ARPACK_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS ${CMAKE_ARGS}
    BUILD_BYPRODUCTS ${ARPACK_LIBRARY}
    WORKING_DIRECTORY ${ARPACK_prefix}
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(arpack cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${ARPACK_prefix}
  )
  
  ExternalProject_Add_Step(arpack post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(arpack cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS ${CFSDEPS} arpack)

