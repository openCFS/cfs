#-------------------------------------------------------------------------------
# Set paths to lapack sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------

# Note that the proper name would be actually NETLIB, as the blas and lapack from NETLIB are built!

set(lapack_prefix "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/lapack")
set(lapack_source  "${lapack_prefix}/src/lapack")
set(lapack_install  "${lapack_prefix}/install")
file(MAKE_DIRECTORY ${lapack_install})

set(CMAKE_ARGS
  -DCMAKE_INSTALL_PREFIX:PATH=${lapack_install}
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
  -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
  -DBUILD_TESTING:BOOL=OFF
)

if(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE} )
endif()



#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.netlib.org/lapack/${LAPACK_GZ}"
  "https://github.com/Reference-LAPACK/lapack/archive/${LAPACK_GZ}"
  "${CFS_DS_SOURCES_DIR}/lapack/${LAPACK_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/lapack/${LAPACK_GZ}")
SET(MD5_SUM ${LAPACK_MD5})

SET(DLFN "${lapack_prefix}/lapack-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# see below for what we need the post_install (pre-zip)
SET(PI "${lapack_prefix}/lapack-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/lapack/lapack-post_install.cmake.in" "${PI}" @ONLY) 

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/lapack/license/" DESTINATION "${CFS_BINARY_DIR}/license/lapack" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "lapack" "${LAPACK_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${lapack_install}")

SET(ZIPFROMCACHE "${lapack_prefix}/lapack-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${lapack_prefix}/lapack-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The lapack external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(lapack
    PREFIX "${lapack_prefix}"
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
  ExternalProject_Add(lapack
    PREFIX "${lapack_prefix}"
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/lapack
    URL ${LOCAL_FILE}
    URL_MD5 ${LAPACK_MD5}
    CMAKE_ARGS
      ${CMAKE_ARGS}
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(lapack cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${lapack_prefix}
  )
  
  # lapack puts is stuff to lib. copy what we need and remove the rest  
  ExternalProject_Add_Step(lapack post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(lapack cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  lapack
)

#-------------------------------------------------------------------------------
# Determine paths of LAPACK libraries.
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(NETLIB_BLAS_LIBRARY
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_LAPACK_LIBRARY
    ${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX}
    CACHE FILEPATH "Netlib LAPACK library.")
ELSE()
  SET(BLAS_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/libblas.a")
  SET(LAPACK_LIB "${CFS_BINARY_DIR}/${LIB_SUFFIX}/liblapack.a")

  IF(CFS_FORTRAN_COMPILER_NAME STREQUAL "GNU")
    LIST(APPEND BLAS_LIB ${GFORTRAN_LIBRARY})
  ENDIF()

  SET(NETLIB_BLAS_LIBRARY ${BLAS_LIB} CACHE FILEPATH "Netlib BLAS library.")
  SET(NETLIB_LAPACK_LIBRARY ${LAPACK_LIB} CACHE FILEPATH "Netlib LAPACK library.")
ENDIF()

#-------------------------------------------------------------------------------
# Mark paths of LAPACK libraries as advanced.
#-------------------------------------------------------------------------------
MARK_AS_ADVANCED(NETLIB_BLAS_LIBRARY)
MARK_AS_ADVANCED(NETLIB_LAPACK_LIBRARY)

#-------------------------------------------------------------------------------
# Set LAPACK_LIBRARY according to configuration
#-------------------------------------------------------------------------------

# see also External_OpenBLAS and FindIntelMKL
# e.g. BLAS_LIBRARY=...lib64/OPENSUSE_TUMBLEWEED_X86_64/libblas.a;/usr/lib64/gcc/x86_64-suse-linux/7/libgfortran.so
# e.g. LAPACK_LIBRARY=...lib64/OPENSUSE_TUMBLEWEED_X86_64/liblapack.a
if(USE_BLAS_LAPACK STREQUAL "NETLIB")
  set(BLAS_LIBRARY "${NETLIB_BLAS_LIBRARY}")
  set(LAPACK_LIBRARY "${NETLIB_LAPACK_LIBRARY}")
endif()
