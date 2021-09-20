#-------------------------------------------------------------------------------
# Library of Iterative Solvers
#
# Project Homepage
# http://www.ssisc.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to lis sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(lis_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/lis")
set(lis_source  "${lis_prefix}/src/lis")
set(lis_install  "${lis_prefix}/install")

IF(UNIX)
  set(LIS_CONFIG_OPTIONS
     "--prefix=${lis_install}"
     "--libdir=${lis_install}/lib64"
     "--includedir=${lis_install}/include"
     "--enable-saamg"
     # "--enable-quad" #          enable quadruple precision operations [[default=no]]
     "--enable-sse2" #          use Intel Streaming SIMD Extensions [[default=yes]]
     "--enable-fma"  #          use fused multiply add [[default=no]]
     #"--enable-complex"    #    enable complex scalar support [[default=no]]
     "--enable-static" # [=PKGS]  build static libraries [default=yes]
     #  --enable-fast-install[=PKGS]  optimize for fast installation [default=yes]
  )
  # and C compiler info, note that --enable-fortran only gives the interfce we don' need
  list(APPEND LIS_CONFIG_OPTIONS "CC=${CMAKE_C_COMPILER}" "CFLAGS=${CFSDEPS_C_FLAGS}")

  if(USE_OPENMP)
    list(APPEND LIS_CONFIG_OPTIONS "--enable-omp")
  endif()

ELSE()
  set(lis_install  "${CMAKE_CURRENT_BINARY_DIR}")

  SET(CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${lis_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CFSDEPS_C_FLAGS}
    -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_FLAGS:STRING=${CFSDEPS_Fortran_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  )
  IF(CMAKE_TOOLCHAIN_FILE)
    LIST(APPEND CMAKE_ARGS
      -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
    )
  ENDIF()
  IF(OPENMP_FOUND)
    LIST(APPEND CMAKE_ARGS -DLIS_ENABLE_OPENMP:BOOL=ON)
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Windows requires a patch 
#-------------------------------------------------------------------------------
IF(WIN32)
  SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lis/lis-patch.cmake.in")
  SET(PFN "${lis_prefix}/lis-patch.cmake")
  CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY)
#  SET(PF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lis/lis_install.cmake.in")
#  SET(PF "${lis_prefix}/lis_install.cmake")
#  CONFIGURE_FILE("${PF_TEMPL}" "${PF}" @ONLY) 
ENDIF()

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.ssisc.org/lis/dl/${LIS_ZIP}"
  "${CFS_DS_SOURCES_DIR}/lis/${LIS_ZIP}"
  "${LIS_URL}/${LIS_ZIP}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/lis/${LIS_ZIP}")
SET(MD5_SUM ${LIS_MD5})

SET(DLFN "${lis_prefix}/lis-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/lis/license/" DESTINATION "${CFS_BINARY_DIR}/license/lis" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "lis" "${LIS_VER}")
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
IF(WIN32)
  SET(TMP_DIR "${lis_prefix}")
ELSE()
  SET(TMP_DIR "${lis_install}")
ENDIF()

SET(ZIPFROMCACHE "${lis_prefix}/lis-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${lis_prefix}/lis-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/lis/lis-post_install.cmake.in")
SET(PI "${lis_prefix}/lis-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

#-----------------------------------------------------------------------------
# Determine paths of LIS libraries.
#-----------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
SET(LIS_LIBRARY
  "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}lis${CMAKE_STATIC_LIBRARY_SUFFIX}"
  CACHE FILEPATH "LIS library.")
MARK_AS_ADVANCED(LIS_LIBRARY)

#-------------------------------------------------------------------------------
# The lis external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(lis
    PREFIX "${lis_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${LIS_LIBRARY}
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  IF(UNIX)
    find_program(LIS_MAKE_PROGRAM make)
    ExternalProject_Add(lis
      PREFIX "${lis_prefix}"
      SOURCE_DIR "${lis_source}"
      URL ${LOCAL_FILE}
      URL_MD5 ${LIS_MD5}
      PATCH_COMMAND ""
      CONFIGURE_COMMAND "${lis_source}/configure" ${LIS_CONFIG_OPTIONS}
      BUILD_COMMAND ${LIS_MAKE_PROGRAM}
      INSTALL_COMMAND ${LIS_MAKE_PROGRAM} install
      BUILD_BYPRODUCTS ${LIS_LIBRARY}
    )
  ELSE()
    ExternalProject_Add(lis
      PREFIX "${lis_prefix}"
      SOURCE_DIR "${lis_source}"
      URL ${LOCAL_FILE}
      URL_MD5 ${LIS_MD5}
      PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
      CMAKE_ARGS
        ${CMAKE_ARGS}
        -DLIS_BUILD_TEST:BOOL=ON
        -DLIS_ENABLE_FORTRAN:BOOL=ON
        -DLIS_ENABLE_SAAMG:BOOL=ON
      BUILD_BYPRODUCTS ${LIS_LIBRARY}
    )
  ENDIF()
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(lis cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${lis_prefix}
  )

  #-------------------------------------------------------------------------------
  # Execute the stuff from lis-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  IF(UNIX)
    ExternalProject_Add_Step(lis post_install
      COMMAND ${CMAKE_COMMAND} -P "${PI}"
      WORKING_DIRECTORY ${lis_prefix}
      DEPENDEES install
    )
  ENDIF()
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(lis cfsdeps_zipToCache
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
SET(CFSDEPS ${CFSDEPS} lis) 

SET(LIS_INCLUDE_DIR "${CFS_BINARY_DIR}/include")
MARK_AS_ADVANCED(LIS_INCLUDE_DIR)
