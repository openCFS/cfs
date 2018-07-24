#-------------------------------------------------------------------------------
# IPOPT (interior point optimizer) is a general purpose open source optimizer
#
# IPOPT needs a solver. This can be HSL, which can be requested for academic purpose
# for free or MKL-Pardiso.
# * HSL is provied encrypted via ipopt_hsl.zip and requires CFS_KEY_IPOPT_HSL
# * MKL pardiso is only recently possible for IPOPT and not used here. You might 
#   change this file to use MKL pardiso!
#-------------------------------------------------------------------------------


set(IPOPT_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ipopt")
set(IPOPT_SOURCE  "${IPOPT_PREFIX}/src/ipopt")
# we use this as temporary install directory to alllow packing for precompiled cfsdeps
set(IPOPT_INSTALL  "${IPOPT_PREFIX}/install")

# ipopt is build by configure, therefore no cmake args needed

SET(MIRRORS
  "http://www.coin-or.org/download/source/Ipopt/${IPOPT_TGZ}"
  "${CFS_DS_SOURCES_DIR}/ipopt/${IPOPT_TGZ}"
)

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ipopt/${IPOPT_TGZ}")
SET(MD5_SUM ${IPOPT_MD5})

# the encrypted ipopt_hsl.zip is so small and for everyone available on the web if you register, we provide it
# in the cfs/cfsdeps/ipopt itself. You just need the key.
if(NOT CFS_KEY_IPOPT_HSL)
  message(FATAL_ERROR "Key for encrypted ipopt_hsl.zip required in CFS_KEY_IPOPT_HSL. E.g. set in your ~/.cfs_platform_defaults.cmake")
endif() 
# The patch triggers downloading netlib blas and lapack such that ipopt compiles with blas support.
# However netlib blas will be removed after installation and the cfs stuff will be used
# Furthermore ipopt_hsl is extracted 
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ipopt/ipopt-patch.cmake.in")
SET(PFN "${IPOPT_PREFIX}/ipopt-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

# this handles the mirrors and cfsdepscache source stuff.  
SET(DLFN "${IPOPT_PREFIX}/ipopt-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ipopt/ipopt-post_install.cmake.in")
SET(PI "${IPOPT_PREFIX}/ipopt-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

file(MAKE_DIRECTORY ${IPOPT_INSTALL})

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "ipopt" "${IPOPT_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${IPOPT_INSTALL}")

SET(ZIPFROMCACHE "${IPOPT_PREFIX}/ipopt-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${IPOPT_PREFIX}/ipopt-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The IPOPT external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(ipopt
    PREFIX "${IPOPT_PREFIX}"
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
  ExternalProject_Add(ipopt
    PREFIX "${IPOPT_PREFIX}"
    SOURCE_DIR "${IPOPT_SOURCE}"
    # the step cfsdeps_download ensures the file will be found in cfsdepscache
    URL ${LOCAL_FILE}
    URL_MD5 ${IPOPT_MD5}
    # note that when unzipping the source, the Ipopt-3.11.9 directory is omitted
    BUILD_IN_SOURCE 1
    # handle mirror
    # Fill ThirdParty content (blas/ lapack/ HSLold)
    PATCH_COMMAND  ${CMAKE_COMMAND} -P "${PFN}"
    # let it install to the temporay directory where we can remove libcoinblas and libcoinlapack and prepare to copy to precompiled cfsdeps
    CONFIGURE_COMMAND env "CFLAGS=${CFS_C_FLAGS}" "CXXFLAGS=${CFS_CXX_FLAGS}" ${IPOPT_SOURCE}/configure --prefix=${IPOPT_INSTALL} --libdir=${IPOPT_INSTALL}/lib64/${CFS_ARCH_STR} --disable-shared --disable-linear-solver-loader --with-metis-lib=${METIS_LIBRARY} --with-metis-incdir=${CMAKE_CURRENT_BINARY_DIR}/include --disable-pkg-config F77=${CMAKE_Fortran_COMPILER} OPT_FFLAGSS=-O3 CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} OPT_CXXFLAGS=-O3
  )
  
  ExternalProject_Add_Step(ipopt cfsdeps_download 
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}" 
    DEPENDERS download 
    DEPENDS "${DLFN}" 
    WORKING_DIRECTORY ${ipopt_prefix} 
  ) 

  ExternalProject_Add_Step(ipopt post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(ipopt cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make snopt"
SET(CFSDEPS ${CFSDEPS} ipopt)

# Determine paths of IPOPT libraries. We have only libs include/coin
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(IPOPT_LIBRARY "${LD}/libipopt.a;${LD}/libcoinhsl.a" CACHE FILEPATH "IPOPT library.")
SET(IPOPT_INCLUDE_DIR "${CFS_BINARY_DIR}/include" CACHE FILEPATH "IPOPT library.")

MARK_AS_ADVANCED(IPOPT_LIBRARY)
MARK_AS_ADVANCED(IPOPT_INCLUDE_DIR)
