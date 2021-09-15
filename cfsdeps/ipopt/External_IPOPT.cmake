#-------------------------------------------------------------------------------
# IPOPT (interior point optimizer) is a general purpose open source optimizer
#
# IPOPT needs a linear solver, there are several options but almost all have
# own license agreements to be fulfilled before they can be downloaded
# https://coin-or.github.io/Ipopt/INSTALL.html
# Therefore we restrict ourselves to use MKL pardiso, which was not possible in
# (much) older version of ipopt
#-------------------------------------------------------------------------------


set(IPOPT_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/ipopt")
set(IPOPT_SOURCE  "${IPOPT_PREFIX}/src/ipopt")
# we use this as temporary install directory to alllow packing for precompiled cfsdeps
set(IPOPT_INSTALL  "${IPOPT_PREFIX}/install")

# ipopt is build by configure, therefore no cmake args needed

SET(MIRRORS
  "http://www.coin-or.org/download/source/Ipopt/${IPOPT_ZIP}"
  "${CFS_DS_SOURCES_DIR}/ipopt/${IPOPT_ZIP}"
)

SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/ipopt/${IPOPT_ZIP}")
SET(MD5_SUM ${IPOPT_MD5})

# this handles the mirrors and cfsdepscache source stuff.  
SET(DLFN "${IPOPT_PREFIX}/ipopt-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/ipopt/ipopt-post_install.cmake.in")
SET(PI "${IPOPT_PREFIX}/ipopt-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

file(COPY "${CFS_SOURCE_DIR}/cfsdeps/ipopt/license/" DESTINATION "${CFS_BINARY_DIR}/license/ipopt" )

file(MAKE_DIRECTORY ${IPOPT_INSTALL})

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "ipopt" "${IPOPT_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${IPOPT_INSTALL}")

SET(ZIPFROMCACHE "${IPOPT_PREFIX}/ipopt-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${IPOPT_PREFIX}/ipopt-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

# Determine paths of IPOPT libraries. We have only libs include/coin
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
SET(IPOPT_LIBRARY "${LD}/libipopt.a" CACHE FILEPATH "IPOPT library.")

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
    BUILD_BYPRODUCTS ${IPOPT_LIBRARY}
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
    BUILD_IN_SOURCE 1
    
    # if on mac it complains for missing stdio.h, do export SDKROOT=$(xcrun --sdk macosx --show-sdk-path)
    # https://stackoverflow.com/questions/51761599/cannot-find-stdio-h
    
    # ipopt 3.14.2 has a but ignoring --with-lapack-lflags and searching mkl by itself. This works only when mklvars.sh is sourced
    # therefore we give LIBRARY_PATH manually, simulating sourcing mkl, what is not done on most gitlab runners
    CONFIGURE_COMMAND env "CXXFLAGS=${CFS_CXX_FLAGS}" env "LIBRARY_PATH=${MKL_LIB_DIR}" ${IPOPT_SOURCE}/configure --with-lapack-lflags="-L${MKL_LIB_DIR} -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm"  --prefix=${IPOPT_INSTALL} --exec-prefix=${IPOPT_INSTALL} --libdir=${IPOPT_INSTALL}/lib64 --enable-shared=no --enable-static --disable-java --disable-sipopt F77=${CMAKE_Fortran_COMPILER} OPT_FFLAGSS=-O3 CXX=${CMAKE_CXX_COMPILER} OPT_CXXFLAGS=-O3
    BUILD_BYPRODUCTS ${IPOPT_LIBRARY}
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
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make snopt"
SET(CFSDEPS ${CFSDEPS} ipopt)


#The ipopt depends on the metis libraries and the metis library are built by ilupack when its on.
#In this case one needs both metis and metisomp libs to be linked along with ipopt 

SET(IPOPT_INCLUDE_DIR "${CFS_BINARY_DIR}/include" CACHE FILEPATH "IPOPT library.")

MARK_AS_ADVANCED(IPOPT_LIBRARY)
MARK_AS_ADVANCED(IPOPT_INCLUDE_DIR)
