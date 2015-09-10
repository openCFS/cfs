#-------------------------------------------------------------------------------
# SGPP (sparse nonlinear optimizer) is a commercial general purpose optimizer
# which can be used efficiently for structural optimization.
#
# The chair applied mathematics 2 (FAU) has a license for the code.
#
# You need to request the encrypted file sgopt.zip with the source 
# from the optimization group. 
# Place sgopt.zip to the cfsdepscache source/sgpp directory
# and set the password for the file to CFS_KEY_SGPP (e.g. in .cfs_platform_defaults.cmake)
#-------------------------------------------------------------------------------

set(SGPP_PREFIX  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/sgpp")
set(SGPP_SOURCE  "${SGPP_PREFIX}/src/sgpp")
# we use this as temporary install directory to alllow packing for precompiled cfsdeps
#set(SGPP_INSTALL  "${SGPP_PREFIX}/src/sgpp/${SGPP_BASE}_${SGPP_VER}")
set(SGPP_INSTALL  "${SGPP_PREFIX}/src/sgpp/${SGPP_BASE}")

# sgpp is build by configure, therefore no cmake args needed

# no patches are required

# we do not download but assume the file to magically exist there. It is encrypted!
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/sgpp/${SGPP_ZIP}")
# we don't check the md5 here

if(NOT EXISTS "${LOCAL_FILE}")
  message(FATAL_ERROR "Please provide the encypted sgpp source in ${LOCAL_FILE}") 
endif()

if(NOT CFS_KEY_SGPP)
  message(FATAL_ERROR "Key for encrypted ${SGPP_ZIP} required in CFS_KEY_SGPP. Set e.g. in your ~/.cfs_platform_defaults.cmake")
endif()    

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/sgpp/sgpp-post_install.cmake.in")
SET(PI "${SGPP_PREFIX}/sgpp-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY) 

IF(WIN32)
  SET(PRECOMPILED_PCKG_NAME "sgpp_${SGPP_VER}_${CFS_ARCH_STR}_${TOOLSET_ID}_${CMAKE_BUILD_TYPE}.zip")
ELSE(WIN32)
  SET(PRECOMPILED_PCKG_NAME "sgpp_${SGPP_VER}_${CFS_ARCH_STR}_${FC_ID}_${CMAKE_BUILD_TYPE}.zip")
ENDIF(WIN32)
SET(PRECOMPILED_PCKG_FILE "${CFS_DEPS_CACHE_DIR}/precompiled/CFSDEPS/${PRECOMPILED_PCKG_NAME}")
  
SET(PREFIX_DIR "${SGPP_PREFIX}")

SET(ZIPFROMCACHE "${SGPP_PREFIX}/sgpp-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${SGPP_PREFIX}/sgpp-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The sgpp external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(sgpp
    PREFIX "${SGPP_PREFIX}"
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
  ExternalProject_Add(sgpp
    PREFIX "${SGPP_PREFIX}"
    SOURCE_DIR "${SGPP_SOURCE}"
#    BINARY_DIR "${SGPP_SOURCE}/${SGPP_BASE}_${SGPP_VER}"
    BINARY_DIR "${SGPP_SOURCE}/${SGPP_BASE}"
    # don't dowload but simply copy w/o md5 check
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy ${LOCAL_FILE} ${SGPP_SOURCE}
    # the zip file is encrypted!
    PATCH_COMMAND unzip -q -u -P ${CFS_KEY_SGPP} ${SGPP_ZIP}
    INSTALL_COMMAND ""
    CONFIGURE_COMMAND ""
    # the libs will be created in lib/sgpp and we manually copy them to lib64/CFS_ARCH_STR
#    BUILD_COMMAND scons -j 4 -s SG_OPT=yes OMP=yes UMFPACK=yes EIGEN=yes ARMADILLO=yes GMMPP=yes
    BUILD_COMMAND scons -j 4 -s OMP=yes USE_UMFPACK=yes USE_EIGEN=yes USE_ARMADILLO=yes USE_GMMPP=no NO_UNIT_TESTS=yes
  )  
  
  #-------------------------------------------------------------------------------
  # execute the stuff from sgpp-post_install.cmake after installation
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(sgpp post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
  )
  
  IF("${CFS_DEPS_TOCACHE}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(sgpp cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_CACHE}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

# Add project to global list of CFSDEPS, this allows "make sgpp"
SET(CFSDEPS ${CFSDEPS} sgpp)

# Determine paths of SGPP libraries.
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
# Old SGPP (SGpp)
#SET(SGPP_LIBRARY "${LD}/libsgppopt.a;${LD}/libumfpack4sgpp.a;${LD}/libsuitesparseconfig4sgpp.a;${LD}/libarmadillo4sgpp.a;${LD}/libsgppbase.a;${LD}/libsgpppde.a;${LD}/libsgppsolver.a" CACHE FILEPATH "SGPP library.")
# New SGPP (sgopt)
SET(SGPP_LIBRARY 
  ${LD}/libsgppoptimization.a;
  ${LD}/libsgppbase.a;
  libarmadillo.so;
  ${UMFPACK_LIBRARY};
  ${ARPACK_LIBRARY}
  CACHE FILEPATH "SGPP library."
  )

MESSAGE("----SGPP_LIBRARY----")
MESSAGE("${SGPP_LIBRARY}")

MARK_AS_ADVANCED(SGPP_LIBRARY)
