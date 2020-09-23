#-------------------------------------------------------------------------------
# FEAST Eigenvalue Solver
#
# Project Homepage
# http://www.feast-solver.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to FEAST sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
SET(feast_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/feast")
SET(feast_source  "${feast_prefix}/src/feast")
SET(feast_install "${feast_prefix}/install")
SET(feast_include "${feast_source}/${FEAST_VER}/include")

# determine which feast library to copy over in the install step
IF(UNIX)
  IF(USE_FEAST_COMMUNITY_PRECOMPILED)
    SET(FEAST_LIB_DIR "${feast_source}/${FEAST_VER}/lib/x64")
  ELSE()
    SET(FEAST_LIB_DIR "${feast_source}/${FEAST_VER}/lib/${CFS_ARCH_STR}")
  ENDIF()
ELSE()
  SET(FEAST_LIB_DIR "${feast_install}/lib64")
ENDIF()

#-----------------------------------------------------------------------------
# Determine paths of FEAST includes.
#-----------------------------------------------------------------------------
SET(FEAST_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include/feast" CACHE PATH "include path for FEAST")
MARK_AS_ADVANCED(FEAST_INCLUDE_DIR)

#-------------------------------------------------------------------------------
# Configure FEAST by copying the config file.
#-------------------------------------------------------------------------------
IF(UNIX)
  SET(CONF_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/feast/make.inc.in")
  SET(CONF "${feast_prefix}/make.inc"	)
  CONFIGURE_FILE("${CONF_TEMPL}" "${CONF}" @ONLY)
ELSE()
  SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/feast/feast-patch.cmake.in")
  SET(PFN "${feast_prefix}/feast-patch.cmake")
  CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY)
ENDIF()

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://www.ecs.umass.edu/~polizzi/feast/m3-0/feast_3.0.tgz"
  "${FEAST_URL}/${FEAST_GZ}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/feast/${FEAST_GZ}")
SET(MD5_SUM ${FEAST_MD5})

SET(DLFN "${feast_prefix}/feast-download.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

#-------------------------------------------------------------------------------
# After the installation we copy to cfs
#-------------------------------------------------------------------------------
SET(PI_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/feast/feast-post_install.cmake.in")
SET(PI "${feast_prefix}/feast-post_install.cmake")
CONFIGURE_FILE("${PI_TEMPL}" "${PI}" @ONLY)

PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "feast" "${FEAST_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${feast_prefix}")

SET(ZIPFROMCACHE "${feast_prefix}/feast-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${feast_prefix}/feast-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-----------------------------------------------------------------------------
# Determine paths of FEAST libraries.
#-----------------------------------------------------------------------------
IF(UNIX)
  SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
  SET(FEAST_LIBS "feast" "feast_sparse")
  foreach(LIB IN LISTS FEAST_LIBS)
    LIST(APPEND LIBS "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}${LIB}${CMAKE_STATIC_LIBRARY_SUFFIX}")
  endforeach()
  SET(FEAST_LIBRARY ${LIBS} CACHE FILEPATH "FEAST library.")
ELSE()
  SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")
  SET(FEAST_LIBRARY "${LD}/${CMAKE_STATIC_LIBRARY_PREFIX}feast${CMAKE_STATIC_LIBRARY_SUFFIX}" CACHE FILEPATH "FEAST library.")
ENDIF()
MARK_AS_ADVANCED(FEAST_LIBRARY)

IF(WIN32)
  IF(USE_OPENMP)
    SET(FEAST_FORTRAN_FLAGS "${FEAST_FORTRAN_FLAGS} /Qopenmp")
  ENDIF(USE_OPENMP)

  SET(CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX:PATH=${feast_install}
    -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
    -DCMAKE_Fortran_COMPILER:FILEPATH=${CMAKE_Fortran_COMPILER}
    -DCMAKE_Fortran_FLAGS:STRING=${FEAST_FORTRAN_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=Release
    -DLIB_SUFFIX:STRING=${LIB_SUFFIX}
  )

  IF(CMAKE_TOOLCHAIN_FILE)
    LIST(APPEND CMAKE_ARGS
      -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
    )
  ENDIF()
ENDIF()
#-------------------------------------------------------------------------------
# The FEAST external project
#-------------------------------------------------------------------------------

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(feast
    PREFIX "${feast_prefix}"
    DOWNLOAD_COMMAND ${CMAKE_COMMAND} -P "${ZIPFROMCACHE}"
    PATCH_COMMAND ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${FEAST_LIBRARY}
  )
ELSE("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  IF("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    # FEAST does not use CMake but make only: We cannot use the CMake Ninja-generator,
    # we use make instead to build FEAST
    find_program(FEAST_MAKE_PROGRAM make)
  ELSE("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    SET(FEAST_MAKE_PROGRAM ${CMAKE_MAKE_PROGRAM} CACHE FILEPATH "program to build FEAST")
  ENDIF("${CMAKE_GENERATOR}" STREQUAL "Ninja")
  MARK_AS_ADVANCED(FEAST_MAKE_PROGRAM)
  # here we should decide if we need to build feast, or only copy the pre-compiled libs a below
  IF(USE_FEAST_COMMUNITY_PRECOMPILED)
  # only copy over
  ExternalProject_Add(feast
    PREFIX "${feast_prefix}"
    SOURCE_DIR "${feast_source}"
    URL ${LOCAL_FILE}
    URL_MD5 "${FEAST_MD5}"
    BINARY_DIR "${feast_source}/${FEAST_VER}/src"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${FEAST_LIBRARY}
  )
  ELSE() # not USE_FEAST_COMMUNITY_PRECOMPILED
    # compile feast
    IF(UNIX)
      ExternalProject_Add(feast
        PREFIX "${feast_prefix}"
        SOURCE_DIR "${feast_source}"
        URL ${LOCAL_FILE}
        URL_MD5 "${FEAST_MD5}"
        BINARY_DIR "${feast_source}/${FEAST_VER}/src"
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy "${CONF}" "${feast_source}/${FEAST_VER}/src" # copy over file
        BUILD_COMMAND ${FEAST_MAKE_PROGRAM} "ARCH=${CFS_ARCH_STR}" "LIB=feast" "all"
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS ${FEAST_LIBRARY}
    )
    ELSE()
      ExternalProject_Add(feast
        PREFIX "${feast_prefix}"
        DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/feast
        URL ${FEAST_URL}/${FEAST_GZ}
        URL_MD5 ${FEAST_MD5}
        PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
	SOURCE_DIR "${feast_source}"
        CMAKE_ARGS ${CMAKE_ARGS}
        INSTALL_DIR ${feast_install}
        LOG_CONFIGURE 1
        LOG_BUILD 1
        LOG_INSTALL 1
      )
    ENDIF()
  ENDIF() # switch of USE_FEAST_COMMUNITY_PRECOMPILED

  #-------------------------------------------------------------------------------
  # Add custom patch step, needed for windows
  #-------------------------------------------------------------------------------
  #  ExternalProject_Add_Step(feast winpatch
  #  COMMAND ${CMAKE_COMMAND} -P "${PFN}"
  #  DEPENDERS build
  #  DEPENDEES download
  #  DEPENDS "${PFN}"
  #  WORKING_DIRECTORY ${BOOST_source}
  #)

  # post install for precompiled
  ExternalProject_Add_Step(feast post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install
    DEPENDS "${PI}"
  )
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(feast cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${feast_prefix}
  )

  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(feast cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES post_install
      DEPENDS "${ZIPTOCACHE}"
      WORKING_DIRECTORY ${CFS_BINARY_DIR}
    )
  ENDIF()
ENDIF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")

#-------------------------------------------------------------------------------
# Add project to global list of CFSDEPS
#-------------------------------------------------------------------------------
SET(CFSDEPS
  ${CFSDEPS}
  feast
)
