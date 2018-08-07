#-------------------------------------------------------------------------------
# Boost provides free peer-reviewed portable C++ source libraries.
# Needed by Botan and CGAL
#
# Project Homepage
# http://www.boost.org
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set prefix path and path to Boost sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(BOOST_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/boost")
set(BOOST_source  "${BOOST_prefix}/src/boost")
set(BOOST_install "${BOOST_prefix}/install")

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
#SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/boost/boost-patch.cmake.in")
#SET(PFN "${ARPACK_prefix}/arpack-patch.cmake")
#CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/boost/${BOOST_GZ}")

# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${BOOST_prefix}")

SET(ZIPFROMCACHE "${BOOST_prefix}/boost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${BOOST_prefix}/boost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

SET(MIRRORS
  "https://dl.bintray.com/boostorg/release/${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}.0/source/${BOOST_GZ}"
  "${CFS_DS_SOURCES_DIR}/boost/${BOOST_ZIP}"
)

SET(MD5_SUM ${BOOST_MD5})

SET(DLFN "${BOOST_prefix}/boost-download.cmake") 
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
set(PI "${BOOST_prefix}/boost-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/boost/boost-post_install.cmake.in" "${PI}" @ONLY) 

PRECOMPILED_ZIP_NOBUILD(PRECOMPILED_PCKG_FILE "boost" "${BOOST_VER}")
 
set(TMP_DIR "${BOOST_install}")

set(ZIPFROMCACHE "${BOOST_prefix}/boost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${BOOST_prefix}/boost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of BOOST libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")

SET(Boost_LIBRARY_DIR "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}" CACHE PATH "Boost library dir.")

#SET(BOOST_EXTRA_PARAMS toolset=mingw) #CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER}

IF(UNIX)
  SET(BOOST_LIB_SUFFIX "${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(BOOST_LIB_PREFIX "${CMAKE_STATIC_LIBRARY_PREFIX}")
ENDIF(UNIX)

IF(MINGW)
  SET(BOOST_LIB_PREFIX "${CMAKE_STATIC_LIBRARY_PREFIX}")

  STRING(REPLACE "." ";" COMPVER "${CFS_CXX_COMPILER_VER}")
  LIST(GET COMPVER 0 COMPVER_MAJOR)
  LIST(GET COMPVER 1 COMPVER_MINOR)

  SET(BOOST_LIB_SUFFIX "${CMAKE_STATIC_LIBRARY_SUFFIX}")

  SET(BOOST_EXTRA_PARAMS ${BOOST_EXTRA_PARAMS} target-os=windows architecture=x86 address-model=64)
ENDIF(MINGW)


SET(BOOST_DATE_TIME_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_date_time${BOOST_LIB_SUFFIX}")
SET(BOOST_FILESYSTEM_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_filesystem${BOOST_LIB_SUFFIX}")
SET(BOOST_IOSTREAMS_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_iostreams${BOOST_LIB_SUFFIX}")
SET(BOOST_PRG_EXEC_MONITOR_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_prg_exec_monitor${BOOST_LIB_SUFFIX}")
SET(BOOST_PROGRAM_OPTIONS_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_program_options${BOOST_LIB_SUFFIX}")
SET(BOOST_PYTHON_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_python${BOOST_LIB_SUFFIX}")
SET(BOOST_REGEX_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_regex${BOOST_LIB_SUFFIX}")
SET(BOOST_SERIALIZATION_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_serialization${BOOST_LIB_SUFFIX}")
SET(BOOST_SERIALIZATION_HDF5_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_serialization_hdf5${BOOST_LIB_SUFFIX}")
SET(BOOST_SIGNALS_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_signals${BOOST_LIB_SUFFIX}")
SET(BOOST_SYSTEM_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_system${BOOST_LIB_SUFFIX}")
SET(BOOST_TEST_EXEC_MONITOR_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_test_exec_monitor${BOOST_LIB_SUFFIX}")
SET(BOOST_THREAD_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_thread${BOOST_LIB_SUFFIX}")
SET(BOOST_UNIT_TEST_FRAMEWORK_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_unit_test_framework${BOOST_LIB_SUFFIX}")
SET(BOOST_WSERIALIZATION_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_wserialization${BOOST_LIB_SUFFIX}")

#-------------------------------------------------------------------------------
# The BOOST external project
#-------------------------------------------------------------------------------

IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(boost
    PREFIX "${BOOST_prefix}"
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
  ExternalProject_Add(boost
    PREFIX "${BOOST_prefix}"
    URL ${LOCAL_FILE}
    CONFIGURE_COMMAND ""
    PATCH_COMMAND ""
    BINARY_DIR ${BOOST_source}
    BUILD_COMMAND ./bootstrap.sh --without-libraries=python --prefix=${BOOST_install} ${BOOST_EXTRA_PARAMS}
    INSTALL_COMMAND ./b2 threading=multi install
  )

  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(boost cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${BOOST_prefix}
  )
  
  ExternalProject_Add_Step(boost post_install
    COMMAND ${CMAKE_COMMAND} -P "${PI}"
    DEPENDEES install )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(boost cfsdeps_zipToCache
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
SET(CFSDEPS ${CFSDEPS} boost)

