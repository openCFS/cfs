#-------------------------------------------------------------------------------
# Computational Geometry Algorithms Library
#
# Project Homepage
# http://www.cgal.org/
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set paths to cgal sources according to ExternalProject.cmake 
#-------------------------------------------------------------------------------
set(cgal_prefix  "${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/cgal")
set(cgal_source  "${cgal_prefix}/src/cgal")
set(cgal_install  "${CMAKE_CURRENT_BINARY_DIR}")

string(TOLOWER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
IF(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Debug")
ELSE(CMAKE_BUILD_TYPE STREQUAL "debug")
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF(CMAKE_BUILD_TYPE STREQUAL "debug")

# Make sure that static Boost libs are used under Windows.
SET(CGAL_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -DBOOST_THREAD_USE_LIB")

IF(CFS_CXX_COMPILER_NAME STREQUAL "OPEN64")
  SET(CGAL_CXX_FLAGS "${CGAL_CXX_FLAGS} -mieee-fp -fp-accuracy=strict -DCGAL_DISABLE_ROUNDING_MATH_CHECK")
ENDIF()

SET(CMAKE_ARGS
  -DCMAKE_COLOR_MAKEFILE:BOOL=${CMAKE_COLOR_MAKEFILE}
  -DCMAKE_INSTALL_PREFIX:PATH=${cgal_install}
  -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
  -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS:STRING=${CGAL_CXX_FLAGS}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_RANLIB:FILEPATH=${CMAKE_RANLIB}
)

IF(CMAKE_TOOLCHAIN_FILE)
  LIST(APPEND CMAKE_ARGS
    -DCMAKE_TOOLCHAIN_FILE:FILEPATH=${CMAKE_TOOLCHAIN_FILE}
  )
ENDIF()

IF(MINGW)
  IF(CFS_ARCH STREQUAL "X86_64")
    LIST(APPEND CMAKE_ARGS
      -C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/cgal/TryRunResults_GCC45_CENTOS6_WIN64.cmake
      )
  ENDIF()
ENDIF(MINGW)

IF(APPLE AND CMAKE_CROSSCOMPILING)
  IF(CFS_ARCH STREQUAL "X86_64")
    LIST(APPEND CMAKE_ARGS
      -C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/cgal/TryRunResults_MACOSX_X86_64_CENTOS6.cmake
      )
  ELSEIF(CFS_ARCH STREQUAL "I386")
    LIST(APPEND CMAKE_ARGS
      -C ${CMAKE_CURRENT_SOURCE_DIR}/cfsdeps/cgal/TryRunResults_MACOSX_I386_CENTOS6.cmake
      )
  ENDIF()
ENDIF()

#-------------------------------------------------------------------------------
# Set names of patch file and template file.
#-------------------------------------------------------------------------------
SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/cgal-patch.cmake.in")
SET(PFN "${cgal_prefix}/cgal-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

SET(BOOST_SETUP_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/cgal/CGAL_SetupBoost.cmake.in")
SET(BOOST_SETUP "${cgal_prefix}/CGAL_SetupBoost.cmake")
CONFIGURE_FILE("${BOOST_SETUP_TEMPL}" "${BOOST_SETUP}" @ONLY) 

#-------------------------------------------------------------------------------
# Set up a list of publicly available mirrors, since the non-standard port 
# number of the FTP server on the CFS++ development server  may not be
# accessible from behind firewalls.
# Also set name of local file in CFS_DEPS_CACHE_DIR and MD5_SUM which will be
# used to configure the download CMake file for the library.
#-------------------------------------------------------------------------------
SET(MIRRORS
  "http://archive.ubuntu.com/ubuntu/pool/universe/c/cgal/cgal_4.2.orig.tar.bz2"
  "https://gforge.inria.fr/frs/download.php/32361/${CGAL_BZ2}"
  "${CGAL_URL}/${CGAL_BZ2}"
)
SET(LOCAL_FILE "${CFS_DEPS_CACHE_DIR}/sources/cgal/${CGAL_BZ2}")
SET(MD5_SUM ${CGAL_MD5})

SET(DLFN "${cgal_prefix}/cgal-download.cmake")
CONFIGURE_FILE(
  "${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in"
  "${DLFN}"
  @ONLY
)

PRECOMPILED_ZIP_CXX(PRECOMPILED_PCKG_FILE "cgal" "${CGAL_VER}")  
  
# This should be either PREFIX_DIR (install manifest is used for zipping)
# or INSTALL_DIR (install directory will be zipped)
SET(TMP_DIR "${cgal_prefix}")

SET(ZIPFROMCACHE "${cgal_prefix}/cgal-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

SET(ZIPTOCACHE "${cgal_prefix}/cgal-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# The CGAL external project
#-------------------------------------------------------------------------------
IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON" AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  #-------------------------------------------------------------------------------
  # If precompiled package exists copy files from cache
  #-------------------------------------------------------------------------------
  ExternalProject_Add(cgal
    PREFIX "${cgal_prefix}"
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
  ExternalProject_Add(cgal
    DEPENDS boost zlib gmp mpfr
    PREFIX "${cgal_prefix}"
    URL ${LOCAL_FILE}
    URL_MD5 ${CGAL_MD5}
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PFN}"
    CMAKE_ARGS
      ${CMAKE_ARGS}
      -DCGAL_INSTALL_BIN_DIR:PATH=bin/${CFS_ARCH_STR}
      -DCGAL_INSTALL_CMAKE_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}/CGAL
      -DCGAL_INSTALL_LIB_DIR:PATH=${LIB_SUFFIX}/${CFS_ARCH_STR}
      -DBUILD_SHARED_LIBS:BOOL=OFF
#      -DBoost_INCLUDE_DIR:PATH=${cgal_install}/include
#      -DBoost_LIBRARY_DIRS:PATH=${cgal_install}/${LIB_SUFFIX}/$ARCH_STR
#      -DBoost_THREAD_LIBRARY:FILEPATH=${BOOST_THREAD_LIB}
#      -DBoost_THREAD_LIBRARY_RELEASE:FILEPATH=${BOOST_THREAD_LIB_RELEASE}
#      -DBoost_THREAD_LIBRARY_DEBUG:FILEPATH=${BOOST_THREAD_LIB_DEBUG}
      -DWITH_CGAL_Qt3:BOOL=OFF
      -DWITH_CGAL_Qt4:BOOL=OFF
      -DZLIB_INCLUDE_DIR:PATH=${ZLIB_INCLUDE_DIR}
      -DZLIB_LIBRARY:PATH=${ZLIB_LIBRARY}
      -DGMP_INCLUDE_DIR:PATH=${GMP_INCLUDE_DIR}
      -DGMP_LIBRARIES:FILEPATH=${GMP_LIBRARY}
      -DGMP_LIBRARIES_DIR:PATH=${cgal_install}/${LIB_SUFFIX}/${CFS_ARCH_STR}
      -DWITH_GMP:BOOL=ON
      -DMPFR_INCLUDE_DIR:PATH=${MPFR_INCLUDE_DIR}
      -DMPFR_LIBRARIES:FILEPATH=${MPFR_LIBRARY}
      -DWITH_MPFR:BOOL=ON
  )
  
  #-------------------------------------------------------------------------------
  # Add custom download step to be able to download from a list of mirrors
  # instead of just a single URL.
  #-------------------------------------------------------------------------------
  ExternalProject_Add_Step(cgal cfsdeps_download
    COMMAND ${CMAKE_COMMAND} -P "${DLFN}"
    DEPENDERS download
    DEPENDS "${DLFN}"
    WORKING_DIRECTORY ${cgal_prefix}
  )
  
  IF("${CFS_DEPS_PRECOMPILED}" STREQUAL "ON")
    #-------------------------------------------------------------------------------
    # Add custom step to zip a precompiled package to the cache.
    #-------------------------------------------------------------------------------
    ExternalProject_Add_Step(cgal cfsdeps_zipToCache
      COMMAND ${CMAKE_COMMAND} -P "${ZIPTOCACHE}"
      DEPENDEES install
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
  cgal
)

SET(CGAL_INCLUDE_DIR "${CFS_BINARY_DIR}/include")

#-------------------------------------------------------------------------------
# Determine paths of CGAL libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}/${CFS_ARCH_STR}")
SET(CGAL_LIBRARY
  "${LD}/libCGAL.a"
  CACHE FILEPATH "CGAL library.")

MARK_AS_ADVANCED(CGAL_INCLUDE_DIR)
MARK_AS_ADVANCED(CGAL_LIBRARY)
