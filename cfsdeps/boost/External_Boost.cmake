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

SET(PFN_TEMPL "${CFS_SOURCE_DIR}/cfsdeps/boost/boost-patch.cmake.in")
SET(PFN "${BOOST_prefix}/boost-patch.cmake")
CONFIGURE_FILE("${PFN_TEMPL}" "${PFN}" @ONLY) 

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
  "https://boostorg.jfrog.io/artifactory/main/release/${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}.0/source/${BOOST_GZ}"
  "https://sourceforge.net/projects/boost/files/boost/${BOOST_MAJOR_VER}.${BOOST_MINOR_VER}.0/${BOOST_GZ}/download"
  "${CFS_DS_SOURCES_DIR}/boost/${BOOST_GZ}"
)

SET(MD5_SUM ${BOOST_MD5})

SET(DLFN "${BOOST_prefix}/boost-download.cmake") 
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_download.cmake.in" "${DLFN}" @ONLY)

# After the installation we copy to cfs
set(PI "${BOOST_prefix}/boost-post_install.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cfsdeps/boost/boost-post_install.cmake.in" "${PI}" @ONLY) 

#copy license
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/boost/license/" DESTINATION "${CFS_BINARY_DIR}/license/boost" )



PRECOMPILED_ZIP(PRECOMPILED_PCKG_FILE "boost" "${BOOST_VER}")
 
set(TMP_DIR "${BOOST_install}")

set(ZIPFROMCACHE "${BOOST_prefix}/boost-zipFromCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipFromCache.cmake.in" "${ZIPFROMCACHE}" @ONLY)

set(ZIPTOCACHE "${BOOST_prefix}/boost-zipToCache.cmake")
CONFIGURE_FILE("${CFS_SOURCE_DIR}/cmake_modules/cfsdeps_zipToCache.cmake.in" "${ZIPTOCACHE}" @ONLY)

#-------------------------------------------------------------------------------
# Determine paths of BOOST libraries.
#-------------------------------------------------------------------------------
SET(LD "${CFS_BINARY_DIR}/${LIB_SUFFIX}")

IF(UNIX)
  SET(Boost_LIBRARY_DIR "${CFS_BINARY_DIR}/${LIB_SUFFIX}" CACHE PATH "Boost library dir.")
  SET(BOOST_LIB_SUFFIX "${CMAKE_STATIC_LIBRARY_SUFFIX}")
  SET(BOOST_LIB_PREFIX "${CMAKE_STATIC_LIBRARY_PREFIX}")
ELSE()
  SET(Boost_LIBRARY_DIR "${CFS_BINARY_DIR}/${LIB_SUFFIX}/" CACHE PATH "Boost library dir.")

  SET(BOOST_VER_WIN "${BOOST_MAJOR_VER}_${BOOST_MINOR_VER}")
  SET(BOOST_LIB_PREFIX "lib")
  IF(BOOST_VER VERSION_LESS_EQUAL 1.72)
    SET(BOOST_MSVC_VERSION 15.0)
  ELSE()
    SET(BOOST_MSVC_VERSION ${CFS_MSVC_VERSION})
  ENDIF()

  IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    IF(DEBUG)
      SET(BOOST_LIB_SUFFIX "-iw-mt-gd-x64-${BOOST_VER_WIN}.lib")
    ELSE()
      SET(BOOST_LIB_SUFFIX "-iw-mt-x64-${BOOST_VER_WIN}.lib")
    ENDIF(DEBUG)
  ELSE()
    IF(BOOST_MSVC_VERSION VERSION_GREATER_EQUAL 16.0)
      IF(DEBUG)
        SET(BOOST_LIB_SUFFIX "-vc142-mt-gd-x64-${BOOST_VER_WIN}.lib")
      ELSE()
        SET(BOOST_LIB_SUFFIX "-vc142-mt-x64-${BOOST_VER_WIN}.lib")
      ENDIF()
    ELSE()
      IF(DEBUG)
        SET(BOOST_LIB_SUFFIX "-vc141-mt-gd-x64-${BOOST_VER_WIN}.lib")
      ELSE()
        SET(BOOST_LIB_SUFFIX "-vc141-mt-x64-${BOOST_VER_WIN}.lib")
      ENDIF()
    ENDIF()
  ENDIF()
ENDIF()
MARK_AS_ADVANCED(Boost_LIBRARY_DIR)

# boost configuration
set(BOOST_BOOTSTRAP_PARAMS "") 
set(BOOST_PATCH_COMMAND "") # to set user-config.jam for toolset

#SET(BOOST_JAM_PATCH_COMMAND "")
SET(BOOST_B2_PARAMS "")

IF(UNIX)
  # see closed issue Boost build error nach Systemupgrade auf Ubuntu 18.04 in Group Stingl
  set(BOOST_B2_PARAMS cxxflags=-fPIC)
  if(DEBUG)
    set(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} variant=debug)
  else()
    set(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} variant=release)
  endif()  
  get_filename_component(BOOST_DEP_LIBPATH ${ZLIB_LIBRARY} DIRECTORY) # get the library path without filename from library location
  set(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "-sZLIB_INCLUDE=${CFS_BINARY_DIR}/include" "-sZLIB_LIBPATH=${BOOST_DEP_LIBPATH}" "-sNO_ZLIB=0")

  # setting the mac system clang compiler as toolset does not work
  set(BOOST_SYSTEM_COMPILER OFF)
  if(CFS_DISTRO STREQUAL "MACOSX" AND CFS_CXX_COMPILER_NAME STREQUAL "CLANG")
    set(BOOST_SYSTEM_COMPILER ON)
  endif()
  # set the gcc and clang compiler via toolset it might be specified via CXX=g++-8 or CXX=clang++
  if(NOT BOOST_SYSTEM_COMPILER)
    if(CFS_CXX_COMPILER_NAME STREQUAL "GCC")
      set(BOOST_TOOLSET_COMPILER "gcc")
    endif()
    if(CFS_CXX_COMPILER_NAME STREQUAL "CLANG")
      set(BOOST_TOOLSET_COMPILER "clang")
    endif()
    if(BOOST_TOOLSET_COMPILER)
      set(BOOST_PATCH_COMMAND echo "using ${BOOST_TOOLSET_COMPILER} : ${CFS_CXX_COMPILER_VER} : ${CMAKE_CXX_COMPILER} $<SEMICOLON>" > user-config.jam)
      set(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "--user-config=user-config.jam" "toolset=${BOOST_TOOLSET_COMPILER}-${CFS_CXX_COMPILER_VER}")
    endif()  
  endif()
  IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    # Simon made it without toolset and it shall work.
    set(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} --with-toolset=intel-linux )
    set(BOOST_B2_PARAMS ${BOOST_B2_PARAMS})
  ENDIF()
ENDIF(UNIX)

IF(WIN32)

  SET(BOOST_BOOTSTRAP_PARAMS "--without-libraries=python" "--prefix=${BOOST_install}" "--target-os=windows" "--architecture=x86" "--address-model=64" "--with-filesystem")
  SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "target-os=windows" "architecture=x86" "address-model=64" "--prefix=${BOOST_install}")
  SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "variant=debug,release" "threading=multi" "link=static,shared" "runtime-link=shared")

  IF(BOOST_MSVC_VERSION VERSION_GREATER_EQUAL 16.0)
    IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
      IF(CFS_CXX_COMPILER_VER GREATER_EQUAL 20.2)
        SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=intel-2021.2-vc14.2")
        SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=intel-2021.2-vc14.2")
      ELSEIF(CFS_CXX_COMPILER_VER GREATER_EQUAL 20.1)
        SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=intel-2021.1-vc14.2")
        SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=intel-2021.1-vc14.2")
      ELSE()
        SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=intel-19.1-vc14.2")
        SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=intel-19.1-vc14.2")
      ENDIF()
    ELSE()
      SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=msvc-14.2")
      SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=msvc-14.2")
    ENDIF()
  ELSE()
    IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
      SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=intel-19.1-vc14.1")
      SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=intel-19.1-vc14.1")
    ELSE()
      SET(BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS} "--with-toolset=msvc-14.1")
      SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "toolset=msvc-14.1")
    ENDIF()
  ENDIF()

  SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "--with-filesystem" "--with-date_time" "--with-iostreams" "--with-log" "--with-program_options" "--with-regex" "--with-serialization" "--with-thread" "--with-chrono" "--with-test" "--with-random")

  SET(BOOST_B2_PARAMS ${BOOST_B2_PARAMS} "-sZLIB_SOURCE=${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/zlib/src/zlib" "-sBZIP2_SOURCE=${CMAKE_CURRENT_BINARY_DIR}/cfsdeps/bzip2/src/bzip2" "install")
ENDIF(WIN32)

#message("CFS_CXX_COMPILER_NAME ${CFS_CXX_COMPILER_NAME}")
#message("BOOST_SYSTEM_COMPILER ${BOOST_SYSTEM_COMPILER}")
#message("BOOST_PATCH_COMMAND ${BOOST_PATCH_COMMAND}")
#message("BOOST_B2_PARAMS ${BOOST_B2_PARAMS}")
#message("BOOST_BOOTSTRAP_PARAMS ${BOOST_BOOTSTRAP_PARAMS}")

SET(BOOST_DATE_TIME_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_date_time${BOOST_LIB_SUFFIX}")
SET(BOOST_FILESYSTEM_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_filesystem${BOOST_LIB_SUFFIX}")
SET(BOOST_IOSTREAMS_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_iostreams${BOOST_LIB_SUFFIX}")
SET(BOOST_PRG_EXEC_MONITOR_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_prg_exec_monitor${BOOST_LIB_SUFFIX}")
SET(BOOST_PROGRAM_OPTIONS_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_program_options${BOOST_LIB_SUFFIX}")
SET(BOOST_PYTHON_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_python${BOOST_LIB_SUFFIX}")
SET(BOOST_REGEX_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_regex${BOOST_LIB_SUFFIX}")
SET(BOOST_SERIALIZATION_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_serialization${BOOST_LIB_SUFFIX}")
SET(BOOST_SERIALIZATION_HDF5_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_serialization_hdf5${BOOST_LIB_SUFFIX}")
SET(BOOST_TEST_EXEC_MONITOR_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_test_exec_monitor${BOOST_LIB_SUFFIX}")
SET(BOOST_THREAD_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_thread${BOOST_LIB_SUFFIX}")
SET(BOOST_UNIT_TEST_FRAMEWORK_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_unit_test_framework${BOOST_LIB_SUFFIX}")
SET(BOOST_WSERIALIZATION_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_wserialization${BOOST_LIB_SUFFIX}")
SET(BOOST_LOG_LIB "${Boost_LIBRARY_DIR}/${BOOST_LIB_PREFIX}boost_log${BOOST_LIB_SUFFIX}")

# in case someone manages to break cross-compilation again, this is how it works:
# echo "using gcc : mingw : x86_64-w64-mingw32-gcc ;" > user-config.jam
# ./bootstrap.sh --without-libraries=python --prefix=<path> --with-toolset=gcc target-os=windows architecture=x86 address-model=64 release
# ./b2 --user-config=user-config.jam toolset=gcc-mingw target-os=windows architecture=x86 address-model=64 release

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
    BUILD_BYPRODUCTS "${BOOST_DATE_TIME_LIB}"
  )
ELSE()
  #-------------------------------------------------------------------------------
  # If precompiled package does not exist build external project
  #-------------------------------------------------------------------------------
  IF(UNIX)
    ExternalProject_Add(boost
      PREFIX "${BOOST_prefix}"
      URL ${LOCAL_FILE}
      PATCH_COMMAND ${BOOST_PATCH_COMMAND}
      CONFIGURE_COMMAND ""
      BINARY_DIR ${BOOST_source}
      # newer boost has system and signals head only.
      # full list via bootstrap.sh --show-libraries
      # 1.73.0: atomic,chrono,container,context,contract,coroutine,date_time,exception,fiber,filesystem,graph,graph_parallel,headers,iostreams,locale,log,math,mpi,nowide,program_options,random,regex,serialization,stacktrace,system,test,thread,timer,type_erasure,wave  
      BUILD_COMMAND ./bootstrap.sh --with-libraries=date_time,filesystem,iostreams,log,program_options,regex,serialization,thread,chrono,test --prefix=${BOOST_install} ${BOOST_BOOTSTRAP_PARAMS} 
      INSTALL_COMMAND ./b2  ${BOOST_B2_PARAMS} define=BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX link=static threading=multi runtime-link=static install --no-cmake-config
    )
  ELSE()
    # message("BOOST_install: ${BOOST_install}")
    # message("BOOST_BOOTSTRAP_PARAMS:  ${BOOST_BOOTSTRAP_PARAMS}")
    # message("BOOST_B2_PARAMS: ${BOOST_B2_PARAMS}")
    ExternalProject_Add(boost
      PREFIX "${BOOST_prefix}"
      URL ${LOCAL_FILE}
      #      PATCH_COMMAND ${BOOST_JAM_PATCH_COMMAND}
      PATCH_COMMAND ""
      CONFIGURE_COMMAND ""
      BINARY_DIR ${BOOST_source}
      BUILD_COMMAND ./bootstrap.bat ${BOOST_BOOTSTRAP_PARAMS}
      INSTALL_COMMAND ./b2 ${BOOST_B2_PARAMS}
    )
  ENDIF()

  #-------------------------------------------------------------------------------
  # Add custom patch step
  #-------------------------------------------------------------------------------
  IF(CFS_CXX_COMPILER_NAME STREQUAL "ICC")
    ExternalProject_Add_Step(boost winpatch
      COMMAND ${CMAKE_COMMAND} -P "${PFN}"
      DEPENDERS build
      DEPENDEES download
      DEPENDS "${PFN}"
      WORKING_DIRECTORY ${BOOST_source}
    )
  ENDIF()

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

ADD_DEPENDENCIES(boost zlib)
