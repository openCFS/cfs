# An enormous rich set of C++ addones. Mostly header only (> 170 MB includes). Mandatory
# https://www.boost.org

clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "boost")
set(PACKAGE_VER "1.91.0")  
# with boost 1.91.tar.bz2 we have unpacking issues in the RHEL8 pipeline
set(PACKAGE_FILE "boost_1_91_0.tar.gz") # does not reflect PACKAGE_VER style
set(PACKAGE_MD5 "e799ed3e5af9708739fb2e088c670ae1") # 1.91.0
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.
  
# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://archives.boost.io/release/${PACKAGE_VER}/source/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "IntelLLV")
  # see TOOLSET setting below
  force_c_cxx_compiler("cl") # will set FORCE_C_CXX and do NOT set the default compiler options!
  assert_set(FORCE_C_CXX)
else()
  use_c_and_fortran(ON OFF)
endif() 

# annoying bjam ignores CXX, it is a pain in the ass :(
# if you have issues with a specific compiler, usually any compiler works to build the lib - most ist heady only anyway
# see https://www.intel.com/content/www/us/en/developer/articles/technical/building-boost-with-oneapi.html
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang") # matches also AppleClang
  set(TOOLSET clang)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(TOOLSET gcc)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM" AND UNIX) # Linux as there is no oneAPI for APPLE any more
  set(TOOLSET intel-linux) 
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM" AND WIN32) # Linux and Windows
  # works for intel-win32 in theory but not necessarily in reality (error: Supported msvc versions not known for intel win32)
  # we use force_c_cxx_compiler() for proper .zip naming  
  set(TOOLSET msvc) 
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel") # no more IntelLLVM but Linux and Windows
  set(TOOLSET intel) # seems to first check for icpx (oneAPI LLVM based), then icpcp (classic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(TOOLSET msvc)
else()
  message(FATAL "no specific boost toolset for ${CMAKE_CXX_COMPILER_ID} implemented")
endif()

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# generate BOOST_LIBRARY=PACKAGE_LIBARAY with os specific list of static libs. 
set_package_library_list_lib_prefix("boost_atomic;boost_iostreams;boost_program_options;boost_thread;boost_log;boost_regex;boost_unit_test_framework;boost_log_setup;boost_serialization") 

# to not always need to have all libs from BOOST_LIBRARY - no need for cache variables
set(BOOST_SERIALIZATION_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_serialization${CMAKE_STATIC_LIBRARY_SUFFIX})
set(BOOST_THREAD_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_thread${CMAKE_STATIC_LIBRARY_SUFFIX})

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# create user-config.jam from zlib-config.jam.in to make sure boost uses our zlib and does not search the system
assert_set(ZLIB_VER)
assert_set(DEPS_PREFIX)
configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/zlib-config.jam.in" "${DEPS_PREFIX}/user-config.jam" @ONLY)

if(UNIX)
  set(BOOTSTRAP ./bootstrap.sh)
else()
  set(BOOTSTRAP ./bootstrap.bat)
endif()

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)
# create zlib-config.jam. Only necessary for Windows. See https://www.boost.org/doc/libs/1_81_0/tools/build/doc/html/index.html

# filter from DEPS_INSTALL, zip and copy to binary dir 
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# dump_dependencycd _variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()
# if not, build newly and possibly pack the stuff
else()
  if(WIN32)
    # compile cfs with BOOST_ALL_NO_LIB  and make sure _WIN32_WINNT matches cfs 
    set(DEFINE "define=_WIN32_WINNT=0x0A00") # required as we build as system. https://github.com/boostorg/log/issues/172
  else()
    set(DEFINE "define=BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX") 
  endif()    

  # Windows (cl and icx) needs matching deps to prevent linker errors but debugging on Windows sucks anyway
  if(DEBUG AND WIN32)
    set(_BUILD "debug")
    set(_DS "on")
  else()
    set(_BUILD "release")
    set(_DS "off")
  endif()

  assert_unset(PATCHES_SCRIPT)

  # more and more of boost is header only. Add here the libs which are not required for cfs and some cfsdeps which also use boost
  set(WITHOUT --without-python --without-graph_parallel --without-mpi --without-container --without-context --without-json --without-math --without-nowide --without-contract --without-coroutine --without-graph --without-stacktrace --without-fiber --without-wave --without-url --no-cmake-config) 

  # we need to build the package - here in cmake style
  ExternalProject_Add(${PACKAGE_NAME}
    PREFIX ${DEPS_PREFIX}
    BINARY_DIR ${DEPS_SOURCE}
    URL ${PACKAGE_MIRRORS}
    URL_MD5 ${PACKAGE_MD5}
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}
    # in case the mirrors have different file names we always store to the same
    DOWNLOAD_NAME ${PACKAGE_FILE}
    DOWNLOAD_NO_PROGRESS ON 
    PATCH_COMMAND ""
    COMMAND ${CMAKE_COMMAND} -E copy "${DEPS_PREFIX}/user-config.jam" "${DEPS_SOURCE}"
    # build b2 with the system compiler. In the pipeline (only) intel-linux fails for missing pthread
    CONFIGURE_COMMAND ${BOOTSTRAP}
    # on Windows calling b2 might result in security issues (access violation)
    # --threading=multi seems to work even with USE_OPENMP=OFF on all systems?!
    BUILD_COMMAND ./b2 cxxstd=17 --toolset=${TOOLSET} --user-config=user-config.jam ${WITHOUT} --layout=system --prefix=${DEPS_INSTALL} ${DEFINE} link=static address-model=64 threading=multi runtime-link=shared variant=${_BUILD} debug-symbols=${_DS} install
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY}
    # Wrap build step in script to log output, since it's super long and clutters the pipeline
    # see https://cmake.org/cmake/help/v3.0/module/ExternalProject.html
    LOG_BUILD 1) # see <build>/cfsdeps/boost/src/boost-stamp/boost-build-out.log

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})

# we depend on zlib, don't forget to have ZLIB_LIBRARY in LL_TARGET after BOOST_LIBRARY
add_dependencies(boost zlib)
