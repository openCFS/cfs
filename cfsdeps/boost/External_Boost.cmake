# An enormous rich set of C++ addones. Mostly header only (> 170 MB includes). Mandatory
# https://www.boost.org

clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "boost")
# note that any newer version than 1.78.0 causes > 30 wrong test results, e.g. ExpressionHeatSource
# probably in conjunction with muparser 2.2.6. There is a branch upgrade_boost which contains changes to compille 1.84
set(PACKAGE_VER "1.78.0")  
set(PACKAGE_FILE "boost_1_78_0.tar.bz2") # does not reflect PACKAGE_VER style
set(PACKAGE_MD5 "db0112a3a37a3742326471d20f1a186a") # 1.78.0
#set(PACKAGE_VER "1.71.0")  
#set(PACKAGE_FILE "boost_1_71_0.tar.bz2") # does not reflect PACKAGE_VER style
#set(PACKAGE_MD5 "4cdf9b5c2dc01fb2b7b733d5af30e558") # 1.71.0

set(DEPS_VER "-c") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.
  
# the mirrors can point to arbitrary file names. 
#set(PACKAGE_MIRRORS "https://boostorg.jfrog.io/artifactory/main/release/${PACKAGE_VER}/source/${PACKAGE_FILE}")
set(PACKAGE_MIRRORS "https://archives.boost.io/release/${PACKAGE_VER}/source/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# the following generic code could set a specific compiler, but the question is, how important it is to have the
# few libs compiled with the exact compiler? Most of boost is header only and most important is that openCFS compiles :)
# PATCH_COMMAND echo "using <clang/gcc/...> : ${CFS_CXX_COMPILER_VER} : ${CMAKE_CXX_COMPILER} $<SEMICOLON>" > user-config.jam
# in case change cfs to version in zlib-toolset-config.jam.in - bjam is unfortuantely a real pain and horrible compared to cmake!
# set(TOOLSET --user-config=user-config.jam toolset=<clang/gcc/...>)
# see https://www.intel.com/content/www/us/en/developer/articles/technical/building-boost-with-oneapi.html
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang") # matches also AppleClang
  set(TOOLSET toolset=clang)
  # macOS Sonoma 14.4 fails to build without flags 
  set(B2_ARGS "cxxflags=-std=c++14 -Wno-enum-constexpr-conversion")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(TOOLSET toolset=gcc)
  # for APPLE and gcc we would need set(TOOLSET_NAME gcc) but the have the version in zlib-toolset-config.jam.in   
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM") # Linux and Windows
  # https://www.boost.org/doc/libs/1_81_0/tools/build/doc/html/index.html
  set(TOOLSET_NAME clang) # this triggers to use zlib-toolset-config.jam.in instead of zlib-config.jam.in
  set(TOOLSET toolset=clang-cfs) # check zlib-toolset-config.jam.in, the cfs is hard coded there
  set(B2_ARGS "cxxflags=-Wno-enum-constexpr-conversion") # disable compiler error on warning
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel") # no more IntelLLVM but Linux and Windows
  set(TOOLSET toolset=intel) # seems to first check for icpx (oneAPI LLVM based), then icpcp (classic)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  set(TOOLSET toolset=msvc)
else()
  message(STATUS "no specific boost toolset for ${CMAKE_CXX_COMPILER_ID} implemented, could work anyway ...")
endif()

 # we'll disable fortran for lis by not using saamg which is fast, but very sensitive to system condition
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# generate BOOST_LIBRARY=PACKAGE_LIBARAY with os specific list of static libs. 
set_package_library_list_lib_prefix("boost_atomic;boost_iostreams;boost_program_options;boost_thread;boost_chrono;boost_log;boost_regex;boost_unit_test_framework;boost_date_time;boost_log_setup;boost_serialization;boost_filesystem") 

# to not always need to have all libs from BOOST_LIBRARY - no need for cache variables
set(BOOST_SERIALIZATION_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_serialization${CMAKE_STATIC_LIBRARY_SUFFIX})
set(BOOST_THREAD_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_thread${CMAKE_STATIC_LIBRARY_SUFFIX})

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# create user-config.jam which is either zlib-config.jam.in or zlib-toolset-config.jam.in
assert_set(ZLIB_VER)
assert_set(DEPS_PREFIX)
if(TOOLSET_NAME)
  configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/zlib-toolset-config.jam.in" "${DEPS_PREFIX}/user-config.jam" @ONLY)
else()
  configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/zlib-config.jam.in" "${DEPS_PREFIX}/user-config.jam" @ONLY)
endif()

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

# dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # even with the list in boostrap a lot of unnecessary stuff would be build
  # 1.81.0 set(WITHOUT --without-python --without-graph_parallel --without-mpi --without-container --without-context --without-json --without-math --without-nowide --without-contract --without-coroutine --without-graph --without-stacktrace --without-fiber --without-wave --without-url --no-cmake-config) 
  set(WITHOUT --without-python --without-graph_parallel --without-mpi --without-container --without-context --without-math --without-contract --without-coroutine --without-graph --without-stacktrace --without-fiber --without-wave  --no-cmake-config)

  if(WIN32)
    # compile cfs with BOOST_ALL_NO_LIB  and make sure _WIN32_WINNT matches cfs 
    set(DEFINE "define=_WIN32_WINNT=0x0A00") # required as we build as system. https://github.com/boostorg/log/issues/172
  else()
    set(DEFINE "define=BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX") 
  endif()    

  # Windows (cl and icx) needs matching deps to prevent linker errors
  if(DEBUG AND WIN32)
    set(_BUILD "debug")
    set(_DS "on")
  else()
    set(_BUILD "release")
    set(_DS "off")
  endif()

  # some patches are required
  generate_patches_script()

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
    PATCH_COMMAND ${CMAKE_COMMAND} -P "${PATCHES_SCRIPT}"
    COMMAND ${CMAKE_COMMAND} -E copy "${DEPS_PREFIX}/user-config.jam" "${DEPS_SOURCE}"
    # we call bootstap without the system compiler. --with-libraries seems to have no effect (all be default)
    CONFIGURE_COMMAND ${BOOTSTRAP} 
    # on Windows calling b2 might result in security issues (access violation)
    # --threading=multi seems to work even with USE_OPENMP=OFF on all systems?!
    BUILD_COMMAND ./b2 --user-config=user-config.jam ${WITHOUT} --layout=system --prefix=${DEPS_INSTALL} ${TARGET} ${DEFINE} ${B2_ARGS} link=static address-model=64 threading=multi runtime-link=shared variant=${_BUILD} debug-symbols=${_DS} install
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${PACKAGE_LIBRARY}
    # Wrap build step in script to log output, since it's super long and clutters the pipeline
    # see https://cmake.org/cmake/help/v3.0/module/ExternalProject.html
    LOG_BUILD 1 ) # see <build>/cfsdeps/boost/src/boost-stamp/boost-build-out.log

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

# we depend on zlib
add_dependencies(boost zlib)

# with old boost 1.78.0 we have issue with C++17 with the following workaround. Remove with more recent boost! 
if(APPLE)
  # https://stackoverflow.com/questions/77133361/no-template-named-unary-function-in-namespace-std-did-you-mean-unary-fun
  add_compile_definitions(_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION)
endif()
# boost 1.78 by default still uses std::unary_function which is removed in C++17
add_compile_definitions(BOOST_NO_CXX98_FUNCTION_BASE)

# ---------------------------------------------------------------------------
# Boost for preCICE (only when preCICE is built from source: CFS_BUILD_PRECICE).
#
# preCICE 3.x is always a SHARED library and links *shared* boost (libprecice.so
# gets DT_NEEDED libboost_*.so); it also uses find_package(Boost CONFIG). The
# default cfsdeps boost above is static and built with --no-cmake-config, so
# preCICE cannot use it. The default build above is left completely untouched
# (same version 1.78, used by cfs and everything else that pins it).
#
# Detect-or-build: if a suitable *shared* boost (with a CMake config) is already
# available (system / build image), preCICE links that and we build nothing.
# Otherwise build a shared copy of the SAME version here, into its own prefix,
# and install the .so into the cfs build's lib/ so cfs finds it at runtime
# (rpath $ORIGIN/../lib) and it ships in the package. iostreams uses the system
# zlib (universally present) to avoid linking the non-PIC cfsdeps zlib.
#
# NOTE: not validated in this environment - verify on the build machine.
# ---------------------------------------------------------------------------
if(CFS_BUILD_PRECICE AND UNIX)
  # 1) detect an already-available shared boost with the components preCICE needs
  set(_PRECICE_BOOST_SAVE_STATIC "${Boost_USE_STATIC_LIBS}")
  set(Boost_USE_STATIC_LIBS OFF)
  find_package(Boost 1.74 CONFIG QUIET COMPONENTS log log_setup program_options thread)
  set(Boost_USE_STATIC_LIBS "${_PRECICE_BOOST_SAVE_STATIC}")

  if(Boost_FOUND)
    # preCICE's own find_package(Boost CONFIG) will locate the same system boost.
    message(STATUS "preCICE: found shared Boost ${Boost_VERSION} (${Boost_DIR}) - not building a boost copy")
    set(PRECICE_BOOST_ROOT "" CACHE INTERNAL "boost root hint for preCICE (empty = use system)" FORCE)
  else()
    message(STATUS "preCICE: no suitable shared Boost found - building shared boost ${PACKAGE_VER} in cfsdeps")
    set(BOOST_PIC_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/boost-pic")
    set(BOOST_PIC_INSTALL "${BOOST_PIC_PREFIX}/install")
    set(PRECICE_BOOST_ROOT "${BOOST_PIC_INSTALL}" CACHE INTERNAL "boost root hint for preCICE" FORCE)
    set(BOOST_PIC_ZIP "${CFS_DEPS_CACHE_DIR}/precompiled/boost-pic_${PACKAGE_VER}${DEPS_VER}_${CFS_ARCH_STR}_C-${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}.tar.gz")

    if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${BOOST_PIC_ZIP}")
      # reuse the cached shared boost: unpack it into the install prefix
      ExternalProject_Add(boost-pic
        PREFIX ${BOOST_PIC_PREFIX}
        DOWNLOAD_COMMAND "" CONFIGURE_COMMAND "" BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${BOOST_PIC_INSTALL}
        COMMAND ${CMAKE_COMMAND} -E chdir ${BOOST_PIC_INSTALL} ${CMAKE_COMMAND} -E tar xzf ${BOOST_PIC_ZIP})
    else()
      # shared build WITH the boost CMake config (preCICE uses find_package CONFIG).
      # No --user-config: b2 uses the toolset compiler from PATH and the system zlib
      # for iostreams (the cfsdeps zlib is static/non-PIC and cannot go into a .so).
      set(_PIC_WITHOUT --without-python --without-graph_parallel --without-mpi --without-container
                       --without-context --without-math --without-contract --without-coroutine
                       --without-graph --without-stacktrace --without-fiber --without-wave)
      ExternalProject_Add(boost-pic
        PREFIX ${BOOST_PIC_PREFIX}
        BINARY_DIR ${BOOST_PIC_PREFIX}/src/boost-pic
        URL ${PACKAGE_MIRRORS}
        URL_MD5 ${PACKAGE_MD5}
        DOWNLOAD_DIR ${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}
        DOWNLOAD_NAME ${PACKAGE_FILE}
        DOWNLOAD_NO_PROGRESS ON
        CONFIGURE_COMMAND ${BOOTSTRAP}
        BUILD_COMMAND ./b2 ${_PIC_WITHOUT} --layout=system --prefix=${BOOST_PIC_INSTALL} ${TOOLSET}
                      define=BOOST_UUID_RANDOM_PROVIDER_FORCE_POSIX
                      cflags=-fPIC cxxflags=-fPIC link=shared address-model=64
                      threading=multi variant=release install
        INSTALL_COMMAND ""
        LOG_BUILD 1)
      if(${CFS_DEPS_PRECOMPILED})
        ExternalProject_Add_Step(boost-pic store-precompiled
          COMMAND ${CMAKE_COMMAND} -E make_directory ${CFS_DEPS_CACHE_DIR}/precompiled
          COMMAND ${CMAKE_COMMAND} -E chdir ${BOOST_PIC_INSTALL} ${CMAKE_COMMAND} -E tar czf ${BOOST_PIC_ZIP} .
          DEPENDEES install)
      endif()
    endif()

    # make the shared boost .so available in the cfs build lib/ (runtime rpath + packaging)
    ExternalProject_Add_Step(boost-pic copy-so-to-lib
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/lib
      COMMAND sh -c "cp -a ${BOOST_PIC_INSTALL}/lib/libboost_*.so* ${CMAKE_BINARY_DIR}/lib/ 2>/dev/null || true"
      DEPENDEES install)

    set(CFSDEPS ${CFSDEPS} boost-pic)
  endif()
endif()
