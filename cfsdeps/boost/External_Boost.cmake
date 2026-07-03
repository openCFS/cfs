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

# openCFS usually builds all dependencies static for a self-contained package.
# preCICE (libprecice.so, always a shared lib) however links *shared* boost and
# locates it via find_package(Boost CONFIG), which the static build cannot serve.
# So when preCICE is built from source the strategy of THIS single boost build
# switches to shared (implies -fPIC) with the boost CMake config files; cfs and
# preCICE then use the same boost. The .so files are copied to <build>/lib and
# shipped in the package like DLLs on Windows (redistributables.cmake installs
# lib/, the cfs binary finds them via rpath $ORIGIN/../lib).
# DEPS_ID gives the shared variant its own precompiled cache name; DEPS_LIB_TYPE
# "dynamic" packs lib/*.so* and lib/cmake into the precompiled zip.
if(CFS_BUILD_PRECICE AND UNIX)
  set(CFS_BOOST_SHARED ON)
  # -a: bump the letter when only the shared build changed (the first shared
  # builds had a broken empty DT_RUNPATH from an unsurvivable $ORIGIN linkflag)
  set(DEPS_ID "shared-a")
  set(DEPS_LIB_TYPE "dynamic")
  set(_BOOST_LIB_SUFFIX ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
  set(CFS_BOOST_SHARED OFF)
  set(_BOOST_LIB_SUFFIX ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# The cfsdeps *_LIBRARY cache variables and the ExternalProject stamps are
# write-once per build directory. Toggling USE_PRECICE in an already configured
# directory (also the normal ccmake flow: first configure with defaults, then
# enable it) switches the boost strategy (static <-> shared), but the cached
# BOOST_LIBRARY and possibly already built libs silently stay the other variant;
# the build then dies late with "No rule to make target 'lib/libboost_atomic.a'".
# Detect the switch and self-heal: drop the stale cache entries and discard the
# affected dependency builds (stamps + libs) so they are rebuilt in the new
# variant. zlib/libxml2/precice change with the switch too (PIC / shipped
# libz.so) but keep identical lib paths, so boost is the detectable canary.
# The distinct precompiled zips (DEPS_ID) of both variants stay in the cfsdeps
# cache and are reused on a flip back.
if(DEFINED BOOST_LIBRARY)
  string(REPLACE "." "\\." _BOOST_SUFFIX_RE "${_BOOST_LIB_SUFFIX}")
  list(GET BOOST_LIBRARY 0 _BOOST_CACHED_FIRST)
  if(NOT _BOOST_CACHED_FIRST MATCHES "${_BOOST_SUFFIX_RE}$")
    message(STATUS "boost strategy switched to ${DEPS_LIB_TYPE} (USE_PRECICE toggled): "
      "boost/zlib/libxml2/precice will be rebuilt in the new variant")
    unset(BOOST_LIBRARY CACHE)
    unset(BOOST_INCLUDE_DIR CACHE)
    # discard stamps + extracted sources (src) and the installs of the affected
    # deps so their next build starts clean (stale objects with the wrong
    # PIC-ness would otherwise be reused, e.g. by the libxml2 autotools build).
    # NOT the whole cfsdeps/<name> dir: it holds configure-time generated
    # scripts (zlib's were already written - it is included before boost).
    # Sources re-extract from the tarballs in CFS_DEPS_CACHE_DIR, no re-download.
    foreach(_STALE_DEP boost zlib libxml2 precice)
      file(REMOVE_RECURSE
        "${CMAKE_BINARY_DIR}/cfsdeps/${_STALE_DEP}/src"
        "${CMAKE_BINARY_DIR}/cfsdeps/${_STALE_DEP}/install")
    endforeach()
    # remove stale libs of the other variant from <build>/lib (would otherwise
    # linger and even get shipped by redistributables.cmake)
    file(GLOB _BOOST_STALE_LIBS
      "${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_*"
      "${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libprecice*"
      "${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libz${CMAKE_SHARED_LIBRARY_SUFFIX}*")
    if(_BOOST_STALE_LIBS)
      file(REMOVE ${_BOOST_STALE_LIBS})
    endif()
  endif()
endif()

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# generate BOOST_LIBRARY=PACKAGE_LIBRARY with os specific list of libs (static or
# shared, see CFS_BOOST_SHARED above). Always lib prefix, which is what b2 produces
# even on Windows (like set_package_library_list_lib_prefix, but suffix-aware).
foreach(_BOOST_LIB boost_atomic boost_iostreams boost_program_options boost_thread boost_chrono boost_log boost_regex boost_unit_test_framework boost_date_time boost_log_setup boost_serialization boost_filesystem)
  list(APPEND PACKAGE_LIBRARY ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/lib${_BOOST_LIB}${_BOOST_LIB_SUFFIX})
endforeach()

# to not always need to have all libs from BOOST_LIBRARY - no need for cache variables
set(BOOST_SERIALIZATION_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_serialization${_BOOST_LIB_SUFFIX})
set(BOOST_THREAD_LIB ${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libboost_thread${_BOOST_LIB_SUFFIX})

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
  set(WITHOUT --without-python --without-graph_parallel --without-mpi --without-container --without-context --without-math --without-contract --without-coroutine --without-graph --without-stacktrace --without-fiber --without-wave)

  if(CFS_BOOST_SHARED)
    # keep the boost CMake config (preCICE finds boost via find_package(Boost CONFIG);
    # FindBoost module mode is gone in cmake >= 3.30).
    # Deliberately NO -rpath,$ORIGIN linkflags on Linux: a literal $ORIGIN does not
    # survive b2's sh-executed link commands (it ends up as an EMPTY runpath entry,
    # which means CWD-search and, worse, a present DT_RUNPATH disables DT_RPATH
    # inheritance from the executable). The cfs executables carry a classic,
    # inherited DT_RPATH instead (see --disable-new-dtags below), through which
    # libboost_iostreams.so finds the shipped libz.so.
    if(APPLE)
      # macOS has no rpath inheritance - give the dylibs their own loader path
      set(_B2_LINK link=shared "linkflags=-Wl,-rpath,@loader_path")
    else()
      set(_B2_LINK link=shared)
    endif()
  else()
    # the static libs need no cmake config files (cfs links them via BOOST_LIBRARY)
    set(_B2_LINK link=static --no-cmake-config)
  endif()

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
    BUILD_COMMAND ./b2 --user-config=user-config.jam ${WITHOUT} --layout=system --prefix=${DEPS_INSTALL} ${TARGET} ${DEFINE} ${B2_ARGS} ${_B2_LINK} address-model=64 threading=multi runtime-link=shared variant=${_BUILD} debug-symbols=${_DS} install
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

# Boost.Log requires this define when linking the shared lib (the other compiled
# boost libs behave the same either way on gcc/clang). Deliberately NOT
# BOOST_ALL_DYN_LINK: that would also switch Boost.Test conventions, but the
# unittests use the header-only boost/test/included variant.
if(CFS_BOOST_SHARED)
  add_compile_definitions(BOOST_LOG_DYN_LINK)

  # With shared deps the executables must reliably load THEIR OWN libs from
  # lib/ (libprecice.so, libboost_*.so, libz.so): emit classic DT_RPATH instead
  # of DT_RUNPATH. DT_RPATH ranks BEFORE LD_LIBRARY_PATH and is inherited for
  # the whole loaded tree, so
  #  * a developer environment with LD_LIBRARY_PATH=<own precice install>
  #    cannot hijack libprecice.so.3 with a foreign build (linked against a
  #    different boost -> two boost versions in one process -> mixed-ABI heap
  #    corruption, seen as malloc abort at exit), and
  #  * transitive deps (libboost_iostreams.so -> libz.so.1) resolve from the
  #    shipped lib/ although the boost libs carry no rpath of their own.
  # Applies to build tree and installed package alike ($ORIGIN/../lib, see the
  # CMAKE_BUILD_RPATH block in the top-level CMakeLists.txt).
  if(NOT APPLE)
    add_link_options("LINKER:--disable-new-dtags")
  endif()
endif()
