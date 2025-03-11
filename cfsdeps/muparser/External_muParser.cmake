# muParser - a fast math parser library
# http://muparser.sourceforge.net/

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "muparser")
set(PACKAGE_VER "2.2.6.1")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "410d29b4c58d1cdc2fc9ed1c1c7f67fe") # 2.2.6.1
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# this is the original version. There are variants on gitlab 
set(PACKAGE_MIRRORS "https://github.com/beltoforion/muparser/archive/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# pure C
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# we don't want executabls, so go for install_dir
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # use default compiler flags 
# add the specific settings for the packge which comes in cmake style

# whyever we really need on Windows the shared lib?!
set(MUPARSER_SHARED_LIBS OFF)
if(WIN32)
  set(MUPARSER_SHARED_LIBS ON)
endif()

set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_SHARED_LIBS:BOOL=${MUPARSER_SHARED_LIBS}
  -DBUILD_TESTING:BOOL=OFF 
  # muParser 2.2.6.1 has issues with openMP. On Mac it does not find omp.h, on Linux "Could NOT find OpenMP_CXX (missing: OpenMP_CXX_FLAGS OpenMP_CXX_LIB_NAMES)"
  -DENABLE_OPENMP:BOOL=OFF 
  -DENABLE_SAMPLES:BOOL=OFF )

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# no patch needed
assert_unset(PATCHES_SCRIPT)

# generate package creation script. We use manifest as for Win we also have to handle the .dll
generate_packing_script_manifest()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # standard cmake build without patch
  create_external_cmake()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif() 
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
