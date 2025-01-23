# FLANN - Fast Library for Approximate Nearest Neighbors
# https://github.com/flann-lib/flann

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "flann")
# 1.9.2 has the issues that is requires liblz4.a (cmake in build cmake) and that it failes on Windows "Could NOT find PkgConfig (missing: PKG_CONFIG_EXECUTABLE)"
set(PACKAGE_VER "1.9.2")
set(PACKAGE_FILE "${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "9a1f10c0d890a9595f2f4312436af50f") # 1.9.2
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

#https://github.com/flann-lib/flann/pull/507/commits/e6afde06ae8071be96d9c15e9953dd71e5ebee91

if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO_OPENMP")
endif()

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/flann-lib/flann/archive/refs/tags/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# our lib has not the standard name but standard prefix and suffix
set_package_library_list("flann_cpp_s")

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for cmake projects but we don't want to use install_manifest.txt but pick the stuff
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_DOC:BOOL=OFF
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_C_BINDINGS:BOOL=OFF
  -DBUILD_MATLAB_BINDINGS:BOOL=OFF
  -DBUILD_PYTHON_BINDINGS:BOOL=OFF
  -DBUILD_TESTS:BOOL=OFF
  -DUSE_MPI:BOOL=OFF
  -DUSE_OPENMP:BOOL=${USE_OPENMP}
  -Dlz4_DIR:PATH=${CMAKE_BINARY_DIR}) # we add this feature with our patch

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# Generate ${PACKAGE_NAME}-patch.cmake we use for our external project
generate_patches_script() # sets PATCHES_SCRIPT

# generate package ceation script. We get the files from an install directory
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  create_external_cmake_patched()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  else()
    # without manifest (installs directly to binary dir) an without packing, we need to copy manually  
    add_install_dir_to_binary_step()  
  endif()  
endif()

add_dependencies(flann lz4)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
