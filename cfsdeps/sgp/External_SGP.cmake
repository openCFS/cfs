# SGP (Sequqential Global Programming) is an open source optimizer for structural optimization
# Developed by the Chair of Continuous Optimization, Friedrich-Alexander-Universität Erlangen-Nürberg
# https://gitlab.com/MOGLI_BM/sgp-cfs/-/tree/debug-3drot

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "sgp")
set(PACKAGE_VER "23.7.12")
set(SGP_VER ${PACKAGE_VER}) # for Dependencies.cc
# this is currently a branch from Ngoc, the next version might be already from the master (-> mirror)
set(PACKAGE_FILE "sgp-cfs-debug-3drot-${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "186e68b457f21ee786c9dcf7fa1a36f5")
set(DEPS_VER "-a") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# OpenMP
if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO_OPENMP")
endif()

# replace when we can download a tagged version from master by a next revision
set(PACKAGE_MIRRORS "${CFS_FAU_MIRROR}/sources/sgp/${PACKAGE_FILE}") 

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# determine paths of libraries and make it visible (and editable) via ccmake
set_package_library_default()
# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()
# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
assert_set(MKL_ROOT) # should be provided by MKLConfig.cmake of oneAPI
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DMKL_ROOT_DIR=${MKL_ROOT}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DUSE_OPENMP:BOOL=${USE_OPENMP})

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# no patch needed
assert_unset(PATCHES_SCRIPT)

# generate package creation script. We get the files from an install_manifest.txt
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
  # add external project step actually building an cmake package including a patch 
  # also genearate the patch script via generate_patches_script()
  create_external_cmake()  

  # new data just built: shall we pack and store as precompiled?
  if(${CFS_DEPS_PRECOMPILED})
    # add custom step to zip a precompiled package to the cache.
    add_external_storage_step()
  endif()  
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})