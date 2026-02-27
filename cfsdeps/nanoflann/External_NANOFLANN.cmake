# nanoflann: a C++11 header only KD-Trees replacement for FLANN for efficient neighbor search
# https://github.com/jlblancoc/nanoflann
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "nanoflann")
set(PACKAGE_VER "1.9.0")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "ffd8a26b95520011340b8bb696e78194")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.
# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://github.com/jlblancoc/nanoflann/archive/refs/tags/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

set(NANOFLANN_VER ${PACKAGE_VER}) # there is no nice version in nanoflann.hpp

 # we do not compile as nanoflann is header only
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# cmake package building directly to cfs build. precompiled package via manifest
set(DEPS_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
set(DEPS_SOURCE  "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

# the clean-<package> target deletes everything to allow a clean make <package>
add_clean_target()

# Don't use use install_manifest.txt
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project 
set_deps_args_default(ON)
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DNANOFLANN_BUILD_EXAMPLES=OFF
  -DNANOFLANN_BUILD_TESTS=OFF)

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

# filter from DEPS_INSTALL, zip and copy to binary dir 
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

# dump_dependency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  create_external_cmake()  

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
