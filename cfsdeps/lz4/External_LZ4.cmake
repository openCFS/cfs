# LZ4 is lossless compression algorithm
# https://github.com/lz4/lz4
# dependency for flann > 1.9.2
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "lz4")
set(PACKAGE_VER "1.10.0") # to request from github as use in the precompiled zip archive
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "dead9f5f1966d9ae56e1e32761e4e675")
set(PACKAGE_MIRRORS "https://github.com/lz4/lz4/archive/refs/tags/${PACKAGE_FILE}")  
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.


# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a C++ compiler
use_c_and_fortran(ON OFF)

set_precompiled_pckg_file()

set_package_library_default()

set_standard_variables() 

# Don't use use install_manifest.txt
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project and compiler flags
set_deps_args_default(ON) 
list(APPEND DEPS_ARGS "-DBUILD_SHARED_LIBS=OFF")
list(APPEND DEPS_ARGS "-DBUILD_STATIC_LIBS=ON")
list(APPEND DEPS_ARGS "-DLZ4_BUILD_CLI=OFF")

# the expection is here, that CMakeLists.txt is not in root (there is a Makefile)
set(DEPS_SUBDIR "build/cmake")

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

assert_unset(POSTINSTALL_SCRIPT)

# filter from DEPS_INSTALL, zip and copy to binary dir 
generate_packing_script_install_dir()

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
