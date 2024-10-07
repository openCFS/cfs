# topology optimization optimiers MMA implementations by Jérémie Dumas
# https://github.com/jdumas/mma contains MMASolver and GCMMASolver based on https://github.com/topopt/TopOpt_in_PETSc
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "dumas")
set(PACKAGE_VER "aa51333") # to request from github as use in the precompiled zip archive
set(DUMAS_VER "0.99") # for def_use_dumas.hh.in to be used in Dependencies.cc
set(PACKAGE_FILE "${PACKAGE_VER}.zip")
set(PACKAGE_MD5 "22e1eebd6bb65dad2e0694c859050635")
set(PACKAGE_MIRRORS "https://github.com/jdumas/mma/archive/${PACKAGE_FILE}")  
set(DEPS_VER "-a") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# we cannot link a parallel compiled suitesparse with debug without openmp
# we need to set DEPS_ID before calling set_precompiled_pckg_file()
if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO_OPENMP")
endif()

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a C++ compiler
use_c_and_fortran(ON OFF)

set_precompiled_pckg_file()

# libdumas.a
set_package_library_default()

set_standard_variables() 

# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${CMAKE_BINARY_DIR}")

# set DEPS_ARG with defaults for a cmake project and compiler flags
set_deps_args_default(ON) 
list(APPEND DEPS_ARGS "-DUSE_OPENMP=${USE_OPENMP}") # our own CMakeLists.txt trigges by this OpenMP

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# use our own CMakeLists.txt 
generate_patches_script()

assert_unset(POSTINSTALL_SCRIPT)

# generate package creation script. We get the files from an install_manifest.txt
generate_packing_script_manifest()

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
  endif()
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
