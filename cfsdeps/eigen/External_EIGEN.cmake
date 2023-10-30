# Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
# https://eigen.tuxfamily.org/index.php?title=Main_Page#Documentation
#
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "eigen")
set(PACKAGE_VER "3.4.0")
set(PACKAGE_FILE "eigen-${PACKAGE_VER}.tar.bz2")
set(PACKAGE_MD5 "132dde48fe2b563211675626d29f1707")
set(DEPS_VER "-a") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# the mirrors can point to arbitrary file names. 
set(PACKAGE_MIRRORS "https://gitlab.com/libeigen/eigen/-/archive/${PACKAGE_VER}/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we do not compile as eigen is header only
use_c_and_fortran(OFF OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# cmake package building directly to cfs build. precompiled package via manifest
set(DEPS_PREFIX  "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
set(DEPS_SOURCE  "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

# the clean-<package> target deletes everything to allow a clean make <package>
add_custom_target(clean-${PACKAGE_NAME} cmake -E remove_directory ${DEPS_PREFIX}
COMMAND cmake -E remove ${PRECOMPILED_PCKG_FILE}
COMMENT "delete cfsdeps/${PACKAGE_NAME} and precompiled")

# this is the standard target for cmake projects. The files to package come from the install_manifest.txt
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON)

# add flag for using only MPL2 or weaker-licenced parts of eigen
set(EIGEN_CXX_FLAGS "${CFSDEPS_CXX_FLAGS} -DEIGEN_MPL2_ONLY")

# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS ${DEPS_ARGS} -DBUILD_TESTING:BOOL=OFF -DINCLUDE_INSTALL_DIR:PATH=${DEPS_INSTALL}/include/ -DCMAKE_CXX_FLAGS:STRING=${EIGEN_CXX_FLAGS})

# --- it follows generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(PATCHES_SCRIPT)

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
