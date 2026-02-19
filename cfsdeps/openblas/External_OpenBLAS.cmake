# OpenBLAS is a free and efficient and parallel BLAS/ LAPACK implementation.
# It is a successor of the stalled GOTO BLAS
# http://www.openblas.net/

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "openblas")
set(PACKAGE_VER "0.3.31")
set(OPENBLAS_VER ${PACKAGE_VER})
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "05050271d9196f65bc4ac3a89c6a3b05")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# this is the original version. There are variants on gitlab 
set(PACKAGE_MIRRORS "https://github.com/xianyi/OpenBLAS/archive/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# pure C
use_c_and_fortran(ON ON)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# just sets PACKAGE_LIBRARY=libopenblas.a
set_package_library_default()

# There is no separation between BLAS and LAPACK, both are combined. While BLAS is optimized
# only parts of LAPACK are optimized. OpenBLAS needs be build with LAPACK for openCFS 

# all possible BLAS implementations for openCFS need to set BLAS_LIBRARY. 
# LAPACK_LIBLIBRARY does not neet to be set, as we have the content in BLAS_LIBRARY
set(BLAS_LIBRARY ${PACKAGE_LIBRARY}) # works without additional -lpthread;-lm

set_standard_variables() # sets cached OPENBLAS_LIBRARY

# we don't want executabls, so go for install_dir
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) #  set compiler flags, they still add their own stuff 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_TESTING:BOOL=OFF) # not really so much stuff to set

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# no patch needed
assert_unset(PATCHES_SCRIPT)

# the forked metis insists onb bulding executables to bin. Sort them out automatically
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
  # standard cmake build without patch
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
