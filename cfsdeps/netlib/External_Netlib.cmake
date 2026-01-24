# from the possible BLAS/LAPACK implementations, Netlib is the reference implementation
# but also very slow as it is not optimized.
# The lapack package builds blas and lapack
# http://www.netlib.org/

# make sure not to uninetendently use another packages settings. Supports assert_set() checks. Is mandatory!
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "netlib")
set(PACKAGE_VER "3.12.1")
set(NETLIB_VER  ${PACKAGE_VER}) # for Dependencies.cc
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "2f069617e16b42f5eddcfee85768f204")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# this is the original version. There are variants on gitlab 
set(PACKAGE_MIRRORS "https://github.com/Reference-LAPACK/lapack/archive/refs/tags/${PACKAGE_FILE}")

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# pure C
use_c_and_fortran(OFF ON)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# we have BLAS and LAPACK as separate libs in netlib
# BLAS_LIBRARY (and LAPACK_LIBRARY) need to be set for cfs by every BLAS/LAPACK provider
# we assume CFS_FORTRAN_LIBS to be set target libs whenever we use BLAS/LAPACK
set(BLAS_LIBRARY ${CMAKE_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}blas${CMAKE_STATIC_LIBRARY_SUFFIX})
set(LAPACK_LIBRARY ${CMAKE_BINARY_DIR}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lapack${CMAKE_STATIC_LIBRARY_SUFFIX})
set(PACKAGE_LIBRARY ${LAPACK_LIBRARY};${BLAS_LIBRARY})

set_standard_variables() # 

# we don't want executabls, so go for install_dir
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) #  set compiler flags, they still add their own stuff 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_TESTING:BOOL=OFF
  -DBUILD_COMPLEX:BOOL=OFF # single precision complex
  -DBUILD_SINGLE:BOOL=OFF # single precision
  -DTEST_FORTRAN_COMPILER:BOOL=OFF) 

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
