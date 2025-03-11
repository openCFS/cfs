# SuperLU - a direct solver
# https://portal.nersc.gov/project/sparse/superlu/
# this is the serial version, it is recommended to switch to the OpenMP version SuperLU-MT
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "superlu")
set(PACKAGE_VER "6.0.1")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "d15c61705f4ddf0777731d3f388e287f")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# this is the original version. There are variants on gitlab 
set(PACKAGE_MIRRORS "https://github.com/xiaoyeli/superlu/archive/refs/tags/${PACKAGE_FILE}")
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
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -Denable_complex16:BOOL=ON
  -Denable_double:BOOL=ON
  -Denable_complex:BOOL=OFF 
  -Denable_doc:BOOL=OFF
  -Denable_examples:BOOL=OFF
  -Denable_fortran:BOOL=OFF
  -Denable_internal_blaslib:BOOL=OFF
  -Denable_matlabmex:BOOL=OFF
  -Denable_single:BOOL=OFF
  -Denable_tests:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF,
  -DTPL_ENABLE_METISLIB:BOOL=OFF)
  
if(UNIX AND USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  list(APPEND DEPS_ARGS -DTPL_BLAS_LIBRARIES=${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libopenblas.a)
elseif(UNIX AND USE_BLAS_LAPACK STREQUAL "MKL") # we assume properly set up mkl for Windows
  list(APPEND DEPS_ARGS -DTPL_BLAS_LIBRARIES=${MKL_LIB_DIR}/libmkl_intel_lp64.a) 
endif()  

# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# no pre-compile patch needed
assert_unset(PATCHES_SCRIPT)

# we need a post-install patch, so better no manifest
generate_packing_script_install_dir()

generate_postinstall_script()

#dump_depencency_variables()

# do we want to use precompiled and do we already have the package?
if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # copy files from cache
  create_external_unpack_precompiled()

# if not, build newly and possibly pack the stuff
else()
  # standard cmake build without patch
  create_external_cmake()  
#  add_postinstall_step()
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