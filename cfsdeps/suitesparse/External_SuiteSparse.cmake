# SuiteSparse contains the directl solvers CHOLMOD (extremely fast!) and UMFPACK
# https://github.com/DrTimothyAldenDavis/SuiteSparse
#
# SuiteSparse 6.0.1 is now CMake based, but has not root CMakeLists.txt which we add ourselves
# see also https://github.com/Fabian188/SuiteSparse-root-cmake
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "suitesparse")
set(PACKAGE_VER "6.0.1")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "3bb660ac217791c7e9fabac944c8ee07")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

set(PACKAGE_MIRRORS "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/${PACKAGE_FILE}") 
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a fortran compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
if(USE_SUITESPARSE_GPL)
  set(DEPS_ID "GPL")
else()
 set(DEPS_ID "NOGPL")
endif() 

# we cannot link a parallel compiled suitesparse with debug without openmp
if(USE_OPENMP)
  set(DEPS_ID "${DEPS_ID}-OPENMP")
else()
  set(DEPS_ID "${DEPS_ID}-NO_OPENMP")
endif()

set_precompiled_pckg_file()

set_static_cache_lib("AMD_LIBRARY" "amd" "AMD lib from SuiteSparse")

# generate PACKAGE_LIBARAY with os specific list of static libs
set_package_library_list("umfpack;cholmod;camd;ccolamd;colamd;amd;suitesparseconfig")
# creates SUITESPARSE_LIBARAY as CACHE variable, hence it will not be overwritten once in cache!
set_standard_variables() 

# we have no trustworthy install_manifest.txt, hence use install dir
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default() 
# add the specific settings for the packge which comes in cmake style
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_STATIC:BOOL=ON
  -DUSE_OPENMP:BOOL=${USE_OPENMP}
  -DALLOW_64BIT_BLAS:BOOL=ON
  -DALLOW_GPL_EXTENSIONS=${USE_SUITESPARSE_GPL} )
if(USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  set(DEPS_ARGS ${DEPS_ARGS} -DSUGGEST_BLAS_LIBRARIES=${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libopenblas.a) # we assue mkl to be used for Windows
elseif(USE_BLAS_LAPACK STREQUAL "MKL")
  set(DEPS_ARGS ${DEPS_ARGS} -DSUGGEST_BLAS_LIBRARIES=${MKL_LIB_DIR}/libmkl_intel_lp64.a) # add windows
endif()  
  
# --- it follows generic final block for cmake packages with a patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/"
     DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# copies your CMakeLists.txt
generate_patches_script()

# generate package ceation script. Somehow the install_manifest.txt fails for snopt. It is not unacked to lib. Possibly EOL issue?!
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

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
