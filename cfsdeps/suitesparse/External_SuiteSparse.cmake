# SuiteSparse contains the directl solvers CHOLMOD (extremely fast!) and UMFPACK
# https://github.com/DrTimothyAldenDavis/SuiteSparse
# we configure for GLP. Without choldmod is quite slow for 3D
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "suitesparse")
# SuiteSparse >= 6 is CMake based, but has not root CMakeLists.txt which we add ourselves
# see also https://github.com/Fabian188/SuiteSparse-root-cmake
# we compile this GPL with all available modules
set(PACKAGE_VER "7.0.1")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "d31bbe2a26dced338b23e71f7c9b541a")
set(PACKAGE_MIRRORS "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/${PACKAGE_FILE}")  
set(DEPS_VER "-a") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a C compiler
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path

# we cannot link a parallel compiled suitesparse with debug without openmp
if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO-OPENMP")
endif()

set_precompiled_pckg_file()

set_static_cache_lib("AMD_LIBRARY" "amd" "AMD lib from SuiteSparse")

# Windows has non-standard _static.lib ending
assert_unset(PACKAGE_LIBRARY)
set(LIBS "umfpack;cholmod;camd;ccolamd;colamd;amd;suitesparseconfig")
foreach(ITEM ${LIBS})
  if(UNIX)
    list(APPEND PACKAGE_LIBRARY "${CMAKE_BINARY_DIR}/lib/lib${ITEM}.a")
  else()
    list(APPEND PACKAGE_LIBRARY "${CMAKE_BINARY_DIR}/lib/${ITEM}_static.lib")  
  endif()  
endforeach()   

# creates SUITESPARSE_LIBARAY as CACHE variable, hence it will not be overwritten once in cache!
set_standard_variables() 

# we have no trustworthy install_manifest.txt, hence use install dir
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 

# we potentiall build the suitesparse subprojects in parallel but we need to process them sequentially
set(DEPS_BUILD_THREADS 1)

set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_STATIC:BOOL=ON
  -DUSE_OPENMP:BOOL=${USE_OPENMP}
  -DALLOW_64BIT_BLAS:BOOL=ON
  -DALLOW_GPL_EXTENSIONS:BOOL=ON 
  -DBUILD_THREADS=${CFS_DEPS_BUILD_THREADS}) # the global number of build threads (default is system threads)

if(UNIX AND USE_BLAS_LAPACK STREQUAL "OPENBLAS")
  list(APPEND DEPS_ARGS -DSUGGEST_BLAS_LIBRARIES=${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/libopenblas.a)
elseif(UNIX AND USE_BLAS_LAPACK STREQUAL "MKL") # we assume properly set up mkl for Windows
  list(APPEND DEPS_ARGS -DSUGGEST_BLAS_LIBRARIES=${MKL_LIB_DIR}/libmkl_intel_lp64.a) 
endif()  

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# copies your CMakeLists.txt and suitesparse_install.cmake for SuiteSparse 7 as long this is not upstream
generate_patches_script()

assert_unset(POSTINSTALL_SCRIPT)

# no manifest 
generate_packing_script_install_dir()

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
    # this will dump the unnecessary shared variants, but the precompiled package is clean  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# copy the appropriate special license in the license folder
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/LICENSE.binary.GPL" DESTINATION "${CFS_BINARY_DIR}/license/")

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
