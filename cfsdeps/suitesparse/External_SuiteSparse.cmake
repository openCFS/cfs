# SuiteSparse contains the direct solvers CHOLMOD (extremely fast!) and UMFPACK
# https://github.com/DrTimothyAldenDavis/SuiteSparse
# we configure for GLP. Without, cholmod is quite slow for 3D
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "suitesparse")
# SuiteSparse >= 7 has good CMake support, we use it as-is
# we compile the required (sub-)projects including GPL licensed code paths
set(PACKAGE_VER "7.13.0") #
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "704e33348f0ef892c806aa1f480f3ba3") # 7.13.0
set(PACKAGE_MIRRORS "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/${PACKAGE_FILE}")  
set(DEPS_VER "-b") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

 # we only have a C compiler
use_c_and_fortran(ON OFF)

# we cannot link a parallel compiled suitesparse with debug without openmp
if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO-OPENMP")
endif()

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# Windows has non-standard _static.lib ending
assert_unset(PACKAGE_LIBRARY)
set(LIBS "umfpack;cholmod;camd;ccolamd;colamd;amd;suitesparse_config")
foreach(ITEM ${LIBS})
  string(REPLACE "_" "" ITEM "${ITEM}") # remove the underscore in suitesparse_config because the lib is called differently than the project
  if(UNIX)
    list(APPEND PACKAGE_LIBRARY "${CMAKE_BINARY_DIR}/lib/lib${ITEM}.a")
  else()
    list(APPEND PACKAGE_LIBRARY "${CMAKE_BINARY_DIR}/lib/${ITEM}_static.lib")  
  endif()
endforeach()

# creates SUITESPARSE_LIBRAY as CACHE variable, hence it will not be overwritten once in cache!
set_standard_variables() 

# we have no trustworthy install_manifest.txt, hence use install dir
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 

# we potentially build the suitesparse subprojects in parallel but we need to process them sequentially
set(DEPS_BUILD_THREADS 1)

# switch off shared-libs, demos and testing so suitesparse compiles without the need to link BLAS/LAPACK
# enabling umfpack enables all other LIBS since they are dependencies. This avoids weird CMake issues with quotes and semicolons
set(DEPS_ARGS
  ${DEPS_ARGS}
  -DSUITESPARSE_ENABLE_PROJECTS=umfpack
  -DBUILD_STATIC_LIBS:BOOL=ON
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DSUITESPARSE_USE_OPENMP:BOOL=${USE_OPENMP}
  -DSUITESPARSE_DEMOS=OFF
  -DBUILD_TESTING=OFF
  -DSUITESPARSE_USE_64BIT_BLAS:BOOL=OFF
  -DCHOLMOD_GPL:BOOL=ON 
  -DCHOLMOD_USE_CUDA:BOOL=OFF
  )

# make sure SUITESPARSE does not try to find BLAS/LAPACK
# this is not required as we link it in CFS
list(APPEND DEPS_ARGS -DBLA_VENDOR=Generic -DBLAS_LIBRARIES=/dummy -DLAPACK_LIBRARIES=/dummy)

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

assert_unset(POSTINSTALL_SCRIPT)

# no manifest 
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
    # this will dump the unnecessary shared variants, but the precompiled package is clean  
    add_install_dir_to_binary_step()  
  endif()  
endif()

# copy the appropriate special license in the license folder
file(COPY "${CFS_SOURCE_DIR}/cfsdeps/LICENSE.binary.GPL" DESTINATION "${CFS_BINARY_DIR}/license/")

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
