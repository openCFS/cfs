# Ginkgo is a high-performance linear algebra library for manycore systems, with a focus on solution of sparse linear systems. 
# Ginko has iterative solvers including AMG preconditioner
# https://ginkgo-project.github.io
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "ginkgo")
set(PACKAGE_VER "1.9.0")
set(ZLIB_VER ${PACKAGE_VER}) # required by boost 
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "461e21d1f7bdeacdf371af47b9c1b103")
set(DEPS_VER "") # set to "-a", "-b", when dependency changed with same PACKAGE_VER. Reset to "" with new PACKAGE_VER.

set(PACKAGE_MIRRORS "https://github.com/ginkgo-project/ginkgo/archive/refs/tags/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

if(WIN32 AND CMAKE_CXX_COMPILER_ID MATCHES "MSVC" AND USE_OPENMP)
  # msvc has even with /openmp:llvm not all necessary OpenMP features
  force_c_cxx_compiler("icx") # will set FORCE_C_CXX and do NOT set the default compiler options!
  assert_set(FORCE_C_CXX)
else()
  use_c_and_fortran(ON OFF)
endif()  

if(USE_OPENMP)
  set(DEPS_ID "OPENMP")
else()
  set(DEPS_ID "NO-OPENMP")
endif()

if(USE_CUDA)
  set(DEPS_ID "${DEPS_ID}_CUDA")
endif()

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# all executors (reference (=serial), omp, cuda,hip, dpcpp (sycl?!)) have empty default content in case
set_package_library_list("ginkgo;ginkgo_omp;ginkgo_cuda;ginkgo_reference;ginkgo_hip;ginkgo_dpcpp;ginkgo_device")

# set hidden cache variables *_LIBRARY = PACKAGE_LIBRARY, *_INCLUDE and some defaults
set_standard_variables()

# we don't need the share/man stuff and the dynamic lib
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARG with defaults for a cmake project
set_deps_args_default(ON) # set compiler flags 
# add the specific settings for the packge which comes in cmake style

# DGINKGO_BUILD_REFERENCE is the serial executor, we only use, when there is no OPENMP
invert(_NOT_OMP USE_OPENMP)

set(DEPS_ARGS
  ${DEPS_ARGS}
  -DBUILD_GMOCK=OFF
  -DBUILD_SHARED_LIBS=OFF
  -DBUILD_TESTING=OFF
  -DGINKGO_BUILD_BENCHMARKS=OFF
  -DGINKGO_BUILD_CUDA=${USE_CUDA}
  -DGINKGO_BUILD_EXAMPLES=OFF
  -DGINKGO_BUILD_REFERENCE=${_NOT_OMP}
  -DGINKGO_BUILD_OMP=${USE_OPENMP}
  -DGINKGO_BUILD_HIP=OFF
  -DGINKGO_BUILD_MPI=OFF
  -DGINKGO_BUILD_PAPI_SDE=OFF
  -DGINKGO_BUILD_SYCL=OFF
  -DGINKGO_BUILD_TESTS=OFF
  -DGINKGO_CONFIG_LOG_DETAILED=OFF
  -DGINKGO_ENABLE_HALF=OFF
  -DGINKGO_EXTENSION_KOKKOS_CHECK_=OFF
  -DGINKGO_VERBOSE_LEVEL=1
  -DGINKGO_WITH_CCACHE=OFF
  -DINSTALL_GTEST=OFF
  -DMETIS_INCLUDE_DIR=${METIS_INCLUDE_DIR}
  -DMETIS_LIBRARY=${METIS_LIBRARY})
  
if(FORCE_C_CXX)
  list(APPEND DEPS_ARGS -DCMAKE_CXX_COMPILER=${FORCE_C_CXX})
endif()  
    
# --- it follows generic final block for cmake packages w/o patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}" )

# Generate ${PACKAGE_NAME}-patch.cmake we use for our external project
generate_patches_script() # sets PATCHES_SCRIPT

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
  # standard cmake build with patch
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

add_dependencies(ginkgo metis nlohmann_json)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})