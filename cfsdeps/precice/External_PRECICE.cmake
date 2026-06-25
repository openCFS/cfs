# preCICE - a coupling library for partitioned multi-physics simulations
# (fluid-structure interaction, conjugate heat transfer, ...).
# https://precice.org
#
# Built from source in the cfsdeps stage against cfs' own boost / eigen / libxml2
# so the boost ABI is consistent with the rest of openCFS. Consumed by
# source/Utils/preciceAdapter via the PRECICE_LIBRARY / PRECICE_INCLUDE_DIR
# cache variables (the cfsdeps convention - NOT find_package, which cannot run at
# configure time before the dependency is built).
clear_depencency_variables()

# set mandatory variables for the macros in DependencyTools.cmake.
set(PACKAGE_NAME "precice")
set(PACKAGE_VER "3.4.1")
set(PACKAGE_FILE "v${PACKAGE_VER}.tar.gz")
set(PACKAGE_MD5 "22001a546429de058adff723e9c0870e")
set(DEPS_VER "") # set to "-a", "-b", ... when only the build (not PACKAGE_VER) changed

set(PACKAGE_MIRRORS "https://github.com/precice/precice/archive/refs/tags/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# C++ library: build with C (and no fortran)
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# PRECICE_LIBRARY = <build>/lib/libprecice.a , PRECICE_INCLUDE_DIR = <build>/include.
# set_standard_variables() also sets DEPS_PREFIX / DEPS_SOURCE and the clean-precice target.
set_package_library_list("precice")
set_standard_variables()

# install into the package's own prefix and copy to the cfs build dir (no manifest)
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARGS with defaults for a cmake project (install dirs, build type, compilers, flags)
set_deps_args_default(ON)

# preCICE-specific configuration:
#  * static lib, minimal feature set (no PETSc / MPI / Python) to keep the
#    dependency surface small and avoid extra runtime deps in the testsuite,
#  * find boost / eigen / libxml2 from the cfsdeps installs in the cfs build dir.
set(DEPS_ARGS ${DEPS_ARGS}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DPRECICE_FEATURE_PETSC_MAPPING:BOOL=OFF
  -DPRECICE_FEATURE_MPI_COMMUNICATION:BOOL=OFF
  -DPRECICE_FEATURE_PYTHON_ACTIONS:BOOL=OFF
  -DPRECICE_BINDINGS_PYTHON:BOOL=OFF
  -DPRECICE_BINDINGS_C:BOOL=OFF
  -DPRECICE_BINDINGS_FORTRAN:BOOL=OFF
  # CMAKE_INSTALL_PREFIX / _LIBDIR / build type / compilers come from
  # set_deps_args_default() above.
  # let preCICE find the cfsdeps-built boost / eigen / libxml2 (their CMake
  # config / headers are copied into the cfs build dir by the cfsdeps macros)
  -DCMAKE_PREFIX_PATH:PATH=${CMAKE_BINARY_DIR}
  -DBOOST_ROOT:PATH=${CMAKE_BINARY_DIR})

# --- generic final block for cmake packages with no patch and no postinstall ---

# copy "static" license as we configure this dependency. Check if license is still valid!
file(COPY "${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/license/" DESTINATION "${CMAKE_BINARY_DIR}/license/${PACKAGE_NAME}")

# preCICE needs no patch
assert_unset(PATCHES_SCRIPT)

# filter from DEPS_INSTALL, zip and copy to binary dir
generate_packing_script_install_dir()

# we have no postinstall, so don't call generate_postinstall_script()
assert_unset(POSTINSTALL_SCRIPT)

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
    # without manifest (installs directly to binary dir) and without packing, copy manually
    add_install_dir_to_binary_step()
  endif()
endif()

# preCICE is built against cfs' boost / eigen / libxml2 -> build those first
add_dependencies(precice boost eigen libxml2)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
