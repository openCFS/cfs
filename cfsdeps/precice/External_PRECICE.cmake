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

# Optional MPI / PETSc features (USE_PRECICE_MPI / USE_PRECICE_PETSC, decided in
# FindCFSDEPS.cmake; PETSc implies MPI). Encode them into DEPS_VER *before*
# set_precompiled_pckg_file() so feature variants use distinct precompiled caches
# instead of clobbering each other.
if(USE_PRECICE_MPI)
  set(DEPS_VER "${DEPS_VER}-mpi")
endif()
if(USE_PRECICE_PETSC)
  set(DEPS_VER "${DEPS_VER}-petsc")
endif()

set(PACKAGE_MIRRORS "https://github.com/precice/precice/archive/refs/tags/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

# C++ library: build with C (and no fortran)
use_c_and_fortran(ON OFF)

# sets PRECOMPILED_PCKG_FILE to the full precompiled name including path
set_precompiled_pckg_file()

# preCICE 3.x is ALWAYS a shared library (add_library(precice SHARED); BUILD_SHARED_LIBS
# is ignored), so the artifact is libprecice.so, not a static .a. Set the library
# variable manually to the shared lib (set_package_library_list would assume .a).
set(PACKAGE_LIBRARY "${CMAKE_BINARY_DIR}/${LIB_SUFFIX}/${CMAKE_SHARED_LIBRARY_PREFIX}precice${CMAKE_SHARED_LIBRARY_SUFFIX}")
# set_standard_variables() also sets DEPS_PREFIX / DEPS_SOURCE and the clean-precice target
# (PRECICE_LIBRARY = the .so above, PRECICE_INCLUDE_DIR = <build>/include).
set_standard_variables()

# install into the package's own prefix and copy to the cfs build dir (no manifest)
set(DEPS_INSTALL "${DEPS_PREFIX}/install")

# set DEPS_ARGS with defaults for a cmake project (install dirs, build type, compilers, flags)
set_deps_args_default(ON)

# preCICE-specific configuration:
#  * minimal feature set (no Python; MPI and PETSc are optional - see
#    USE_PRECICE_MPI / USE_PRECICE_PETSC) to keep the dependency surface small,
#  * find boost / eigen / libxml2 from the cfsdeps installs in the cfs build dir.
# preCICE is always shared, so we do NOT pass BUILD_SHARED_LIBS (it is ignored
# and only warns). The required deps are boost / eigen / libxml2 (+ Threads, and
# MPI / PETSc when the corresponding feature is enabled).
if(USE_PRECICE_MPI)
  set(_PRECICE_MPI ON)
else()
  set(_PRECICE_MPI OFF)
endif()
if(USE_PRECICE_PETSC)
  set(_PRECICE_PETSC ON)
else()
  set(_PRECICE_PETSC OFF)
endif()
set(DEPS_ARGS ${DEPS_ARGS}
  -DBUILD_TESTING:BOOL=OFF
  -DPRECICE_FEATURE_PETSC_MAPPING:BOOL=${_PRECICE_PETSC}
  -DPRECICE_FEATURE_MPI_COMMUNICATION:BOOL=${_PRECICE_MPI}
  -DPRECICE_FEATURE_PYTHON_ACTIONS:BOOL=OFF
  -DPRECICE_BINDINGS_C:BOOL=OFF
  -DPRECICE_BINDINGS_FORTRAN:BOOL=OFF
  -DPRECICE_BUILD_TOOLS:BOOL=OFF
  # CMAKE_INSTALL_PREFIX / _LIBDIR / build type / compilers come from set_deps_args_default().
  # eigen is header-only and always taken from the normal cfsdeps install.
  -DEigen3_ROOT:PATH=${CMAKE_BINARY_DIR})

# When PETSc mapping is on, point preCICE's FindPETSc at openCFS' PETSc:
#  * USE_PETSC -> the cfsdeps PETSc prefix install (full prefix with lib/petsc/conf,
#    populated whether freshly built or restored from the precompiled cache),
#  * otherwise rely on the system PETSC_DIR / PETSC_ARCH environment (pkg-config).
# PETSC_ARCH is empty for a prefix install.
if(USE_PRECICE_PETSC AND USE_PETSC)
  list(APPEND DEPS_ARGS
    -DPETSC_DIR:PATH=${CMAKE_BINARY_DIR}/cfsdeps/petsc/install
    -DPETSC_ARCH:STRING=)
endif()

# boost + libxml2 follow the detect-or-build model (see cfsdeps/boost + cfsdeps/libxml2):
# if a system/shared copy was detected, *_ROOT is empty and preCICE finds it itself;
# if we built a copy, point preCICE at its prefix (single-path, package-case *_ROOT,
# honored via CMP0074). preCICE links *shared* boost, so no Boost_USE_STATIC_LIBS.
if(PRECICE_BOOST_ROOT)
  list(APPEND DEPS_ARGS -DBoost_ROOT:PATH=${PRECICE_BOOST_ROOT})
endif()
if(LIBXML2_PIC_ROOT)
  list(APPEND DEPS_ARGS -DLibXml2_ROOT:PATH=${LIBXML2_PIC_ROOT})
endif()

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

# preCICE needs (header-only) eigen and, when we build them ourselves, the
# boost / libxml2 copies first. Under the detect-or-build model those targets
# only exist when a system copy was NOT found, so guard with TARGET checks.
add_dependencies(precice eigen)
if(TARGET boost-pic)
  add_dependencies(precice boost-pic)
endif()
if(TARGET libxml2-pic)
  add_dependencies(precice libxml2-pic)
endif()
# with PETSc mapping against the cfsdeps PETSc, build preCICE after petsc
if(USE_PRECICE_PETSC AND TARGET petsc)
  add_dependencies(precice petsc)
endif()

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
