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
# pack lib/*.so* (versioned + symlink) and lib/cmake into the precompiled zip
set(DEPS_LIB_TYPE "dynamic")
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

# With MPI communication on, preCICE runs its own find_package(MPI). The cfsdeps
# build uses the plain (non-mpi) compilers, and the CI _max image ships OpenMPI
# but keeps it OFF the default PATH (so other deps do not silently build MPI
# variants - see cfsdeps:gcc9max in .gitlab-ci.yml). preCICE's find_package(MPI)
# then finds nothing and dies with "Could NOT find MPI (missing: MPI_CXX_FOUND)".
# Locate the MPI wrappers the same way OpenFOAM does (find_program with the RHEL/
# Fedora openmpi HINTS, falling back to PATH) and hand them to preCICE only, so
# the "OpenMPI off PATH" invariant stays intact for every other dependency.
if(USE_PRECICE_MPI)
  find_program(PRECICE_MPICC  NAMES mpicc          HINTS /usr/lib64/openmpi/bin /usr/lib/openmpi/bin)
  find_program(PRECICE_MPICXX NAMES mpicxx mpic++  HINTS /usr/lib64/openmpi/bin /usr/lib/openmpi/bin)
  mark_as_advanced(PRECICE_MPICC PRECICE_MPICXX)
  if(PRECICE_MPICC AND PRECICE_MPICXX)
    message(STATUS "preCICE: using MPI wrappers ${PRECICE_MPICXX} / ${PRECICE_MPICC}")
    list(APPEND DEPS_ARGS
      -DMPI_C_COMPILER:FILEPATH=${PRECICE_MPICC}
      -DMPI_CXX_COMPILER:FILEPATH=${PRECICE_MPICXX})
  else()
    message(WARNING
      "preCICE: USE_PRECICE_MPI is ON but no MPI wrapper (mpicc/mpicxx) was found "
      "in the default PATH or /usr/lib64/openmpi/bin - preCICE's find_package(MPI) "
      "will likely fail. Install openmpi-devel / libopenmpi-dev, or put the MPI bin "
      "dir on PATH before configuring.")
  endif()
endif()

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

# boost + libxml2 are the single cfsdeps builds; their strategy switches when
# preCICE is built (shared boost with CMake config / static -fPIC libxml2, see
# their External_*.cmake). Both install into the cfs build dir, so point preCICE
# there (package-case *_ROOT, honored via CMP0074). Boost is found via its CMake
# config in <build>/lib/cmake, libxml2 via the FindLibXml2 module.
list(APPEND DEPS_ARGS
  -DBoost_ROOT:PATH=${CMAKE_BINARY_DIR}
  -DLibXml2_ROOT:PATH=${CMAKE_BINARY_DIR})

# libprecice.so must find the shipped libboost_*.so next to itself in lib/ on its
# own: RUNPATH (the modern DT_RPATH replacement the linker emits for the cfs
# binary) is not inherited by dependencies, so the cfs rpath does not help here.
if(APPLE)
  list(APPEND DEPS_ARGS -DCMAKE_INSTALL_RPATH:STRING=@loader_path)
elseif(UNIX)
  list(APPEND DEPS_ARGS -DCMAKE_INSTALL_RPATH:STRING=$ORIGIN)
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

# preCICE builds against the cfsdeps installs: (header-only) eigen, the shared
# boost and the static -fPIC libxml2 (both always cfsdeps targets when
# CFS_BUILD_PRECICE, see FindCFSDEPS.cmake include conditions).
add_dependencies(precice eigen boost libxml2)
# with PETSc mapping against the cfsdeps PETSc, build preCICE after petsc
if(USE_PRECICE_PETSC AND TARGET petsc)
  add_dependencies(precice petsc)
endif()

# pkg-config file for consumers outside the cfs cmake world (e.g. the preCICE
# OpenFOAM adapter, see cfsdeps/openfoam-adapter). Points at the cfs build dir,
# where the precice install lands both freshly built and precompiled-restored.
configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/libprecice.pc.in"
               "${CMAKE_BINARY_DIR}/lib/pkgconfig/libprecice.pc" @ONLY)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
