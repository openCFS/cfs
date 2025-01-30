# preCICE OpenFOAM adapter - loaded by OpenFOAM solvers as a function object
# ('libs (preciceAdapterFunctionObject);' in a controlDict), turning them into
# preCICE participants for coupled openCFS<->OpenFOAM simulations.
#
# Build/test-only dependency like cfsdeps/openfoam (GPL, runs as external program,
# nothing shipped in the cfs package). Build order per precice.org:
# dependencies -> preCICE -> OpenFOAM -> OpenFOAM adapter (see add_dependencies below).
clear_depencency_variables()

set(PACKAGE_NAME "openfoam-adapter")

# ---- adapter source - keep easily editable ---------------------------------
# pinned to the openCFS fork/branch with the speaker implementations
set(OPENFOAM_ADAPTER_GIT "https://github.com/DM-TUG/openfoam-adapter.git")
set(OPENFOAM_ADAPTER_TAG "speaker_implementations") # branch, tag or commit
# -----------------------------------------------------------------------------

# OPENFOAM_DIR and OPENFOAM_MPI_BINDIR (may be empty) come from
# cfsdeps/openfoam/External_OpenFOAM.cmake, included right before this file
assert_set(OPENFOAM_DIR)

set(DEPS_PREFIX "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
set(DEPS_SOURCE "${DEPS_PREFIX}/src/${PACKAGE_NAME}")

# The adapter's Allwmake discovers preCICE via pkg-config: External_PRECICE.cmake
# generates <build>/lib/pkgconfig/libprecice.pc pointing at the cfs build dir
# (with an rpath, so the adapter .so finds libprecice.so.3 at runtime without any
# environment). The adapter lib itself is installed into OpenFOAM's own lib dir
# (FOAM_LIBBIN via PRECICE_OPENFOAM_TARGET_DIR), so solvers find it with the
# plain OpenFOAM environment (source <build>/openfoam-env.sh).
set(OPENFOAM_ADAPTER_BUILD_SCRIPT "${DEPS_PREFIX}/openfoam-adapter-build.sh")
configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/openfoam-adapter-build.sh.in"
               "${OPENFOAM_ADAPTER_BUILD_SCRIPT}" @ONLY)

# a moving branch has no stable hash -> fetched via git (shallow), not the
# tarball+md5 route of the other cfsdeps
ExternalProject_Add(openfoam-adapter
  PREFIX "${DEPS_PREFIX}"
  SOURCE_DIR "${DEPS_SOURCE}"
  GIT_REPOSITORY "${OPENFOAM_ADAPTER_GIT}"
  GIT_TAG "${OPENFOAM_ADAPTER_TAG}"
  GIT_SHALLOW ON
  BUILD_IN_SOURCE ON
  CONFIGURE_COMMAND ""
  BUILD_COMMAND bash "${OPENFOAM_ADAPTER_BUILD_SCRIPT}"
  INSTALL_COMMAND ""
  LOG_BUILD 1) # see <build>/cfsdeps/openfoam-adapter/src/openfoam-adapter-stamp/*build*.log

# build order: precice and OpenFOAM must be finished first
add_dependencies(openfoam-adapter openfoam precice)

add_custom_target(clean-${PACKAGE_NAME} cmake -E remove_directory "${DEPS_PREFIX}"
  COMMENT "delete cfsdeps/${PACKAGE_NAME}")

# convenience environment for running coupled cases (testsuite or manually):
#   source <build>/openfoam-env.sh
configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/openfoam-env.sh.in"
               "${CMAKE_BINARY_DIR}/openfoam-env.sh" @ONLY)

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
