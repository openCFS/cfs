# OpenFOAM (openfoam.com / ESI release line) - the CFD partner solver for real
# coupled openCFS<->OpenFOAM simulations via preCICE (see cfsdeps/openfoam-adapter).
#
# Pure build/test dependency: openCFS does NOT link OpenFOAM (it is GPL and is only
# run as an external program by the coupled testsuite). Hence, unlike other cfsdeps,
# there is no PACKAGE_LIBRARY and nothing of it is shipped in the cfs package -
# which also keeps the GPL code out of the closed distribution.
#
# OpenFOAM is a "large dependency" (like ginkgo/vtk on Windows): the CI builds it
# in a separate largedeps job and passes the precompiled package to the regular
# cfsdeps/build jobs, which just unpack it (see .gitlab-ci.yml). The package is a
# .tar.gz (NOT the usual .zip: the lnInclude trees are symlink farms and zip would
# store file copies, breaking relative links and tripling the size). The tree is
# relocatable because etc/bashrc self-locates and the binaries resolve their libs
# via the sourced environment, not via rpaths. The openfoam-adapter is deliberately
# NOT part of the package: it rebuilds in ~2 min against the unpacked tree, which
# also re-links it with the rpath of the consuming build dir.
clear_depencency_variables()

set(PACKAGE_NAME "openfoam")

# ---- OpenFOAM version - keep easily editable -------------------------------
# source tarballs: https://dl.openfoam.com/source/<version>/
set(OPENFOAM_VER "v2512")
set(OPENFOAM_MD5 "8a1abe3864851902bb77809efa94ba3e") # OpenFOAM-v2512.tgz
# -----------------------------------------------------------------------------

set(PACKAGE_VER "${OPENFOAM_VER}")
set(PACKAGE_FILE "OpenFOAM-${OPENFOAM_VER}.tgz")
set(PACKAGE_MD5 "${OPENFOAM_MD5}")
set(PACKAGE_MIRRORS "https://dl.openfoam.com/source/${OPENFOAM_VER}/${PACKAGE_FILE}")
# add default mirrors to PACKAGE_MIRRORS or replace all with LOCAL_PACKAGE_FILE if we already have it
add_standard_mirrors_or_set_local()

set(DEPS_PREFIX "${CMAKE_BINARY_DIR}/cfsdeps/${PACKAGE_NAME}")
set(DEPS_SOURCE "${DEPS_PREFIX}/src/${PACKAGE_NAME}") # becomes WM_PROJECT_DIR

# consumed by cfsdeps/openfoam-adapter and the coupled testsuite
set(OPENFOAM_DIR "${DEPS_SOURCE}" CACHE INTERNAL "OpenFOAM project dir (WM_PROJECT_DIR)" FORCE)

# OpenFOAM builds with the system gcc and the system OpenMPI (WM_MPLIB=SYSTEMOPENMPI,
# the same MPI cfs/preCICE use). mpicc is not necessarily in the default PATH
# (e.g. Fedora installs it to /usr/lib64/openmpi/bin).
find_program(OPENFOAM_MPICC mpicc HINTS /usr/lib64/openmpi/bin /usr/lib/openmpi/bin)
mark_as_advanced(OPENFOAM_MPICC)
if(OPENFOAM_MPICC)
  get_filename_component(OPENFOAM_MPI_BINDIR "${OPENFOAM_MPICC}" DIRECTORY)
else()
  set(OPENFOAM_MPI_BINDIR "")
  message(STATUS "openfoam: no mpicc found - the OpenFOAM build relies on the caller's environment")
endif()

# standard precompiled cache name, but as .tar.gz (symlinks, see header comment).
# use_c_and_fortran only feeds the naming; OpenFOAM builds with the system gcc.
use_c_and_fortran(ON OFF)
set_precompiled_pckg_file()
string(REGEX REPLACE "\\.zip$" ".tar.gz" PRECOMPILED_PCKG_FILE "${PRECOMPILED_PCKG_FILE}")

# the build is multi-step (source etc/bashrc, Allwmake, sanity check): generate a
# script instead of fighting cmake COMMAND quoting; it can also be re-run manually.
set(OPENFOAM_BUILD_SCRIPT "${DEPS_PREFIX}/openfoam-build.sh")
configure_file("${CMAKE_SOURCE_DIR}/cfsdeps/${PACKAGE_NAME}/openfoam-build.sh.in"
               "${OPENFOAM_BUILD_SCRIPT}" @ONLY)

if(${CFS_DEPS_PRECOMPILED} AND EXISTS "${PRECOMPILED_PCKG_FILE}")
  # unpack the prebuilt OpenFOAM tree instead of building (CI: passed as artifact
  # from the largedeps job / found in the cfsdeps cache)
  ExternalProject_Add(openfoam
    PREFIX "${DEPS_PREFIX}"
    DOWNLOAD_COMMAND ""
    UPDATE_COMMAND ""
    PATCH_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND "")
  ExternalProject_Add_Step(openfoam unpack_precompiled
    COMMAND ${CMAKE_COMMAND} -E echo "unpacking ${PRECOMPILED_PCKG_FILE}"
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DEPS_SOURCE}
    COMMAND ${CMAKE_COMMAND} -E chdir ${DEPS_SOURCE} ${CMAKE_COMMAND} -E tar xzf ${PRECOMPILED_PCKG_FILE}
    DEPENDERS build)
else()
  ExternalProject_Add(openfoam
    PREFIX "${DEPS_PREFIX}"
    SOURCE_DIR "${DEPS_SOURCE}"
    URL "${PACKAGE_MIRRORS}"
    URL_MD5 "${PACKAGE_MD5}"
    # DOWNLOAD_DIR is ignored, if URL contains not the file mirror
    DOWNLOAD_DIR "${CFS_DEPS_CACHE_DIR}/sources/${PACKAGE_NAME}"
    DOWNLOAD_NAME "${PACKAGE_FILE}"
    DOWNLOAD_NO_PROGRESS ON
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND ""
    BUILD_COMMAND bash "${OPENFOAM_BUILD_SCRIPT}"
    INSTALL_COMMAND ""
    # wmake output is very long - keep it in the stamp log
    LOG_BUILD 1) # see <build>/cfsdeps/openfoam/src/openfoam-stamp/openfoam-build-*.log

  if(${CFS_DEPS_PRECOMPILED})
    # pack the runtime + adapter-devel part of the tree: platforms (bin/lib/tools),
    # src (headers + lnInclude symlinks, the adapter compiles against them), wmake
    # (toolchain to build the adapter), etc, bin, META-INFO. Deliberately NOT
    # build/ (~2 GB intermediates), tutorials/, applications/, modules/, doc/.
    ExternalProject_Add_Step(openfoam store_precompiled
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CFS_DEPS_CACHE_DIR}/precompiled
      COMMAND ${CMAKE_COMMAND} -E chdir ${DEPS_SOURCE} ${CMAKE_COMMAND} -E tar czf ${PRECOMPILED_PCKG_FILE} Allwmake bin etc META-INFO platforms src wmake
      DEPENDEES build)
  endif()
endif()

# no set_standard_variables() here (it requires PACKAGE_LIBRARY), so provide the
# clean-<package> convenience target manually (like the standard one it also
# drops the precompiled package to force a real rebuild)
add_custom_target(clean-${PACKAGE_NAME} cmake -E remove_directory "${DEPS_PREFIX}"
  COMMAND cmake -E remove -f "${PRECOMPILED_PCKG_FILE}"
  COMMENT "delete cfsdeps/${PACKAGE_NAME} (sources + build) and its precompiled package")

# add project to global list of CFSDEPS
set(CFSDEPS ${CFSDEPS} ${PACKAGE_NAME})
