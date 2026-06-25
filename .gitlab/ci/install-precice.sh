#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# install-precice.sh - provision preCICE on a build/test runner.
#
# Builds preCICE from source with a minimal feature set (no PETSc, no MPI, no
# Python actions) against the system boost / eigen / libxml2, and installs it to
# a prefix that openCFS' `find_package(precice CONFIG)` can locate.
#
# Intended uses:
#   * baking preCICE into the openCFS build-environment Docker images
#     (https://gitlab.com/openCFS/build-environment), or
#   * a one-off provisioning step on a developer machine / CI runner.
#
# Usage:
#   PRECICE_VERSION=v3.1.2 PRECICE_PREFIX=/opt/precice ./install-precice.sh
#   ./install-precice.sh --no-deps        # skip system package installation
#
# After installation, configure openCFS with:
#   -DUSE_PRECICE=ON -DCMAKE_PREFIX_PATH=$PRECICE_PREFIX
# (or -Dprecice_DIR=$PRECICE_PREFIX/lib/cmake/precice).
#
# IMPORTANT (reproducibility): the preCICE version used here must match the
# version used to generate the committed test reference (.h5ref). The
# Testsuite/TESTSUIT/preCICE reference was produced with preCICE 3.x; pin
# PRECICE_VERSION accordingly and regenerate references if you change it.
# ---------------------------------------------------------------------------
set -euo pipefail

PRECICE_VERSION="${PRECICE_VERSION:-v3.1.2}"
PRECICE_PREFIX="${PRECICE_PREFIX:-/opt/precice}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 2)}"
INSTALL_DEPS=1
[ "${1:-}" = "--no-deps" ] && INSTALL_DEPS=0

echo "=== preCICE ${PRECICE_VERSION} -> ${PRECICE_PREFIX} (jobs=${JOBS}) ==="

# --- 1) system build dependencies (boost, eigen, libxml2, build tools) -------
if [ "${INSTALL_DEPS}" = "1" ]; then
  if command -v apt-get >/dev/null 2>&1; then           # Debian / Ubuntu
    export DEBIAN_FRONTEND=noninteractive
    apt-get update
    apt-get install -y --no-install-recommends \
      git cmake g++ make ca-certificates \
      libboost-all-dev libeigen3-dev libxml2-dev
  elif command -v dnf >/dev/null 2>&1; then             # Fedora / Rocky / RHEL
    dnf install -y \
      git cmake gcc-c++ make \
      boost-devel eigen3-devel libxml2-devel
  else
    echo "WARN: unknown package manager - assuming boost/eigen/libxml2 are present" >&2
  fi
fi

# --- 2) fetch + build preCICE (minimal features) -----------------------------
workdir="$(mktemp -d)"
trap 'rm -rf "${workdir}"' EXIT
git clone --depth 1 --branch "${PRECICE_VERSION}" https://github.com/precice/precice.git "${workdir}/precice"

cmake -S "${workdir}/precice" -B "${workdir}/precice/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PRECICE_PREFIX}" \
  -DPRECICE_FEATURE_PETSC_MAPPING=OFF \
  -DPRECICE_FEATURE_MPI_COMMUNICATION=OFF \
  -DPRECICE_FEATURE_PYTHON_ACTIONS=OFF \
  -DPRECICE_BINDINGS_PYTHON=OFF \
  -DBUILD_TESTING=OFF

cmake --build "${workdir}/precice/build" -j "${JOBS}" --target install

echo "=== preCICE installed ==="
echo "  precice_DIR = ${PRECICE_PREFIX}/lib/cmake/precice"
ls "${PRECICE_PREFIX}/lib/cmake/precice" 2>/dev/null || \
  ls "${PRECICE_PREFIX}/lib64/cmake/precice" 2>/dev/null || true
