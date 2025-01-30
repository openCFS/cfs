# preCICE-coupled tests. These run cfs together with the lightweight mock
# partner (mock_fluid) instead of a heavy external solver (OpenFOAM), see
# Testsuite/CoupledTest.cmake.in. Only added when cfs is built with USE_PRECICE
# (gated in the parent TESTSUIT/CMakeLists.txt).
add_subdirectory(flapFSI2D)
add_subdirectory(SmoothPrescribedFlap)
add_subdirectory(CharacteristicBC1Dcoupled)
add_subdirectory(PerpendicularFlapSmooth)

# Tests coupling cfs with a REAL OpenFOAM solver need the cfsdeps OpenFOAM +
# openfoam-adapter, i.e. a USE_OPENFOAM build. Gate on that flag (like the
# USE_PRECICE gate of this directory in ../CMakeLists.txt): whenever it is ON,
# the build provides the solvers and generates openfoam-env.sh at configure
# time, so the test is always registered - never silently skipped.
# A standalone Testsuite configure has no USE_OPENFOAM (it is not a property of
# the cfs binary, so `cfs --dependencies` cannot report it) and skips them;
# pass -DUSE_OPENFOAM=ON to opt in. The CI gcc9max test stage does exactly that:
# it unpacks the OpenFOAM runtime handed on from the build job as the
# openfoam-runtime.tar.gz artifact (OpenFOAM is GPL and never part of the cfs
# package itself) - see .gitlab-ci.yml. Other CI lanes skip them.
if(USE_OPENFOAM)
  add_subdirectory(PerpendicularFlapFoam)
else()
  message(STATUS "preCICE: USE_OPENFOAM not enabled - skipping OpenFOAM-coupled tests")
endif()
