# preCICE-coupled tests. These run cfs together with the lightweight mock
# partner (mock_fluid) instead of a heavy external solver (OpenFOAM), see
# Testsuite/CoupledTest.cmake.in. Only added when cfs is built with USE_PRECICE
# (gated in the parent TESTSUIT/CMakeLists.txt).
add_subdirectory(flapFSI2D)
add_subdirectory(SmoothPrescribedFlap)
add_subdirectory(CharacteristicBC1Dcoupled)
add_subdirectory(PerpendicularFlapSmooth)

# Tests coupling cfs with a REAL OpenFOAM solver need the cfsdeps OpenFOAM +
# openfoam-adapter (USE_OPENFOAM builds). openfoam-env.sh only exists in such
# builds - e.g. the CI test stage (packed cfs install, no OpenFOAM) skips them.
if(EXISTS "${CFS_BINARY_DIR}/openfoam-env.sh")
  add_subdirectory(PerpendicularFlapFoam)
else()
  message(STATUS "preCICE: no ${CFS_BINARY_DIR}/openfoam-env.sh (build without USE_OPENFOAM) - skipping OpenFOAM-coupled tests")
endif()
