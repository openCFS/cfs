# preCICE-coupled tests. These run cfs together with the lightweight mock
# partner (mock_fluid) instead of a heavy external solver (OpenFOAM), see
# Testsuite/CoupledTest.cmake.in. Only added when cfs is built with USE_PRECICE
# (gated in the parent TESTSUIT/CMakeLists.txt).
add_subdirectory(flapFSI2D)
add_subdirectory(SmoothPrescribedFlap)
add_subdirectory(CharacteristicBC1Dcoupled)
