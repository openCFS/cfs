#-------------------------------------------------------------------------------
# Standalone regression test for the acoustic characteristicCouplingBC (no preCICE,
# no OpenFOAM): a 1-cycle 5 kHz burst in a short duct is absorbed at S_out by the
# BC driven with w_in = 0 (must reduce to a first-order absorbing BC). Geometry/time
# are downsized from the original 1.5 m / 480-step validation case so the run is
# under 1 s and the reference under 1 MB. See CharacteristicBC1D.xml / make_mesh.py.
# The preCICE-coupled (two-participant) variant lives under TESTSUIT/preCICE.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DEPSILON:STRING=1e-4
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  # Absolute difference: acouCharacteristic at S_out is machine noise (< 1e-14)
  # until the burst arrives (~step 25), so a relative comparison of those steps
  # is solver/platform noise. 1e-4 absolute still locks the physical signal
  # (acouPressure L2 up to ~600, acouCharacteristic peak ~0.5) tightly.
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_STANDARD_TEST}
)
