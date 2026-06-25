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
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)
