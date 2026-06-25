#-------------------------------------------------------------------------------
# Standalone test of the SmoothPDE "prescribed displacement" mode (no preCICE):
# an analytic whole-domain displacement is injected into SMOOTH_DISPLACEMENT every
# step WITHOUT solving. A non-zero result matching the prescribed field (rather than
# u = 0 from a real solve with no BC/load) proves the field is injected, not solved.
# See SmoothPrescribed1D.xml. The preCICE-fed variant lives under TESTSUIT/preCICE.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)
