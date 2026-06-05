#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  # We only test the absolute difference since the solution for the Lagrange multiplier is so small, numerical effects cause problems
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_STANDARD_TEST}
)
#set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken") # the test fails sporadically on serial and parallel runs in serveral build configs
