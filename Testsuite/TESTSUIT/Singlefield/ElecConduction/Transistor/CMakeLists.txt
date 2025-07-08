#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=1e-6
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  # We only test the absolute difference since at the beginning the current density is quite small leading to large relative errors. The standard tolerance of 1e-6 will still adequately check the main features of the test when transistor gets activated, since the current density reaches values up to 3e-3 A/m^2.
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_STANDARD_TEST}
)
