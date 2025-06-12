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
  # We only test the absolute difference since the smoothVelocity computation from the previous sequence step should give us a 0 solution. Due to numerics, the solution is in the range of e-14 leading to large relative errors.
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_STANDARD_TEST}
)
