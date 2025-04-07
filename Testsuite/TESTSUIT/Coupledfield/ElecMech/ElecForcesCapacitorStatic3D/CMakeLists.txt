#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  # relax tolerance for netlib  
  -DEPSILON:STRING=1e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)
