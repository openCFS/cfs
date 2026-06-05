#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=1e-2
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)

# Was unstable but with MKL_CBWR=COMPATIBLE it always fail. 
# It works with MKL_CBWR=AUTO, hence the code is good, 
# but the numerics is extremely sensitive as a changed MKL path makes the difference
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken")
