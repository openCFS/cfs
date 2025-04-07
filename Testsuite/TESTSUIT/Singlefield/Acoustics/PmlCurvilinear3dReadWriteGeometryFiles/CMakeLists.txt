#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
add_test(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_HIST_FILE=${COMPARE_HIST_FILE}
  -DTEST_HIST_FILE:STRING="ON"
  -DTEST:STRING=${TEST_FILE_BASENAME}
  -DEPSILON:STRING=1e-10
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)