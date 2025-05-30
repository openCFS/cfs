#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML}
  -DTEST_INFO_XML:STRING="ON"
  -DADD_TEST_FILES:STRING=../mtimesx-master
  -P ${CFS_STANDARD_TEST}
)
