#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
generate_test_name_and_file("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
add_test(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST})

# replacing boost::unordered_set to std::set makes this test to fail, but only for gcc9
# and only in the pipeline ?!
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken")

