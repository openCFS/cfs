#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

# h5 would only test the test strains but we want to compare the tensor
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML} 
  -DEPSILON=1e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST_INFO_XML:STRING="ON"
  -P ${CFS_STANDARD_TEST}
)

# this test fails sporadically in parallel runs, hence, we set it to unstable
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "unstable")
