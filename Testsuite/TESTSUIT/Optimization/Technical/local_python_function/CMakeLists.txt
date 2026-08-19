#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for the last iteration in info.xml and not h5. Only the info.xml has
# max_abs_slope and infeas_count_slope, which is what the python local function does.
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML}
  -DEPSILON:STRING=1e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DLAST=--last
  -DTEST_INFO_XML:STRING="ON"
  -P ${CFS_STANDARD_TEST}
)
