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
  -P ${CFS_STANDARD_TEST}
)
# The LighthillVector test sometimes gives a wrong result, e.g. on Windows.
# As of writing is always fails iff gcc-13 MUSL stable serial.
# See more information on https://gitlab.com/openCFS/cfs/-/issues/67#note_2088857350
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken")
