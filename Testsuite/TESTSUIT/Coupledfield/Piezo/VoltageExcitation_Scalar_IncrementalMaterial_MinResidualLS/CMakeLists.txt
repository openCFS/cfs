#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=5e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)

# fails sometimes for gcc7 max parallel and almost any time with gcc12/13 - even serial
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken")
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "slow")
