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
  -DCFSTOOL_MODE=absL2diff
  -P ${CFS_STANDARD_TEST}
)

# unstable for multi-threading
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "unstable")
