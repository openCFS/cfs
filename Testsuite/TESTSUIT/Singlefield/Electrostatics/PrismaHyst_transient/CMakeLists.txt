#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DEPSILON:STRING=1e-5
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)

# Jiles-Atherton model is known to be unstable for multi-threading, even if pass, don't trust parallel results
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "unstable")
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "slow")
