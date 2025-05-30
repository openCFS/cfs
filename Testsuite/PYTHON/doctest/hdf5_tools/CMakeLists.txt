#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define tests for this directory
#-------------------------------------------------------------------------------
set(TEST_SCRIPT "${CFS_BINARY_DIR}/share/python/hdf5_tools.py")
ADD_TEST(NAME "${TEST_NAME}"
  COMMAND ${Python_EXECUTABLE} "${TEST_SCRIPT}" -v
  WORKING_DIRECTORY "${TESTSUITE_DIR}"
)
