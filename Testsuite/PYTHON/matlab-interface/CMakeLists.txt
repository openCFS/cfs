#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define tests for this directory
#-------------------------------------------------------------------------------
# The test command is simply python followed by the module. The module must invoke doctest in it's __main__ method. All examples given in docstrings will be tested automatically.
ADD_TEST("${TEST_NAME}"
  ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/compute_square.py" -v
  WORKING_DIRECTORY "${TESTSUITE_DIR}"
)
