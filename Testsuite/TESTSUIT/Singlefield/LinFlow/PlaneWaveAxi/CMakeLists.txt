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

# This test fails sporadically in parallel runs on Linux, on Windows (MSVC2019,ifort2021) it only passes in 1/8 times -> broken
if(WIN32)
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "broken")
else()
set_property(TEST ${TEST_NAME} APPEND PROPERTY LABELS "unstable")
endif()
