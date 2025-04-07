#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for info.xml and not h5 as the modes are not uniquely defined
# The problem solves for every compiler and solver differently, a very ill-posed problem
# relax the precission and ignore the actual number of iterations (--last)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML} 
  -DEPSILON=0.05
  -DLAST=--last
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST_INFO_XML:STRING="ON"
  -P ${CFS_STANDARD_TEST}
)
