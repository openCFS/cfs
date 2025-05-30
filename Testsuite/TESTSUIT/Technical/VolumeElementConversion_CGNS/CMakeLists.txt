#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  set(ALLOW_CFS_FAIL ON)
  message(STATUS "setting PROCEED_AFTER_SIMULATION_CRASH=${ALLOW_CFS_FAIL} for test '${TEST_NAME}' on '${CMAKE_SYSTEM_NAME}' because it fails on Linux-GCC13-MSUL build (segfault in destructor)")
endif()
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCFS_ARGS:STRING=-gG
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST_INFO_XML:STRING="ON"
  -DSKIP_H5REF:STRING="ON"
  -DPROCEED_AFTER_SIMULATION_CRASH:BOOL=${ALLOW_CFS_FAIL} # test will continue with result check even if simulation run has non-zero exit code
  -P ${CFS_STANDARD_TEST}
)