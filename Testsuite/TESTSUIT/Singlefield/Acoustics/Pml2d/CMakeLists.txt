#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Test: original base test (defaults to autoConsolidated mode)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CFS_STANDARD_TEST}
)

#-------------------------------------------------------------------------------
# Test: autoConsolidated mode (frequency-dependent parts consolidated in SYSTEM_UPDATE)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}_autoConsolidated
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DCFS_ARGS=-pPml2d_autoConsolidated.xml
  -P ${CFS_STANDARD_TEST}
)

#-------------------------------------------------------------------------------
# Test: autoSplit mode (frequency-dependent parts in separate UPDATE matrices)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}_autoSplit
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DCFS_ARGS=-pPml2d_autoSplit.xml
  -P ${CFS_STANDARD_TEST}
)

#-------------------------------------------------------------------------------
# Test: noReassembly mode (testing only, expects warning)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}_noReassembly
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DCFS_ARGS=-pPml2d_noReassembly.xml
  -P ${CFS_STANDARD_TEST}
)

#-------------------------------------------------------------------------------
# Test: fullReassembly mode (legacy, always reassemble into SYSTEM)
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}_fullReassembly
  ${CMAKE_COMMAND} 
  -DEPSILON:STRING=${CFS_DEFAULT_EPSILON}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DCFS_ARGS=-pPml2d_fullReassembly.xml
  -P ${CFS_STANDARD_TEST}
)
