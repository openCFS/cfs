#------------------------------------------------------------------------------
# Generate the test name from the directory name.
#------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

set(TEST_SCRIPT "${CFS_BINARY_DIR}/share/python/extrude_mesh.py")

ADD_TEST(NAME ${TEST_NAME}
  COMMAND ${CMAKE_COMMAND}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}
  -DPYTHON_ARGS:STRING="${TEST_SCRIPT} Mesh2DLinElem.cfs ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}.h5 -1.0 -0.5 0.5 1.0"
  -DTEST_H5DIFF:BOOL=ON
  -DH5DIFF_OBJ:STRING=/Mesh
  -P ${PYTHON_TEST}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
