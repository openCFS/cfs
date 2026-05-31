#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

set(TEST_SCRIPT "${CFS_BINARY_DIR}/share/python/postprocess_2d5_to_3d.py")

ADD_TEST(NAME ${TEST_NAME}
  COMMAND ${CMAKE_COMMAND}
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}
  -DPYTHON_ARGS:STRING="${TEST_SCRIPT} -o ${TESTSUITE_BIN_DIR}/${CURRENT_TEST_SUBDIR}/results_hdf5/${TEST_FILE_BASENAME}.h5 --z '0.1,0.2,0.3,0.4,0.5' --frequency 100.0 SphericalWave100HzSpectrum1.cfs SphericalWave100HzSpectrum2.cfs SphericalWave100HzSpectrum3.cfs"
  -DTEST_H5DIFF:BOOL=ON
  -DEPSILON=1e-4
  -P ${PYTHON_TEST}
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
