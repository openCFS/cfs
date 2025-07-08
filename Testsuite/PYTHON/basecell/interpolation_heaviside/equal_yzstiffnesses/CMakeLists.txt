#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML} 
  -DEPSILON=1e-4
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST_INFO_XML:STRING="ON"
  -DTEST:STRING=${TEST_FILE_BASENAME}
  -DBASECELL_ARGS:STRING=" --res 40 --target volume_mesh --to_info_xml --x1 0.2 --x2 0.2 --y1 0.01 --y2 0.01 --z1 0.01 --z2 0.01 --save ${TEST_FILE_BASENAME} --interpolation heaviside"
  -DCFS_HOMOGENIZE_ARGS:STRING=" -d -m ${TEST_FILE_BASENAME}.mesh -p ${TEST_FILE_BASENAME}_homogenize.xml ${TEST_FILE_BASENAME}_homogenize" 
  -P ${PYTHON_TEST}
)
