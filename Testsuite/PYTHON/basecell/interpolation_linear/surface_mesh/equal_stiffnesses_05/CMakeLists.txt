#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}
  ${CMAKE_COMMAND}
  -DCOMPARE_INFO_XML=${COMPARE_INFO_XML} 
  -DEPSILON=1e-4
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST_INFO_XML:STRING="OFF"
  -DTEST:STRING=${TEST_FILE_BASENAME}
  -DBASECELL_ARGS:STRING=" --res 40 --target surface_mesh --x1 0.5 --y1 0.5 --z1 0.5 --interpolation linear --save_vtp  --beta 7 --eta 0.6 --to_info_xml --save ${TEST_FILE_BASENAME} --tets --smooth_iter 25 --smooth_lambda 0.6"
  -DCFS_HOMOGENIZE_ARGS:STRING=" -d -m ${TEST_FILE_BASENAME}.mesh -p ${TEST_FILE_BASENAME}_homogenize.xml ${TEST_FILE_BASENAME}_homogenize" 
  -P ${PYTHON_TEST}
)
