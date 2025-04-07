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
  -DMATVIZ_ARGS:STRING="  ${CMAKE_CURRENT_SOURCE_DIR}/HomRect_C1_3D.h5ref --show hom_ortho_3d --save ${TEST_FILE_BASENAME}.vtp  --bc_res 20 --bc_interpolation heaviside --hom_samples 2,2,2 --mesh v --h5_nondes None --bc_smooth 5"
  -DNO_COMPARE="TRUE"
  -DINTERPRETATION_ARGS:STRING="ON"
  -DTHREADS:STRING=4
  -P ${PYTHON_TEST}
)
