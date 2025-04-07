#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

# not verified VTK vtp output
ADD_TEST(${TEST_NAME}_vtp
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DNO_COMPARE="TRUE" # a comparison is not yet implemented
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_vtp
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_vtp.vtp --show hom_rect --hom_samples 2 --hom_grad linear ${CMAKE_CURRENT_SOURCE_DIR}/lattice_2x2x2.h5"
  -P ${PYTHON_TEST}
)
