#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

#-------------------------------------------------------------------------------
# Define a test for this directory
#-------------------------------------------------------------------------------
ADD_TEST(${TEST_NAME}_mandel
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_mandel
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_mandel.png --res 1500 --sampling 60 --notation mandel ${CMAKE_CURRENT_SOURCE_DIR}/../stiff1_stiff2.h5ref"
  -P ${PYTHON_TEST}
)

ADD_TEST(${TEST_NAME}_voigt
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_voigt
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_voigt.png --res 1500 ${CMAKE_CURRENT_SOURCE_DIR}/../stiff1_stiff2.h5ref"
  -P ${PYTHON_TEST}
)

ADD_TEST(${TEST_NAME}_scaled
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_scaled
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_scaled.png --res 1500 --sampling 60 --scale 12 ${CMAKE_CURRENT_SOURCE_DIR}/../stiff1_stiff2.h5ref"
  -P ${PYTHON_TEST}
)


