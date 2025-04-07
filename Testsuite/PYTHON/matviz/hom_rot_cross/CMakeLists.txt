#-------------------------------------------------------------------------------
# Generate the test name from the directory name.
#-------------------------------------------------------------------------------
GENERATE_TEST_NAME_AND_FILE("${CMAKE_CURRENT_SOURCE_DIR}")

ADD_TEST(${TEST_NAME}_laminates
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_laminates
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_laminates.png --res 1500 --scale .8 --show hom_rot_cross --hom_grad none --sampling 20 ${CMAKE_CURRENT_SOURCE_DIR}/../laminates.h5ref"
  -P ${PYTHON_TEST}
)


# regression :( The cmake/python based matviz test cannot handle pdf yet :( Fabian, 30.1.2015
ADD_TEST(${TEST_NAME}_dxlaminates_pdf
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DNO_COMPARE="TRUE" # regresssion
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_dxlaminates_pdf
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_dxlaminates.pdf --res 1500 --scale .8 --show hom_rot_cross --hom_grad none --sampling 20 ${CMAKE_CURRENT_SOURCE_DIR}/../dxlaminates.h5ref"
  -P ${PYTHON_TEST}
)

ADD_TEST(${TEST_NAME}_ortho
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_ortho
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_ortho.png --res 1500 --scale .8 --show hom_rot_cross --hom_grad none --sampling 20 --parametrization ortho ${CMAKE_CURRENT_SOURCE_DIR}/../ortho.h5ref"
  -P ${PYTHON_TEST}
)

ADD_TEST(${TEST_NAME}_trans-iso
  ${CMAKE_COMMAND}
  # at the moment we cannot compare pdfs any more but we at least check for runtime errors in matviz :(
  -DEPSILON=1e-3
  -DNO_COMPARE="TRUE"
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_trans-iso
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_trans-iso.pdf --res 1500 --scale .8 --show hom_rot_cross --hom_grad none --sampling 20 --parametrization trans-iso ${CMAKE_CURRENT_SOURCE_DIR}/../trans-iso.h5ref"
  -P ${PYTHON_TEST}
)

ADD_TEST(${TEST_NAME}_dxtrans-iso
  ${CMAKE_COMMAND}
  -DEPSILON=1e-3
  -DCURRENT_TEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -DTEST:STRING=${TEST_FILE_BASENAME}_dxtrans-iso
  -DMATVIZ_ARGS:STRING="--save ${TEST_FILE_BASENAME}_dxtrans-iso.png --res 1500 --scale .8 --show hom_rot_cross --hom_grad none --sampling 20 --color gist_heat --parametrization trans-iso ${CMAKE_CURRENT_SOURCE_DIR}/../dxtrans-iso.h5ref"
  -P ${PYTHON_TEST}
)

