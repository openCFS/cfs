# This script is run by nightly_test.cmake after the corresponding test.
# The corresponding test is determined based on the filename: 
#   see DEPLOY_TEST in test_macros.cmake
#   "${TEST_NAME}_deploy.cmake" -> "${TEST_NAME}.ctest"
# The following variables are submitted with the script execution
#   CMAKE_COMMAND, NIGHTLY_ARCHIVES_DIR, HOSTNAME, TEST_NAME
include("${CMAKE_CURRENT_LIST_DIR}/deploy_draco.cmake")
