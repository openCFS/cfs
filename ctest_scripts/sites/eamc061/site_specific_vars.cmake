# shall be called on every test.
# loaded by SET_SITE_SPECIFIC_VARS() which is called by nightly_test.cmake
# and also be each master testcase. 
message("execute site_specific_vars.cmake")

message("CTEST_SCRIPTS_DIR = ${CTEST_SCRIPTS_DIR}")

# shall be already set in the test to be able to find this file!
if(NOT CTEST_SCRIPTS_DIR)
  GET_FILENAME_COMPONENT(CTEST_SCRIPTS_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
  message("set CTEST_SCRIPTS_DIR = ${CTEST_SCRIPTS_DIR}")
endif()  

# for the master we test the master at TU-Wien and the copy at group_stingl
if(CTEST_SCRIPTS_DIR MATCHES "master_wien")
  set(MASTER "master_wien")
elseif(CTEST_SCRIPTS_DIR MATCHES "master_wien")
  set(MASTER "master_stingl")
else()
  message("MASTER not set, CTEST_SCRIPTS_DIR = ${CTEST_SCRIPTS_DIR}")  
endif()

find_program(CTEST_GIT_COMMAND NAMES git) 
set(CTEST_UPDATE_TYPE "git")
