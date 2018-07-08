# shall be called on every test.
# loaded by SET_SITE_SPECIFIC_VARS() which is called by nightly_test.cmake
# and also be each master testcase. 
message("execute site_specific_vars.cmake")

# for the master we test the master at TU-Wien and the copy at group_stingl
if(CTEST_SCRIPTS_DIR MATCHES "master_wien")
  set(MASTER "master_wien")
  set("set MASTER=${MASTER}")
elseif(CTEST_SCRIPTS_DIR MATCHES "master_stingl")
  set(MASTER "master_stingl")
  set("set MASTER=${MASTER}")
else()
  message("MASTER not set, CTEST_SCRIPTS_DIR = ${CTEST_SCRIPTS_DIR}")  
endif()

find_program(CTEST_GIT_COMMAND NAMES git) 
set(CTEST_UPDATE_TYPE "git")
