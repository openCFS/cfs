# shall be called on every test.
# loaded by SET_SITE_SPECIFIC_VARS() which is called by nightly_test.cmake
message("execute site_specific_vars.cmake")
message("SITE_DIR=${SITE_DIR}")

# for the master we test the master at TU-Wien and the copy at group_stingl
if(SITE_DIR MATCHES "master_wien")
  set(MASTER "master_wien")
if(SITE_DIR MATCHES "master_wien")
  set(MASTER "master_stingl")
endif()

message("find and set git")
find_program(CTEST_GIT_COMMAND NAMES git) 
set(CTEST_UPDATE_TYPE "git")
