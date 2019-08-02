# shall be called on every test.
# loaded by SET_SITE_SPECIFIC_VARS() which is called by nightly_test.cmake
# and also be each master testcase. 
message("execute site_specific_vars.cmake for ${HOSTNAME}")

# for the master we test the master at TU-Wien and the copy at group_stingl
if(CTEST_SCRIPTS_DIR MATCHES "master_wien")
  set(MASTER "master_wien")
elseif(CTEST_SCRIPTS_DIR MATCHES "master_stingl")
  set(MASTER "master_stingl")
else()
  message("neither master_wien nor master_stingle, don't set MASTER")  
endif()
