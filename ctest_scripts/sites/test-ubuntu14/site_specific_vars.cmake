# This script is loaded by `nightly_test.cmake`, it sets source and build
# directories, the hostname and the update method.
# To make the tests run individually also it must be included by
# the single test files too.

# these variables must be set:
set(NIGHTLY_TESTSUITE_DIR "$ENV{HOME}/cfs/testsuite") # source of the testsuite
set(NIGHTLY_CACHE_DIR "$ENV{HOME}/cfs/build/cache") # cache dir for builds
set(NIGHTLY_SOURCE_DIR "$ENV{HOME}/cfs/master") # source of cfs
set(NIGHTLY_BUILD_DIR "$ENV{HOME}/cfs/build") # where the builds should be made (a sub-path is appended for each test)

# set the name of the machine (Site)
if(NOT CFS_BUILD_HOST)
  site_name(CFS_BUILD_HOST) # sets CFS_BUILD_HOST to the name of the computer
  string(REPLACE "." ";" CFS_BUILD_HOST_LIST "${CFS_BUILD_HOST}") # puts it into list
  list(GET CFS_BUILD_HOST_LIST 0 CFS_BUILD_HOST) # gets first element
endif(NOT CFS_BUILD_HOST)
set(CTEST_SITE "mdmt.${CFS_BUILD_HOST}") # use mdmt in front to allow for sorting in cdash

# use git to update
find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_TYPE "git") # use git
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
