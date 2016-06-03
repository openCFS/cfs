# here we only update the testsuite
set(CTEST_BUILD_NAME "Update Testsuite trunk")

# don't delete!!
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

set(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/cfs/testsuite")
set(CTEST_BINARY_DIRECTORY "/tmp/nightly")# must be set (here the output gets written) ...

# set the name of the machine
if(NOT CFS_BUILD_HOST)
  site_name(CFS_BUILD_HOST) # sets CFS_BUILD_HOST to the name of the computer
  string(REPLACE "." ";" CFS_BUILD_HOST_LIST "${CFS_BUILD_HOST}") # puts it into list
  list(GET CFS_BUILD_HOST_LIST 0 CFS_BUILD_HOST) # gets first element
endif(NOT CFS_BUILD_HOST)
set(CTEST_SITE "${CFS_BUILD_HOST}")

MESSAGE("\n=============================================================================")
MESSAGE(" Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")
MESSAGE("=============================================================================\n")

# use git to update
find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_TYPE "git") # use git
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_UPDATE_CUSTOM "${CTEST_GIT_COMMAND};svn;rebase") # use git-svn

CTEST_START(Nightly )
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()
