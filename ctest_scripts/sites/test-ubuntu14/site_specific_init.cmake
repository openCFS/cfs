# here we only update the testsuite
set(CTEST_BUILD_NAME "Update Testsuite Trunk")

# don't delete!!
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

set(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/cfs/testsuite")
set(CTEST_BINARY_DIRECTORY "/tmp/nightly")# must be set (here the output gets written) ...

# write CTestConfig
configure_file($ENV{HOME}/cfs/trunk/ctest_scripts/shared/CTestConfig.cmake.in ${CTEST_BINARY_DIRECTORY}/CTestConfig.cmake @ONLY)

# set the name of the machine
if(CFS_BUILD_HOST STREQUAL "unknown" OR NOT CFS_BUILD_HOST)
  set(CFS_BUILD_HOST "${HOSTNAME}")
endif(CFS_BUILD_HOST STREQUAL "unknown" OR NOT CFS_BUILD_HOST)
set(CTEST_SITE "mdmt.${CFS_BUILD_HOST}")

MESSAGE("\n=============================================================================")
MESSAGE(" Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")
MESSAGE("=============================================================================\n")

# use git to update
find_program(CTEST_GIT_COMMAND NAMES git)
set(CTEST_UPDATE_TYPE "git") # use git
set(CTEST_UPDATE_COMMAND "${CTEST_GIT_COMMAND}")
set(CTEST_GIT_UPDATE_CUSTOM "${CTEST_GIT_COMMAND};svn;rebase") # use git-svn

set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
CTEST_START(Nightly)
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()
