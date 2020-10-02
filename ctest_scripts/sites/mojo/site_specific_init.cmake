# update the testsuites prior to all tests

# needs to set again for the test cases
set(CTEST_SITE "mojo")
#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(ENV{LC_MESSAGES} "C")
SET(ENV{LC_ALL} "C")
SET(ENV{LANG} "C")
SET(ENV{LANGUAGE} "C")

IF(SITE_DIR MATCHES "master")
  # there group_stingl/master and cfs/master need to share the test suite
  SET(CTEST_BUILD_NAME "Update Testsuite master")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/master_cfs-test")
ELSE()
  SET(CTEST_BUILD_NAME "Update Testsuite shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/shared_cfs-test")
ENDIF()

SET(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")

# additionally we create and remove the directories in the ctest files
# note that originally this was false, whyever?!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

MESSAGE("Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")

# note that the cfs tests need to perform this again
FIND_PROGRAM(CTEST_GIT_COMMAND NAMES git)
SET(CTEST_UPDATE_TYPE "git")

CTEST_START(Nightly )
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()
MESSAGE("Finished updating testsuite")

# now we commonly update CFS source
IF(SITE_DIR MATCHES "master")
  # there group_stingl/master and cfs/master need to share the test suite
  SET(CTEST_BUILD_NAME "Update master")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/master")
ELSE()
  SET(CTEST_BUILD_NAME "Update shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/shared")
ENDIF()

message("Update ${CTEST_SOURCE_DIRECTORY}")
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()
MESSAGE("Finished updating ${CTEST_SOURCE_DIRECTORY}")
