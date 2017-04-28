# update the shared_opt testsuite prior to all tests

set(CTEST_SITE "eamc061")
#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(ENV{LC_MESSAGES} "C")
SET(ENV{LC_ALL} "C")
SET(ENV{LANG} "C")
SET(ENV{LANGUAGE} "C")

IF(SITE_DIR MATCHES "trunk")
  SET(CTEST_BUILD_NAME "Update Testsuite trunk")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/trunk_cfs-test")
ELSE()
  SET(CTEST_BUILD_NAME "Update Testsuite shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/shared_cfs-test")
ENDIF()

SET(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")

# additionally we create and remove the directories in the ctest files
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

MESSAGE("\n---------------------------------------------------")
MESSAGE("Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")
MESSAGE("---------------------------------------------------\n")

FIND_PROGRAM(CTEST_SVN_COMMAND NAMES svn)
SET(CTEST_UPDATE_TYPE "svn")

CTEST_START(Nightly )
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()

