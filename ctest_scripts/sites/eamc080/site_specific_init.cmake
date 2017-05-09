# update the shared_opt testsuite prior to all tests

# needs to set again for the test cases
set(CTEST_SITE "eamc080")
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
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/cfstest-trunk")
ELSE()
  SET(CTEST_BUILD_NAME "Update Testsuite shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/cfstest-shared_opt")
ENDIF()

# this is copy pasted to the ctest cases
set(CTEST_BUILD_NAME "Update Testsuite shared_opt")
SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/cfstest-shared_opt")
SET(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")

# don't delete!!
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

MESSAGE("\n---------------------------------------------------")
MESSAGE("Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")
MESSAGE("---------------------------------------------------\n")

FIND_PROGRAM(CTEST_SVN_COMMAND NAMES svn)
SET(CTEST_UPDATE_TYPE "svn")

CTEST_START(Nightly)
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()

