# update the shared_opt testsuite prior to all tests

# needs to set again for the test cases
set(CTEST_SITE "momo")
#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(ENV{LC_MESSAGES} "C")
SET(ENV{LC_ALL} "C")
SET(ENV{LANG} "C")
SET(ENV{LANGUAGE} "C")

IF(${SITE_DIR} MATCHES "master_stingl")
  SET(CTEST_BUILD_NAME "Update Testsuite stingl master")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/testsuites/test-master_stingl")
ELSE()
  SET(CTEST_BUILD_NAME "Update Testsuite shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/testsuites/test-shared_opt")
ENDIF()

SET(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}")

# don't know if this even works!!
# create new empty directory in test specific cmake files, not here!
# SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

MESSAGE("\n---------------------------------------------------")
MESSAGE("Update testsuite ${CTEST_SOURCE_DIRECTORY} ...")
MESSAGE("---------------------------------------------------\n")

FIND_PROGRAM(CTEST_GIT_COMMAND NAMES git)
SET(CTEST_UPDATE_TYPE "git")

CTEST_START(Nightly)
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()

# now we commonly update CFS source
IF(SITE_DIR MATCHES "master")
  # there group_stingl/master and cfs/master need to share the test suite
  SET(CTEST_BUILD_NAME "Update master_stingl")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/master_stingl")
ELSE()
  SET(CTEST_BUILD_NAME "Update shared_opt")
  SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/code/shared_opt")
ENDIF()

message("Update ${CTEST_SOURCE_DIRECTORY}")
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()

MESSAGE("Finished updating ${CTEST_SOURCE_DIRECTORY}")

