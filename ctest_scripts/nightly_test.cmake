# This script is the main entry point for automated nightly builds of CFS++.
# It gets either called from /etc/crontab or through the provisioner scripts
# for Vagrant VBoxes.
#
# To call this script by hand, it is advisable to simulate a crontab
# environment e.g. in the following way:
# env -i HOME=/ PATH=/bin:/usr/bin:/sbin:/usr/sbin \
#        /opt/pckg/cmake-2.8.9/bin/ctest -S $HOME/Documents/dev/nightly_test.cmake
#
# vagrant ssh -c 'env -i HOME=/Users/simon PATH=/bin:/usr/bin:/sbin:/usr/sbin env && uname -a && /opt/pckg/cmake-2.8.10.2-Linux-i386/bin/cmake --version'

# We need at least CMake 2.8.9 for this to work
CMAKE_MINIMUM_REQUIRED(
  VERSION "2.8.9"
  FATAL_ERROR
)

# List command no longer ignores empty elements.
CMAKE_POLICY(SET CMP0007 NEW)

# Get base path of current script in order to include additional macros.
GET_FILENAME_COMPONENT(CTEST_SCRIPTS_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

# Include informations about development server.
INCLUDE("${CTEST_SCRIPTS_DIR}/../cmake_modules/DevelopmentServer.cmake")

# Include further macros required for testing.
INCLUDE("${CTEST_SCRIPTS_DIR}/shared/test_macros.cmake")

# Determine date and time.
EXECUTE_PROCESS(
  COMMAND date "+%F %H:%M:%S"
  OUTPUT_VARIABLE DATE_OUT
)
STRING(REPLACE "\n" "" DATE_OUT ${DATE_OUT})
STRING(STRIP ${DATE_OUT} DATE_OUT)
MESSAGE(
"
=============================================================================
 Starting nightly tests on ${DATE_OUT}...
=============================================================================
"
)

# Set global variables, e.g. path to ctest exe, site base dir, host name etc.
SET_GLOBAL_VARS()

# Set site specific variables, e.g. test user home dir, svn user password, etc.
SET_SITE_SPECIFIC_VARS()

MESSAGE(
"
=============================================================================
 Global variables on ${HOSTNAME}...
=============================================================================
"
)
MESSAGE("CMAKE_COMMAND: ${CMAKE_COMMAND}")
MESSAGE("CMAKE_EXECUTABLE_SUFFIX: ${CMAKE_EXECUTABLE_SUFFIX}")
MESSAGE("CTEST_COMMAND: ${CTEST_COMMAND}")
MESSAGE("CMAKE_CURRENT_LIST_FILE: ${CMAKE_CURRENT_LIST_FILE}")
MESSAGE("SITE_BASE_DIR: ${SITE_BASE_DIR}")
MESSAGE("NIGHTLY_ARCHIVES_DIR: ${NIGHTLY_ARCHIVES_DIR}")
MESSAGE("HOSTNAME: ${HOSTNAME}")
MESSAGE("SITE_DIR: ${SITE_DIR}")
MESSAGE("TESTUSER: ${TESTUSER}")
MESSAGE("HOME: ${HOME}")
MESSAGE("DAYOFWEEK: ${DAYOFWEEK}")

# Perform site specific initialization tasks, e.g. update working copies,
# start VBoxes, etc.
SITE_SPECIFIC_INIT()

# Get list of tests to be performed on this site.
GET_TEST_NAMES(TEST_NAMES)

# Iterate over list of tests.
FOREACH(TEST_NAME IN ITEMS ${TEST_NAMES})

  GET_CTEST_BINARY_DIRECTORY(${TEST_NAME})

  # Actually run the test
  PERFORM_TEST(${TEST_NAME})
  
  # E.g. pack generated binaries to a globally accessible zip archive.
  DEPLOY_TEST(${TEST_NAME})

  # Perform some test specific cleanup chores.
  CLEANUP_TEST(${TEST_NAME})

ENDFOREACH()

# Perform site specific finalization tasks, such as unpacking nightly
# binaries to a global network share.
SITE_SPECIFIC_FINISH()

