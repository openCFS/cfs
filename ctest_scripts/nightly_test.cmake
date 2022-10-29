# This script is the main entry point for automated nightly builds of openCFS.
# It gets either called from /etc/crontab or through the provisioner scripts
# for Vagrant VBoxes.
#
# It runs all tests from hostname -> <site>/generate_test_name_list.cmake 
# a ctest -R <site>/<test>.ctest runs the single test.
#
# To call this script by do:
#
# ctest -S linux64_shared_opt_gcc_release.ctest.
#
# Simon advises the following way to simulate a crontab  environment e.g. in the following way:
#
# env -i HOME=/ PATH=/bin:/usr/bin:/sbin:/usr/sbin \
#        /opt/pckg/cmake-2.8.9/bin/ctest -S $HOME/Documents/dev/nightly_test.cmake
#
# Out of a virtual box: 
#
# vagrant ssh -c 'env -i HOME=/Users/simon PATH=/bin:/usr/bin:/sbin:/usr/sbin env && uname -a && /opt/pckg/cmake-2.8.10.2-Linux-i386/bin/cmake --version'
#
# To embed this script via crontab:
#
#  1  22   *   *   * /usr/bin/ctest -S /path/to/script/script.cmake -V \
#                                   > /dash/logs/tests.log 2>&1
#
#  "1  22   *   *   *"  specifies when to run the scheduled task. The columns correspond to minutes,
# hours, days, months, day of the week. This gets translated to 10:01 PM every day, every month.
# List and edit your crontab via crontab -l and crontab -e respectively
#
# See https://cfs.mdmt.tuwien.ac.at/trac/wiki/NightlyBuilds for more help
#
# For further information concerning CTest and CDash please refer to:
# http://www.cdash.org/cdash/resources/software.html
# http://vtk.org/Wiki/CMake_Testing_With_CTest
# http://vtk.org/Wiki/CMake_Scripting_Of_CTest

# We need at least CMake 2.8.9 for this to work
CMAKE_MINIMUM_REQUIRED(VERSION "2.8.9" FATAL_ERROR)

# List command no longer ignores empty elements.
CMAKE_POLICY(SET CMP0007 NEW)

# Get base path of current script in order to include additional macros.
GET_FILENAME_COMPONENT(CTEST_SCRIPTS_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

include("${CTEST_SCRIPTS_DIR}/../cmake_modules/CFS_macros.cmake")
# Include information about development server and further macros required for testing.
INCLUDE("${CTEST_SCRIPTS_DIR}/shared/test_macros.cmake")

# Determine date and time.
EXECUTE_PROCESS(
  COMMAND date "+%F %H:%M:%S"
  OUTPUT_VARIABLE DATE_OUT
)
STRING(REPLACE "\n" "" DATE_OUT ${DATE_OUT})
STRING(STRIP ${DATE_OUT} DATE_OUT)
headline("Starting nightly tests on ${DATE_OUT}...")

# Set global variables, e.g. path to ctest exe, site base dir, host name etc.
SET_GLOBAL_VARS()

# Set site specific variables, e.g. test user home dir, svn user password, etc.
SET_SITE_SPECIFIC_VARS()

headline("Selected global variables on ${HOSTNAME}...")

message("CMAKE_COMMAND: ${CMAKE_COMMAND}")
message("CMAKE_EXECUTABLE_SUFFIX: ${CMAKE_EXECUTABLE_SUFFIX}")
message("CMAKE_CURRENT_LIST_FILE: ${CMAKE_CURRENT_LIST_FILE}")
message("CTEST_COMMAND: ${CTEST_COMMAND}")
message("CTEST_SCRIPT_DIRECTORY: ${CTEST_SCRIPT_DIRECTORY}")
message("CTEST_SOURCE_DIRECTORY: ${CTEST_SOURCE_DIRECTORY}")
message("SITE_BASE_DIR: ${SITE_BASE_DIR}")
message("NIGHTLY_ARCHIVES_DIR: ${NIGHTLY_ARCHIVES_DIR}")
message("HOSTNAME: ${HOSTNAME}")
message("SITE_DIR: ${SITE_DIR}")
message("TESTUSER: ${TESTUSER}")
message("HOME: ${HOME}")
message("DAYOFWEEK: ${DAYOFWEEK}")
message("")

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

