#-----------------------------------------------------------------------------
# Run this script using: ctest -S /path/to/script/script.cmake -V
#
# Embedding this script into crontab:
#  1  22   *   *   * /usr/bin/ctest -S /path/to/script/script.cmake -V \
#                                   > /dash/logs/tests.log 2>&1
#
# The part of the line:
#
#  1  22   *   *   *
#
# specifies when to run the scheduled task. The columns correspond to minutes,
# hours, days, months, day of the week. This gets translated to 10:01 PM every
# day, every month.
#
# For further information concerning CTest and CDash please refer to:
# http://www.cdash.org/cdash/resources/software.html
# http://vtk.org/Wiki/CMake_Testing_With_CTest
# http://vtk.org/Wiki/CMake_Scripting_Of_CTest
# Mastering CMake (ISBN 1-930934-16-5)
#
# Notes for setting up a CDash server on openSUSE 11.x:
# - Install schema 'Web and LAMP server' in yast (contains Apache2 and MySQL)
# - Install php5-xsl module
# - Change into /srv/www/htdocs
# - Issue svn co https://www.kitware.com:8443/svn/CDash/Release-1-2-1 CDash
# - Issue chown -R wwwrun:www CDash/backup CDash/rss
# - Start and add apache and mysql services in yast runlevel editor
# - Point your web browser to http://yourhostname/CDash - Done!
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Set site name and build name
#-----------------------------------------------------------------------------
STRING(REPLACE "." ";" CFS_BUILD_HOST_LIST "${CFS_BUILD_HOST}")
LIST(GET CFS_BUILD_HOST_LIST 0 CFS_BUILD_HOST)
set(CTEST_SITE "${CFS_BUILD_HOST}")
set(CTEST_BUILD_NAME "${BUILDNAME}")

#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(ENV{LC_MESSAGES} "C")
SET(ENV{LC_ALL} "C")
SET(ENV{LANG} "C")
SET(ENV{LANGUAGE} "C")

#-----------------------------------------------------------------------------
# Set CTest command
#-----------------------------------------------------------------------------
SET(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")

#-----------------------------------------------------------------------------
# Use CMake (cmake) executable corresponding to CTest executable used to run
# this script.
#-----------------------------------------------------------------------------
SET(CTEST_CMAKE_COMMAND  "\"${CMAKE_EXECUTABLE_NAME}\"")

#-----------------------------------------------------------------------------
# Since this is the first test in the night, we have to make sure that
# the source directory is available by checking it out if the directory does
# not previously exist.
#-----------------------------------------------------------------------------

FIND_PROGRAM(CTEST_SVN_COMMAND NAMES svn)

#-----------------------------------------------------------------------------
# If the source directory does not exist perform a checkout
#-----------------------------------------------------------------------------
IF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/.svn")
  EXECUTE_PROCESS(COMMAND ${CMAKE_EXECUTABLE_NAME} -E make_directory ${CTEST_SOURCE_DIRECTORY})
  EXECUTE_PROCESS(COMMAND ${CMAKE_EXECUTABLE_NAME} -E copy ${SITE_DIR}/../../shared/CTestConfig.cmake ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)
  SET(CTEST_CHECKOUT_COMMAND "${CTEST_SVN_COMMAND} --non-interactive co ${REPO} ${CTEST_SOURCE_DIRECTORY}")
ENDIF(NOT EXISTS "${CTEST_SOURCE_DIRECTORY}/.svn")

#-----------------------------------------------------------------------------
# Either way, we have to update the working copy
#-----------------------------------------------------------------------------
SET(CTEST_UPDATE_TYPE "svn")
SET(CTEST_UPDATE_COMMAND "${CTEST_SVN_COMMAND}")
SET(CTEST_UPDATE_OPTIONS "--non-interactive ${CTEST_SOURCE_DIRECTORY}")

#-----------------------------------------------------------------------------
# Copy CDash server configuration file to source dir.
#-----------------------------------------------------------------------------
EXECUTE_PROCESS(COMMAND ${CMAKE_EXECUTABLE_NAME} -E copy_if_different ${SITE_DIR}/../../shared/CTestConfig.cmake ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)

#-----------------------------------------------------------------------------
# Start out with an empty binary directory.
#-----------------------------------------------------------------------------
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY FALSE)

CTEST_START(Nightly)
CTEST_UPDATE(SOURCE "${CTEST_SOURCE_DIRECTORY}" RETURN_VALUE res)
CTEST_SUBMIT()
