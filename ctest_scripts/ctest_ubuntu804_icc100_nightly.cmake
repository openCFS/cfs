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
# Set source and binary directories on rom
#-----------------------------------------------------------------------------
SET(CTEST_SOURCE_DIRECTORY "$ENV{HOME}/Documents/dev/NIGHTLY/CFS_TRUNK_NIGHTLY")
SET(CTEST_BINARY_DIRECTORY "$ENV{HOME}/Documents/dev/NIGHTLY/CFS_BUILD_NIGHTLY")

#-----------------------------------------------------------------------------
# Place CTestConfig.cmake file for project CFS on CDash server rom into
# CTEST_SOURCE_DIRECTORY.
#-----------------------------------------------------------------------------
FILE(WRITE "${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake"
  "
  ## This file should be placed in the root directory of your project.
  ## Then modify the CMakeLists.txt file in the root directory of your
  ## project to incorporate the testing dashboard.
  ## # The following are required to uses Dart and the Cdash dashboard
  # ENABLE_TESTING()
  # INCLUDE(Dart)
  set(CTEST_PROJECT_NAME \"CFS\")
  set(CTEST_NIGHTLY_START_TIME \"00:00:00 CET\")

  set(CTEST_DROP_METHOD \"http\")
  set(CTEST_DROP_SITE \"rom\")
  set(CTEST_DROP_LOCATION \"/cdash/submit.php?project=CFS\")
  set(CTEST_DROP_SITE_CDASH TRUE)"
)

#-----------------------------------------------------------------------------
# Specify that we want to do an experimental build without updating the CFS++
# working copy. I.e. we leave away the ExperimentalUpdate step between
# ExperimentalStart and ExperimentalConfigure. The Subversion update gets
# done by simon on mac before the the test scripts on rom get executed.
#-----------------------------------------------------------------------------
SET(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")
SET(CTEST_COMMAND "${CTEST_COMMAND} -D NightlyStart")
SET(CTEST_COMMAND "${CTEST_COMMAND} -D NightlyConfigure")
SET(CTEST_COMMAND "${CTEST_COMMAND} -D NightlyBuild")
SET(CTEST_COMMAND "${CTEST_COMMAND} -D NightlyTest")
SET(CTEST_COMMAND "${CTEST_COMMAND} -D NightlySubmit")
SET(CTEST_COMMAND "${CTEST_COMMAND} -A ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt")
SET(CTEST_COMMAND "${CTEST_COMMAND} -E optimization")
#SET(CTEST_COMMAND "${CTEST_COMMAND} -R torque3d")

#-----------------------------------------------------------------------------
# Use CMake (cmake) executable corresponding to CTest executable used to run
# this script.
#-----------------------------------------------------------------------------
SET(CTEST_CMAKE_COMMAND  "\"${CMAKE_EXECUTABLE_NAME}\"")

#-----------------------------------------------------------------------------
# Start out with an empty binary directory.
#-----------------------------------------------------------------------------
SET(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY TRUE)

#-----------------------------------------------------------------------------
# Enter the following values into the initial cache
# ${CTEST_BINARY_DIRECTORY}/CMakeCache.txt before starting the configure
# run.
#-----------------------------------------------------------------------------
SET(CTEST_INITIAL_CACHE
  "BUILD_TESTING:BOOL=ON
   DEBUG:BOOL=OFF
   TESTSUITE_DIR:STRING=$ENV{HOME}/Documents/dev/NIGHTLY/CFS_TESTSUITE_NIGHTLY
   CFS_DEPS_ROOT:PATH=$ENV{HOME}/Documents/dev/NIGHTLY/CFSDEPS_NIGHTLY
   CFS_DEPS_CACHE_DIR:PATH=$ENV{HOME}/Documents/dev/NIGHTLY/CFSDEPSCACHE
   MKL_ROOT_DIR:PATH=/opt/intel/mkl/10.0.5.025
   CFS_PARDISO:STRING=MKL
   USE_GMV_INPUT:BOOL=ON
   USE_SCPIP:BOOL=ON")

#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(CTEST_ENVIRONMENT
  "CC=icc"
  "CXX=icpc"
  "FC=ifort"
  "LC_MESSAGES=C"
  "LC_ALL=C"
  "LANG=C"
  "LANGUAGE=C"
  )
