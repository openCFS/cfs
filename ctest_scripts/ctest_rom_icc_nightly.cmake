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
# Specify that we want to do an experimental build without updating the CFS++
# working copy. I.e. we leave away the ExperimentalUpdate step between
# ExperimentalStart and ExperimentalConfigure. The Subversion update gets
# done by simon on mac before the the test scripts on rom get executed.
#-----------------------------------------------------------------------------
SET(CTEST_COMMAND  "\"${CTEST_EXECUTABLE_NAME}\"")

#-----------------------------------------------------------------------------
# Use CMake (cmake) executable corresponding to CTest executable used to run
# this script.
#-----------------------------------------------------------------------------
SET(CTEST_CMAKE_COMMAND  "\"${CMAKE_EXECUTABLE_NAME}\"")

#-----------------------------------------------------------------------------
# Set the following environment variables for the test run. This can be used
# to specifiy the compilers and that all messages should be output in English
# language, so that CTest may properly parse them.
#-----------------------------------------------------------------------------
SET(ENV{CC} "icc")
SET(ENV{CXX} "icpc")
SET(ENV{FC} "ifort")
SET(ENV{LC_MESSAGES} "C")
SET(ENV{LC_ALL} "C")
SET(ENV{LANG} "C")
SET(ENV{LANGUAGE} "C")
SET(ENV{CPLREADER_PERF_SUITE} "/media/CFD_Data/cplreader_performance_suite")

EXEC_PROGRAM("rm -rf ${CTEST_BINARY_DIRECTORY} && mkdir -p ${CTEST_BINARY_DIRECTORY}/CMakeFiles && perl ${CTEST_SOURCE_DIRECTORY}/share/scripts/identify_compiler.pl $ENV{CXX} ${CTEST_SOURCE_DIRECTORY}/share/scripts/IdentifyCXXCompiler.cpp cmake > ${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.cmake"
  ARGS
  OUTPUT_VARIABLE CXX_COMPILER_INFO
  RETURN_VALUE RETVAL)

INCLUDE(${CTEST_BINARY_DIRECTORY}/CMakeFiles/out.cmake)

EXEC_PROGRAM("${CTEST_SOURCE_DIRECTORY}/share/scripts/distro.sh -u"
  ARGS
  OUTPUT_VARIABLE ARCH_STR
  RETURN_VALUE RETVAL)

#MESSAGE("${CXX_ID} ${CXX_VERSION} ${CXX_GCC_VERSION}")
#MESSAGE("${ARCH_STR}")

SET(BUILDTYPE "RELEASE")
SET(BUILDNAME "${ARCH_STR} ${CXX_ID} ${CXX_VERSION} ${BUILDTYPE}")

SITE_NAME(CFS_BUILD_HOST)
set(CTEST_SITE "${CFS_BUILD_HOST}")
set(CTEST_BUILD_NAME "${BUILDNAME}")

#-----------------------------------------------------------------------------
# Copy CDash server configuration file to source dir.
#-----------------------------------------------------------------------------
EXECUTE_PROCESS(COMMAND ${CMAKE_EXECUTABLE_NAME} -E copy_if_different CTestConfig.cmake ${CTEST_SOURCE_DIRECTORY}/CTestConfig.cmake)

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
   CMAKE_BUILD_TYPE:STRING=${BUILDTYPE}   
   TESTSUITE_DIR:STRING=$ENV{HOME}/Documents/dev/NIGHTLY/CFS_TESTSUITE_NIGHTLY
   CFS_DEPS_ROOT:PATH=$ENV{HOME}/Documents/dev/NIGHTLY/CFSDEPS_NIGHTLY
   CFS_DEPS_CACHE_DIR:PATH=$ENV{HOME}/Documents/dev/NIGHTLY/CFSDEPSCACHE
   USE_GMV_INPUT:BOOL=ON
   USE_GMSH:BOOL=ON   
   USE_ANSYSRST:BOOL=ON
   USE_PYTHON:BOOL=ON
   USE_TCL:BOOL=ON
   USE_INTERPOLATION:BOOL=ON
   CPLREADER:BOOL=ON
   CPLREADER_OPENFOAM:BOOL=ON
   CPLREADER_CGNS:BOOL=ON
   CPLREADER_CFX:BOOL=ON
   CPLREADER_FIELDVIEW:BOOL=ON
   CPLREADER_ENSIGHT:BOOL=ON
   CPLREADER_FLUENT:BOOL=ON   
   USE_SCPIP:BOOL=ON
   USE_CHOLMOD:BOOL=ON
   BUILDNAME:STRING=${BUILDNAME}")

SET(CTEST_CMAKE_GENERATOR "Unix Makefiles")

message("Start dashboard...")
ctest_start(Nightly)
FILE(WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" ${CTEST_INITIAL_CACHE})
message("  Configure")
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
# Configure twice so that CFX and OpenFOAM dependencies get built
ctest_configure(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
#message("read custom files after configure")
#ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
message("  Build")
ctest_build(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Test")
ctest_test(BUILD "${CTEST_BINARY_DIRECTORY}" RETURN_VALUE res)
message("  Submit")
ctest_submit(RETURN_VALUE res)
message("  All done")
