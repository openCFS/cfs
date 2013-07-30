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
SITE_NAME(CFS_BUILD_HOST)
SET(BUILDNAME "Update CFSDEPS Trunk")

#-----------------------------------------------------------------------------
# Set source and binary directories on rom/sedici
#-----------------------------------------------------------------------------
SET(CTEST_SOURCE_DIRECTORY "${HOME}/Documents/dev/NIGHTLY/CFSDEPS_TRUNK_NIGHTLY")
SET(CTEST_BINARY_DIRECTORY "${HOME}/Documents/dev/NIGHTLY/CFSDEPS_TRUNK_NIGHTLY")

SET(REPO "https://lse17.e-technik.uni-erlangen.de:2001/svn/cfsdeps/trunk")
SET(USER ${CFS_TESTUSER})

INCLUDE(${SITE_DIR}/ctest_update.cmake)
