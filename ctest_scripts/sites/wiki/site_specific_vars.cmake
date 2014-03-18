MESSAGE(
"
=============================================================================
 Setting site-specific variables on ${HOSTNAME}...
=============================================================================
"
)

# Make HOME also an environment variable. It is necessary since HOME is set
# to / by /etc/crontab. Setting the HOME environment variable is possible
# since this script is run by CTest instead of CMake.
SET(ENV{HOME} ${HOME})

# Set the path to the Vagrant config folder to a path on a different
# partition since Vagrant unpacks temporary files there, which may become
# very large. The testuser home partition has very little space left!
SET(ENV{VAGRANT_HOME} "/var/local/testuser/.vagrant.d")

# User name of test user for SVN
SET(CFS_TESTUSER testuser-klu)

# Read testuser PW for Subversion server from ${HOME}/.testuserpw
FILE(READ "${HOME}/.testuserpw" CFS_TESTUSER_PW)
STRING(REPLACE "=" ";" CFS_TESTUSER_PW_LIST ${CFS_TESTUSER_PW})
LIST(GET CFS_TESTUSER_PW_LIST 1 CFS_TESTUSER_PW)
STRING(REPLACE "\"" "" CFS_TESTUSER_PW ${CFS_TESTUSER_PW})
STRING(STRIP ${CFS_TESTUSER_PW} CFS_TESTUSER_PW)

# Directory for nightly CFS++ binaries
SET(CFS_NIGHTLY_DIR "/opt/pckg/cfs_nightly")

# Directory for nightly CFS++ binaries on web server
SET(APACHE_NIGHTLY_DIR "/var/www/html/cfs_nightly")

# Directory for nightly documenation on web server
SET(CFS_DOCU_DIR "/var/www/html/cfs_docu")

