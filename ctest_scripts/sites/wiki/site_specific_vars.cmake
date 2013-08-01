MESSAGE(
"
=============================================================================
 Setting site-specific variables on ${HOSTNAME}...
=============================================================================
"
)

# Set user name
SET(TESTUSER testuser)

# Set home directory
SET(HOME /home/testuser)

# Make HOME also an environment variable. This is possible since this script
# is run by CTest instead of CMake.
SET(ENV{HOME} ${HOME})

# User name of test user for SVN
SET(CFS_TESTUSER testuser-klu)

# Read testuser PW for Subversion server from ${HOME}/.testuserpw
FILE(READ "${HOME}/.testuserpw" CFS_TESTUSER_PW)
STRING(REPLACE "=" ";" CFS_TESTUSER_PW_LIST ${CFS_TESTUSER_PW})
LIST(GET CFS_TESTUSER_PW_LIST 1 CFS_TESTUSER_PW)
STRING(REPLACE "\"" "" CFS_TESTUSER_PW ${CFS_TESTUSER_PW})
STRING(STRIP ${CFS_TESTUSER_PW} CFS_TESTUSER_PW)
