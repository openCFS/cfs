# Set user name
SET(TESTUSER testuser)

# Set home directory
SET(HOME /home/testuser)

# User name of test user for SVN
SET(CFS_TESTUSER testuser-klu)

# Read testuser PW for Subversion server from ${HOME}/.testuserpw
FILE(READ "${HOME}/.testuserpw" CFS_TESTUSER_PW)
STRING(REPLACE "=" ";" CFS_TESTUSER_PW_LIST ${CFS_TESTUSER_PW})
LIST(GET CFS_TESTUSER_PW_LIST 1 CFS_TESTUSER_PW)
STRING(REPLACE "\"" "" CFS_TESTUSER_PW ${CFS_TESTUSER_PW})
STRING(STRIP ${CFS_TESTUSER_PW} CFS_TESTUSER_PW)
