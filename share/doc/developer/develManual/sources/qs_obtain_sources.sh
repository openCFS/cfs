# Create directory structure for CFS++ development
mkdir -p $HOME/Documents/dev
cd $HOME/Documents/dev

# Set variable for repository base path
REPO_BASE=https://lse17.e-technik.uni-erlangen.de:2001/svn

# Set variable with your user name
USER=PLACE_YOUR_USER_NAME_HERE

# Checkout trunk of CFS++ to local directory CFS_TRUNK
svn checkout --username $USER $REPO_BASE/CFS++/trunk CFS_TRUNK

# Checkout trunk of test suite to local directory CFS_TESTSUITE (optional)
svn checkout --username $USER $REPO_BASE/CFS++_TEST/trunk CFS_TESTSUITE
