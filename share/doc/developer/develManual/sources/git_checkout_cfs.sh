# Create directory structure for CFS++ development
mkdir -p $HOME/Documents/dev
cd $HOME/Documents/dev

# Set variable for repository base path
REPO_BASE=@CFS_DS_SVN@

# Fetch complete history of CFS++ trunk SVN repository (may take a long time)
git svn clone $REPO_BASE/CFS++/trunk CFS_TRUNK_GIT

# Make some changes to sources
...

# Stash local changes
git stash

# Update local repository with commits from central repository
git svn rebase

# Apply stashed changes
git stash apply

# Make a local commit
git commit -a -m "My local commit"

# Commit all local commits to central SVN repository
git svn dcommit
