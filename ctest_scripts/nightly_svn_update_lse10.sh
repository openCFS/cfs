#!/bin/sh

# This script is run by the following entry in /etc/crontab on lse10
#
# 0     1 * * *   root /home/data/programs/cfs/nightly_svn_update_lse10.sh > /home/data/programs/cfs/nightly_svn_update_lse10.sh 2>&1

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

# Set path to home directory of strieben, just to make sure...
export HOME=/home/users/strieben

# Set path to strieben's development directory
DEVDIR="/home/data/programs/cfs"

# Set path to directory for nightly builds
TESTDIR="$DEVDIR/nightly"

# The generated release executables should end up here.
# Directory is writeable by strieben.
#DESTDIR=/opt/pckg/cfs_nightly

# Set ssh environment for Subversion
#export SVN_SSH="ssh -i $DEVDIR/svn_id_dsa_cfstestuser_klu -l svn -p 22022"

# Extend system path by /opt/pckg/bin
#export PATH=/opt/pckg/bin:$PATH

# Path where Pardiso (Schenk) license file may be found
# export PARDISO_LIC_PATH="$HOME/Documents/dev/pardiso/keys/sedici"

# Determine day of week since we only want to do a full rebuild on sunday (7)
DAYOFWEEK=$(date +%u)
echo "DAYOFWEEK $DAYOFWEEK"

# Create nightly build directory
mkdir -p $TESTDIR
cd $TESTDIR

function UpdateWorkingCopy {
    REPOS=$1
    WCDIR=$2
    
    # Get previous revision for CFSDEPS
    if [ -f CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt ]; then
        PREV_REV=$(svn st --xml --verbose $WCDIR | xsltproc CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    fi
    if [ ! "$?" = 0 ]; then PREV_REV="-1"; fi

    if [ ! -d "$WCDIR/.svn" ]; then
	svn co $REPOS $WCDIR
    else
	svn up $WCDIR
    fi

    # Get current revision for CFS++
    if [ -f CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt ]; then
        CURRENT_REV=$(svn st --xml --verbose $WCDIR | xsltproc CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    fi
    if [ ! "$?" = 0 ]; then CURRENT_REV="-1"; fi
}

# Either checkout out or update CFS++ trunk working copy. Subversion update
# had some strange problems if performed by CTest.
TRUNK_REPO="file:///home/data/svn-master/repositories/software/CFS++/trunk"
if [ ! -d CFS_TRUNK_NIGHTLY ]; then
    UpdateWorkingCopy "$TRUNK_REPO" "CFS_TRUNK_NIGHTLY"
    CFS_PREV_REV=$PREV_REV;
    CFS_CURRENT_REV=$CURRENT_REV;
    echo "CFS_PREV_REV $CFS_PREV_REV"
    echo "CFS_CURRENT_REV $CFS_CURRENT_REV"
fi 

DEPS_REPO="file:///home/data/svn-master/repositories/software/cfsdeps/trunk"
UpdateWorkingCopy "$DEPS_REPO" "CFSDEPS_NIGHTLY"
CFSDEPS_PREV_REV=$PREV_REV;
CFSDEPS_CURRENT_REV=$CURRENT_REV;
echo "CFSDEPS_PREV_REV $CFSDEPS_PREV_REV"
echo "CFSDEPS_CURRENT_REV $CFSDEPS_CURRENT_REV"

# Either checkout out or update CFS++ test suite trunk working copy.
SUITE_REPO="file:///home/data/svn-master/repositories/software/CFS++_TEST/trunk"
UpdateWorkingCopy "$SUITE_REPO" "CFS_TESTSUITE_NIGHTLY"
TESTSUITE_PREV_REV=$PREV_REV;
TESTSUITE_CURRENT_REV=$CURRENT_REV;
echo "TESTSUITE_PREV_REV $TESTSUITE_PREV_REV"
echo "TESTSUITE_CURRENT_REV $TESTSUITE_CURRENT_REV"

mkdir -p CFS_BUILD_NIGHTLY
if [ -d CFS_BUILD_NIGHTLY/Testing ]; then
    rm -rf CFS_BUILD_NIGHTLY/Testing
fi

# Update CFS++ source the CTest way.
ctest -VV -S CFS_TRUNK_NIGHTLY/ctest_scripts/ctest_lse10_start_update.cmake

tar cvzf nightly.tgz CFSDEPS_NIGHTLY CFS_TRUNK_NIGHTLY CFS_TESTSUITE_NIGHTLY CFS_BUILD_NIGHTLY/Testing
