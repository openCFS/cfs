#!/bin/sh

# This script is run by the following entry in /etc/crontab on rom
#
# 30    1 * * *   strieben /home/strieben/Documents/dev/nightly_cfs_test.sh > /home/strieben/Documents/dev/nightly_cfs_test.log 2>&1


echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

# Set path to home directory of strieben, just to make sure...
export HOME=/home/strieben

# Set path to strieben's development directory
DEVDIR="$HOME/Documents/dev"

# Set path to directory for nightly builds
TESTDIR="$DEVDIR/NIGHTLY"

# The generated release executables should end up here.
# Directory is writeable by strieben.
DESTDIR=/opt/pckg/cfs_nightly

# Set ssh environment for Subversion
export SVN_SSH="ssh -i $DEVDIR/svn_id_dsa_cfstestuser_klu -l svn -p 22022"

# Extend system path by /opt/pckg/bin
export PATH=/opt/pckg/bin:$PATH

# Create nightly build directory
mkdir -p $TESTDIR
cd $TESTDIR

# Remove previous CFSDEPSCACHE directory
rm -rf $TESTDIR/CFSDEPSCACHE
mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml

# Either checkout out or update CFS++ trunk working copy. Subversion update
# had some strange problems if performed by CTest.
TRUNK_REPO="svn+ssh://131.188.140.10/software/CFS++/trunk"
if [ ! -d "CFS_TRUNK_NIGHTLY/.svn" ]; then
    svn co $TRUNK_REPO CFS_TRUNK_NIGHTLY
else
    svn up $TRUNK_REPO CFS_TRUNK_NIGHTLY
fi

# Account for 90 seconds reconnect timeout of the lse10 firewall.
sleep 120

# Either checkout out or update CFSDEPS trunk working copy.
DEPS_REPO="svn+ssh://131.188.140.10/software/cfsdeps/trunk"
if [ ! -d "CFSDEPS_NIGHTLY/.svn" ]; then
    svn co $DEPS_REPO CFSDEPS_NIGHTLY
else
    svn up CFSDEPS_NIGHTLY
fi

# Copy sources of SCPIP to CFSDEPS_NIGHTLY/scpip
cp $DEVDIR/scpip.tar.bz2 CFSDEPS_NIGHTLY/scpip

# Copy whole CFSDEPS tree to /opt/pckg/CFSDEPS. This should make sure that
# we allways have a current CFSDEPS present on the machine.
rm -rf /opt/pckg/CFSDEPS/* /opt/pckg/CFSDEPS/.svn
cp -a CFSDEPS_NIGHTLY/* CFSDEPS_NIGHTLY/.svn /opt/pckg/CFSDEPS

# Account for 90 seconds reconnect timeout of the lse10 firewall.
sleep 120

# Either checkout out or update CFS++ test suite trunk working copy.
SUITE_REPO="svn+ssh://131.188.140.10/software/CFS++_TEST/trunk"
if [ ! -f "CFS_TESTSUITE_NIGHTLY/CMakeLists.txt" ]; then
    svn co $SUITE_REPO CFS_TESTSUITE_NIGHTLY
else
    svn up CFS_TESTSUITE_NIGHTLY
fi

# Pack directories for Ubuntu VBox
tar cvzf CFS_NIGHTLY.tgz CFS_TESTSUITE_NIGHTLY CFSDEPS_NIGHTLY CFS_TRUNK_NIGHTLY

# Start Ubuntu 8.04 VBox in background
VBOX_LOG=$HOME/Documents/dev/ubuntu804.log
echo "Ubuntu 8.04 VBox started on `date`" > $VBOX_LOG
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo >> $VBOX_LOG
VBoxHeadless -startvm "Ubuntu 8.04"  >> $VBOX_LOG 2>&1 &
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo "Ubuntu 8.04 VBox started in background on `date`" >> $VBOX_LOG
echo >> $VBOX_LOG

# Change into CFS++ source directory and execute CTest for GCC.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_rom_gcc_nightly.cmake
if [ -f $TESTDIR/CFS_BUILD_NIGHTLY/bin/cfs ]; then
    rm -rf $DESTDIR/trunk_gcc
    mkdir $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_gcc
fi

# Change into CFS++ source directory and execute CTest for ICC.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_rom_icc_nightly.cmake
if [ -f $TESTDIR/CFS_BUILD_NIGHTLY/bin/cfs ]; then
    rm -rf $DESTDIR/trunk_icc
    mkdir $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_icc
fi

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
