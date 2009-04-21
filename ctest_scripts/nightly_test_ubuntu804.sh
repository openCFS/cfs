#!/bin/bash

# This script is run by the following entry in /etc/crontab
# on the Ubuntu 8.04 VBox
#
# 30 2    * * *   user    /home/user/Documents/dev/nightly_cfs_test.sh > /home/user/Documents/dev/nightly_cfs_test.log 2>&1

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

# Set path to home directory of user, just to make sure...
export HOME=/home/user

# Set path to strieben's development directory
DEVDIR="$HOME/Documents/dev"

# Set path to directory for nightly builds
TESTDIR="$DEVDIR/NIGHTLY"

# Extend system path by /opt/pckg/bin
export PATH=/opt/pckg/bin:$PATH

# Check if the archive provided by rom is present on the VBox share
if [ ! -f "/mnt/NIGHTLY/CFS_NIGHTLY.tgz" ]; then
    echo "Cannot find TGZ archive on share!"
    exit 1
fi

# Remove and recreate nightly build directory
rm -rf $TESTDIR
mkdir -p $TESTDIR
cd $TESTDIR

# Unpack archive of updated working copies
tar xvzf "/mnt/NIGHTLY/CFS_NIGHTLY.tgz"

# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY

# The following Intel compilers (9.1 and 10.0) only support GCC up to
# version 4.1. Therefore we create links to the proper compilers and
# make sure that they appear first on path.
rm -rf $TESTDIR/bin
mkdir $TESTDIR/bin
ln -s /usr/bin/gcc-4.1 $TESTDIR/bin/gcc
ln -s /usr/bin/g++-4.1 $TESTDIR/bin/g++
ln -s /usr/bin/gfortran-4.1 $TESTDIR/bin/gfortran
export PATH=$TESTDIR/bin:$PATH

. /opt/intel/cc/9.1.043/bin/iccvars.sh
. /opt/intel/fc/9.1.043/bin/ifortvars.sh
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
ctest -V -S ctest_scripts/ctest_ubuntu804_icc91_nightly.cmake

. /opt/intel/cc/10.0.025/bin/iccvars.sh
. /opt/intel/fc/10.0.025/bin/ifortvars.sh
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
ctest -V -S ctest_scripts/ctest_ubuntu804_icc100_nightly.cmake

ctest -V -S ctest_scripts/ctest_ubuntu804_gcc41_nightly.cmake

# The following Intel compilers (10.1 and 11.0) support GCC 4.2.
# Note however that Intel C/C++ 10.1 does not support GCC 4.3!
rm -rf $TESTDIR/bin
mkdir $TESTDIR/bin
ln -s /usr/bin/gcc-4.2 $TESTDIR/bin/gcc
ln -s /usr/bin/g++-4.2 $TESTDIR/bin/g++
ln -s /usr/bin/gfortran-4.2 $TESTDIR/bin/gfortran

ctest -V -S ctest_scripts/ctest_ubuntu804_gcc_nightly.cmake

. /opt/intel/cc/10.1.018/bin/iccvars.sh
. /opt/intel/fc/10.1.018/bin/ifortvars.sh
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
ctest -V -S ctest_scripts/ctest_ubuntu804_icc101_nightly.cmake

. /opt/intel/Compiler/11.0/081/bin/iccvars.sh ia32
. /opt/intel/Compiler/11.0/081/bin/ifortvars.sh ia32
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
ctest -V -S ctest_scripts/ctest_ubuntu804_icc110_nightly.cmake

rm -rf $TESTDIR/bin
mkdir $TESTDIR/bin
ln -s /usr/bin/gcc-3.4 $TESTDIR/bin/gcc
ln -s /usr/bin/g++-3.4 $TESTDIR/bin/g++
ln -s /usr/bin/g77 $TESTDIR/bin/g77
ctest -V -S ctest_scripts/ctest_ubuntu804_gcc34_nightly.cmake

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
