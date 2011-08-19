#!/bin/bash

# This script is run by the last line in /etc/crontab
# on the Ubuntu 8.04 VBox. The first two lines make sure
# that the shared folders for the nightly sources from rom
# and the CFD data from qnap get mounted inside the VBox.
# Please note: It is absolutely necessary to have the VBox
# guest additions installed for this to properly work!
#
# 0  2    * * *   root    /sbin/mount.vboxsf -o uid=1000,gid=1000 dev /mnt > /home/user/Documents/dev/mount_share.log 2>&1
# 10 2    * * *   root    /sbin/mount.vboxsf CFD_Data /media/CFD_Data >> /home/user/Documents/dev/mount_share.log 2>&1
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

rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
mv $TESTDIR/CFSDEPSCACHE $DEVDIR

# Remove and recreate nightly build directory
rm -rf $TESTDIR
mkdir -p $TESTDIR
cd $TESTDIR

mv $DEVDIR/CFSDEPSCACHE $TESTDIR

# Unpack archive of updated working copies
tar xvzf "/mnt/NIGHTLY/CFS_NIGHTLY.tgz"

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
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_ICC91 ]; then
    mv CFS_BUILD_ICC91 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_icc91_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_ICC91

. /opt/intel/cc/10.0.025/bin/iccvars.sh
. /opt/intel/fc/10.0.025/bin/ifortvars.sh
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_ICC100 ]; then
    mv CFS_BUILD_ICC100 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_icc100_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_ICC100

rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_GCC41 ]; then
    mv CFS_BUILD_GCC41 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_gcc41_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_GCC41

# The following Intel compilers (10.1 and 11.0) support GCC 4.2.
# Note however that Intel C/C++ 10.1 does not support GCC 4.3!
rm -rf $TESTDIR/bin
mkdir $TESTDIR/bin
ln -s /usr/bin/gcc-4.2 $TESTDIR/bin/gcc
ln -s /usr/bin/g++-4.2 $TESTDIR/bin/g++
ln -s /usr/bin/gfortran-4.2 $TESTDIR/bin/gfortran

rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_GCC42 ]; then
    mv CFS_BUILD_GCC42 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_gcc_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_GCC42

. /opt/intel/cc/10.1.018/bin/iccvars.sh
. /opt/intel/fc/10.1.018/bin/ifortvars.sh
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_ICC101 ]; then
    mv CFS_BUILD_ICC101 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_icc101_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_ICC101

. /opt/intel/Compiler/11.0/081/bin/iccvars.sh ia32
. /opt/intel/Compiler/11.0/081/bin/ifortvars.sh ia32
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_ICC110 ]; then
    mv CFS_BUILD_ICC110 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_icc110_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_ICC110

rm -rf $TESTDIR/bin
mkdir $TESTDIR/bin
ln -s /usr/bin/gcc-3.4 $TESTDIR/bin/gcc
ln -s /usr/bin/g++-3.4 $TESTDIR/bin/g++
ln -s /usr/bin/g77 $TESTDIR/bin/g77
rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
cd $TESTDIR
if [ -d CFS_BUILD_NIGHTLY ]; then
    rm -rf CFS_BUILD_NIGHTLY
fi
if [ -d CFS_BUILD_GCC34 ]; then
    mv CFS_BUILD_GCC34 CFS_BUILD_NIGHTLY
fi
# Change into CFS++ trunk source dir.
cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_ubuntu804_gcc34_nightly.cmake
cd $TESTDIR && mv CFS_BUILD_NIGHTLY CFS_BUILD_GCC34

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
