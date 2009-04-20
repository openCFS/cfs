#!/bin/bash

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

HOME=/home/user
TGZDIR="/mnt/Documents/dev"
DEVDIR="$HOME/Documents/dev"
TESTDIR="$DEVDIR/NIGHTLY"

if [ ! -f "$TGZDIR/CFS_TRUNK_NIGHTLY.tgz" ]; then
    echo "Share sedici/strieben is not mounted!"
    exit 1
fi

. $HOME/.bashrc

mkdir -p $TESTDIR
cd $TESTDIR

. $TGZDIR/CFS_NIGHTLY_BUILD.sh
export PATH=/opt/pckg/bin:$PATH

echo "CFS $CFS"
echo "CFSDEPS $CFSDEPS"
echo "TESTSUITE $TESTSUITE"

if [ "$CFS" = "OK" ]; then
    tar xvzf $TGZDIR/CFS_TRUNK_NIGHTLY.tgz
fi

if [ "$CFSDEPS" = "OK" ]; then
    tar xvzf $TGZDIR/CFSDEPS_NIGHTLY.tgz
fi

if [ "$TESTSUITE" = "OK" ]; then
    tar xvzf $TGZDIR/CFS_TESTSUITE_NIGHTLY.tgz
fi

rm -rf $TESTDIR/CFSDEPSCACHE

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
