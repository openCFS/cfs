#!/bin/sh

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

export HOME=/home/users/strieben
TGZDIR="$HOME/Documents/dev"
DEVDIR="$HOME/Documents/dev"
TESTDIR="$DEVDIR/NIGHTLY"
DESTDIR=/opt/pckg/cfs_nightly
# export PARDISO_LIC_PATH="$HOME/Documents/dev/pardiso/keys/sedici"

. $HOME/.bashrc

mkdir -p $TESTDIR
cd $TESTDIR

. $TGZDIR/CFS_NIGHTLY_BUILD.sh

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
mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml
cp /opt/pckg/CFSDEPSCACHE/sources/acml/*.tgz $TESTDIR/CFSDEPSCACHE/sources/acml

cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_sedici_gcc_nightly.cmake
if [ -f $TESTDIR/CFS_BUILD_NIGHTLY/bin/cfs ]; then
    rm -rf $DESTDIR/trunk_gcc
    mkdir $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_gcc
fi
#cd $TESTDIR/CFS_BUILD_NIGHTLY && make NightlyTest
#make NightlySubmit

cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_sedici_icc_nightly.cmake
if [ -f $TESTDIR/CFS_BUILD_NIGHTLY/bin/cfs ]; then
    rm -rf $DESTDIR/trunk_icc 
    mkdir $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_icc
fi
#cd $TESTDIR/CFS_BUILD_NIGHTLY && make NightlyTest
#make NightlySubmit

cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_sedici_gcc_schenk_nightly.cmake
ctest -V -S ctest_scripts/ctest_sedici_gcc_acml_nightly.cmake
ctest -V -S ctest_scripts/ctest_sedici_gcc_gotoblas_nightly.cmake

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
