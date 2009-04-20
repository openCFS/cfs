#!/bin/sh

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

export HOME=/home/strieben
TGZDIR="$HOME/Documents/strieben/Documents/dev"
DEVDIR="$HOME/Documents/dev"
TESTDIR="$DEVDIR/NIGHTLY"
DESTDIR=/opt/pckg/cfs_nightly

if [ ! -f "$TGZDIR/CFS_TRUNK_NIGHTLY.tgz" ]; then
    echo "Share sedici/strieben is not mounted!"
    exit 1
fi

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

# Start Ubuntu 8.04 VBox
VBOX_LOG=$HOME/Documents/dev/ubuntu804.log
echo "Ubuntu 8.04 VBox started on `date`" > $VBOX_LOG
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo >> $VBOX_LOG
VBoxHeadless -startvm "Ubuntu 8.04"  >> $VBOX_LOG 2>&1 &
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo "Ubuntu 8.04 VBox started in background on `date`" >> $VBOX_LOG
echo >> $VBOX_LOG


cd $TESTDIR/CFS_TRUNK_NIGHTLY
ctest -V -S ctest_scripts/ctest_rom_gcc_nightly.cmake
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
ctest -V -S ctest_scripts/ctest_rom_icc_nightly.cmake
if [ -f $TESTDIR/CFS_BUILD_NIGHTLY/bin/cfs ]; then
    rm -rf $DESTDIR/trunk_icc
    mkdir $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_icc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_icc
fi
#cd $TESTDIR/CFS_BUILD_NIGHTLY && make NightlyTest
#make NightlySubmit

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
