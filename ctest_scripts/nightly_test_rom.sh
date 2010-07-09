#!/bin/sh

# This script is run by the following entry in /etc/crontab on rom
#
# 30    1 * * *   strieben /home/strieben/Documents/dev/nightly_cfs_test.sh > /home/strieben/Documents/dev/nightly_cfs_test.log 2>&1

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

# Set path to home directory of strieben, just to make sure...
export HOME=/home/users/strieben

# Set path to strieben's development directory
DEVDIR="$HOME/Documents/dev"

# Set path to directory for nightly builds
TESTDIR="$DEVDIR/NIGHTLY"

# The generated release executables should end up here.
# Directory is writeable by strieben.
DESTDIR=/opt/pckg/cfs_nightly

# Extend system path by /opt/pckg/bin
export PATH=/opt/pckg/bin:$PATH

# Determine day of week since we only want to do a full rebuild on sunday (7)
DAYOFWEEK=$(date +%u)
echo "DAYOFWEEK $DAYOFWEEK"

# Create nightly build directory
mkdir -p $TESTDIR
cd $TESTDIR

function GetWorkingCopyRev {
    WCDIR=$1
    
    # Get previous revision for CFSDEPS
    WC_REV=$(svn st --xml --verbose $WCDIR | xsltproc $TESTDIR/CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    if [ ! "$?" = 0 ]; then WC_REV="-1"; fi

    change-svn-wc-format.py $WCDIR "1.5" # --verbose
}

function PerformTest {
    TESTNAME=$1;

    cd $TESTDIR/CFS_TRUNK_NIGHTLY

    cat ctest_scripts/ctest_$TESTNAME.cmake > ctest_scripts/ctest_script.cmake

    ctest -V -S ctest_scripts/ctest_script.cmake

    cd $TESTDIR
    tar cvzf cfs_build_$TESTNAME.tgz CFS_BUILD_NIGHTLY
}


# Remove previous CFSDEPSCACHE directory and remove prebuilt archives
if [ "$DAYOFWEEK" = "7" ]; then
    rm -rf $TESTDIR/CFSDEPSCACHE
    rm -rf $TESTDIR/cfs_build_*.tgz
fi

# Checkout or update CFSDEPS
cd $TESTDIR/CFS_TRUNK_NIGHTLY
GetWorkingCopyRev $TESTDIR/CFSDEPS_NIGHTLY
CFSDEPS_PREV_REV=$WC_REV;
ctest -V -S ctest_scripts/ctest_rom_update_cfsdeps.cmake
GetWorkingCopyRev $TESTDIR/CFSDEPS_NIGHTLY
CFSDEPS_CURRENT_REV=$WC_REV;

# Checkout or update test suite
cd $TESTDIR/CFS_TRUNK_NIGHTLY
GetWorkingCopyRev $TESTDIR/CFS_TESTSUITE_NIGHTLY
TESTSUITE_PREV_REV=$WC_REV;
ctest -V -S ctest_scripts/ctest_rom_update_cfsdeps.cmake
GetWorkingCopyRev $TESTDIR/CFS_TESTSUITE_NIGHTLY
TESTSUITE_CURRENT_REV=$WC_REV;

# Remove cache directory if CFSDEPS have changed
if [ ! "$CFSDEPS_CURRENT_REV" = "$CFSDEPS_PREV_REV" ]; then
    rm -rf $TESTDIR/CFSDEPSCACHE
    mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml
    cp /opt/pckg/CFSDEPSCACHE/sources/acml/*.tgz $TESTDIR/CFSDEPSCACHE/sources/acml
fi

cd $TESTDIR

# Copy sources of SCPIP to CFSDEPS_NIGHTLY/scpip
# cp $DEVDIR/scpip.tar.bz2 CFSDEPS_NIGHTLY/scpip

# Copy whole CFSDEPS tree to /opt/pckg/CFSDEPS. This should make sure that
# we allways have a current CFSDEPS present on the machine.
rm -rf /opt/pckg/CFSDEPS/* /opt/pckg/CFSDEPS/.svn
cp -a CFSDEPS_NIGHTLY/* CFSDEPS_NIGHTLY/.svn /opt/pckg/CFSDEPS

# Pack directories for Ubuntu VBox
tar cvzf CFS_NIGHTLY.tgz CFS_TESTSUITE_NIGHTLY CFSDEPS_NIGHTLY CFS_TRUNK_NIGHTLY

# Start Ubuntu 8.04 VBox in background
VBOX_LOG=$HOME/Documents/dev/ubuntu804.log
echo "Ubuntu 8.04 VBox started on `date`" > $VBOX_LOG
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo >> $VBOX_LOG
VBoxHeadless -startvm "Ubuntu 8.04" -vrdp on >> $VBOX_LOG 2>&1 &
PID=$!
echo "PID of VBox: $PID" >> $VBOX_LOG
echo "-----------------------------------------------------------------------------" >> $VBOX_LOG
echo "Ubuntu 8.04 VBox started in background on `date`" >> $VBOX_LOG
echo >> $VBOX_LOG


# Change into CFS++ source directory and execute CTest for GCC.
DISTRO=$($TESTDIR/CFS_TRUNK_NIGHTLY/share/scripts/distro.sh -u)

PerformTest "rom_gcc_nightly"

CFSBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfsbin"
CFSTOOLBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfstoolbin"
CPLREADER="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cplreader"
if [ -f $CFSBIN ] && [ -f $CFSTOOLBIN ] && [ -f $CPLREADER ]; then
    # Copy binaries to /opt/pckg/cfs_nightly
    rm -rf $DESTDIR/trunk_gcc
    mkdir $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_gcc
    cp -a $TESTDIR/CFS_TRUNK_NIGHTLY/share/matlab $DESTDIR/trunk_gcc/share
fi

# Update CFS docu for rom webpage
cd "$TESTDIR/CFS_BUILD_NIGHTLY"      && \
make doc-devel > /dev/null 2>&1      && \
make doc-user > /dev/null 2>&1       && \
make doc-user-html > /dev/null 2>&1;

cp -a "$TESTDIR/CFS_TRUNK_NIGHTLY/share/doc/developer/html" "/srv/www/htdocs/cfsDocu/devel"
cp -a "$TESTDIR/CFS_TRUNK_NIGHTLY/share/doc/user/xmlFile/html" "/srv/www/htdocs/cfsDocu/user"

cd /srv/www/htdocs/cfsDocu/xmlManual/                  && \
svn up > /dev/null 2>&1                                && \
pdflatex -halt-on-error cfsManual.tex > /dev/null 2>&1 && \
pdflatex -halt-on-error cfsManual.tex > /dev/null 2>&1 && \
./buildhtml > /dev/null 2>&1


# Change into CFS++ source directory and execute CTest for ICC.
source /opt/intel/Compiler/11.1/069/bin/iccvars.sh intel64
source /opt/intel/Compiler/11.1/069/bin/ifortvars.sh intel64

PerformTest "rom_icc_nightly"

if [ -f $CFSBIN ] && [ -f $CFSTOOLBIN ] && [ -f $CPLREADER ]; then
    # Copy binaries to /opt/pckg/cfs_nightly
    ICC_ROM_DIR=$DESTDIR/trunk_icc_rom
    rm -rf $ICC_ROM_DIR
    mkdir $ICC_ROM_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $ICC_ROM_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $ICC_ROM_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $ICC_ROM_DIR
    cp -a $TESTDIR/CFS_TRUNK_NIGHTLY/share/matlab $ICC_ROM_DIR/share

    rm -rf /opt/pckg/CFSDEPSCACHE/precompiled/*
    cp -a $TESTDIR/CFSDEPSCACHE/precompiled/* /opt/pckg/CFSDEPSCACHE/precompiled
fi


echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
