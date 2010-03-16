#!/bin/sh

# This script is run by the following entry in /etc/crontab on sedici
#
# 0     2 * * *   strieben /home/users/strieben/Documents/dev/nightly_test_sedici.sh > /home/users/strieben/Documents/dev/nightly_test_sedici.log 2>&1

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

# Set ssh environment for Subversion
export SVN_SSH="ssh -i $DEVDIR/svn_id_dsa_cfstestuser_klu -l svn -p 22022 -q"

# Extend system path by /opt/pckg/bin
export PATH=/opt/pckg/bin:$PATH

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
    
    # Get previous revision for working copy 
    PREV_REV=$(svn st --xml --verbose $WCDIR | xsltproc CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    if [ ! "$?" = 0 ]; then PREV_REV="-1"; fi

    change-svn-wc-format.py $WCDIR "1.5" --verbose

    if [ ! -d "$WCDIR/.svn" ]; then
	svn co $REPOS $WCDIR
    else
	svn up $WCDIR
    fi

    # Get current revision for working copy
    CURRENT_REV=$(svn st --xml --verbose $WCDIR | xsltproc CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    if [ ! "$?" = 0 ]; then CURRENT_REV="-1"; fi

    change-svn-wc-format.py $WCDIR "1.5" --verbose

    # Account for 90 seconds reconnect timeout of the lse10 firewall.
    sleep 120
}

function PerformTest {
    TESTNAME=$1;

    if [ "$DAYOFWEEK" = "3" ]; then
# && [ "$TESTNAME" = "sedici_gcc_nightly" ]; then
        ADDITIONAL_SED_ARG="s/CPLREADER:BOOL\=ON/CPLREADER:BOOL\=ON\n   BUILD_PARAVIEW=ON/g"
    else
        ADDITIONAL_SED_ARG="s///g"
    fi

    cd $TESTDIR/CFS_TRUNK_NIGHTLY

    if [ ! "$CFS_CURRENT_REV" = "$CFS_PREV_REV" ]; then
        # If revision for CFS++ has changed start with an empty build directory
	cat ctest_scripts/ctest_$TESTNAME.cmake > ctest_scripts/ctest_script.cmake
    else
        # If revision for CFS++ has not changed reuse old build directory
	cat ctest_scripts/ctest_$TESTNAME.cmake | sed \
	-e 's/EMPTY_BINARY_DIRECTORY TRUE/EMPTY_BINARY_DIRECTORY FALSE/g' \
	-e "$ADDITIONAL_SED_ARG" > ctest_scripts/ctest_script.cmake
	cd $TESTDIR

        # Make sure old build directory gets removed
	rm -rf CFS_BUILD_NIGHTLY
	if [ -r cfs_build_$TESTNAME.tgz ]; then
	    tar xvzf cfs_build_$TESTNAME.tgz

	    # Remove old libs if CFSDEPS have changed
	    if [ ! "$CFSDEPS_CURRENT_REV" = "$CFSDEPS_PREV_REV" ]; then
		rm -rf CFS_BUILD_NIGHTLY/include
		rm -rf CFS_BUILD_NIGHTLY/lib64
	    fi
	fi

	
	cd $TESTDIR/CFS_TRUNK_NIGHTLY
    fi

    ctest -V -S ctest_scripts/ctest_script.cmake

    if [ -d $TESTDIR/CFS_BUILD_NIGHTLY/paraview ]; then
       cp -a $TESTDIR/CFS_BUILD_NIGHTLY/paraview /opt/pckg/paraview-3.6.1
       rm -rf $TESTDIR/CFS_BUILD_NIGHTLY/paraview
    fi

    cd $TESTDIR
    tar cvzf cfs_build_$TESTNAME.tgz CFS_BUILD_NIGHTLY
}


# Remove previous CFSDEPSCACHE directory and make sure that ACML sources
# get copied to the newly created one.
if [ "$DAYOFWEEK" = "7" ]; then
    rm -rf $TESTDIR/CFSDEPSCACHE
    mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml
    cp /opt/pckg/CFSDEPSCACHE/sources/acml/*.tgz $TESTDIR/CFSDEPSCACHE/sources/acml

    rm -rf $TESTDIR/cfs_build_*.tgz
fi

# Either checkout out or update CFS++ trunk working copy. Subversion update
# had some strange problems if performed by CTest.
TRUNK_REPO="svn+ssh://131.188.140.10/software/CFS++/trunk"
UpdateWorkingCopy "$TRUNK_REPO" "CFS_TRUNK_NIGHTLY"
CFS_PREV_REV=$PREV_REV;
CFS_CURRENT_REV=$CURRENT_REV;
echo "CFS_PREV_REV $CFS_PREV_REV"
echo "CFS_CURRENT_REV $CFS_CURRENT_REV"

DEPS_REPO="svn+ssh://131.188.140.10/software/cfsdeps/trunk"
UpdateWorkingCopy "$DEPS_REPO" "CFSDEPS_NIGHTLY"
CFSDEPS_PREV_REV=$PREV_REV;
CFSDEPS_CURRENT_REV=$CURRENT_REV;
echo "CFSDEPS_PREV_REV $CFSDEPS_PREV_REV"
echo "CFSDEPS_CURRENT_REV $CFSDEPS_CURRENT_REV"

# Either checkout out or update CFS++ test suite trunk working copy.
SUITE_REPO="svn+ssh://131.188.140.10/software/CFS++_TEST/trunk"
UpdateWorkingCopy "SUITE_REPO" "CFS_TESTSUITE_NIGHTLY"
TESTSUITE_PREV_REV=$PREV_REV;
TESTSUITE_CURRENT_REV=$CURRENT_REV;
echo "TESTSUITE_PREV_REV $TESTSUITE_PREV_REV"
echo "TESTSUITE_CURRENT_REV $TESTSUITE_CURRENT_REV"

# Remove cache directory if CFSDEPS have changed
if [ ! "$CFSDEPS_CURRENT_REV" = "$CFSDEPS_PREV_REV" ]; then
    rm -rf $TESTDIR/CFSDEPSCACHE
    mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml
    cp /opt/pckg/CFSDEPSCACHE/sources/acml/*.tgz $TESTDIR/CFSDEPSCACHE/sources/acml
fi

# Copy sources of SCPIP to CFSDEPS_NIGHTLY/scpip
# cp $DEVDIR/scpip.tar.bz2 CFSDEPS_NIGHTLY/scpip

# Copy sources of GotoBLAS to CFSDEPSCACHE/source/gotoblas
mkdir -p CFSDEPSCACHE/sources/gotoblas
cp $DEVDIR/LOCAL_LIBS/GotoBLAS-1.26.tar.gz CFSDEPSCACHE/sources/gotoblas

# Copy PARDISO libs to CFSDEPSCACHE/source/pardiso
mkdir -p CFSDEPSCACHE/sources/pardiso
cp $DEVDIR/LOCAL_LIBS/libpardiso_* CFSDEPSCACHE/sources/pardiso

# Copy whole CFSDEPS tree to /opt/pckg/CFSDEPS. This should make sure that
# we allways have a current CFSDEPS present on the machine.
# rm -rf /opt/pckg/CFSDEPS/* /opt/pckg/CFSDEPS/.svn
# cp -a CFSDEPS_NIGHTLY/* CFSDEPS_NIGHTLY/.svn /opt/pckg/CFSDEPS


# Change into CFS++ source directory and execute CTest for GCC.
DISTRO=$($TESTDIR/CFS_TRUNK_NIGHTLY/share/scripts/distro.sh -u)

PerformTest "sedici_gcc_nightly"

CFSBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfsbin"
CFSTOOLBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfstoolbin"
CPLREADER="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cplreader"
#if [ -f $CFSBIN ] && [ -f $CFSTOOLBIN ] && [ -f $CPLREADER ]; then
#    # Copy binaries to /opt/pckg/cfs_nightly
#    rm -rf $DESTDIR/trunk_gcc
#    mkdir $DESTDIR/trunk_gcc
#    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $DESTDIR/trunk_gcc
#    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $DESTDIR/trunk_gcc
#    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $DESTDIR/trunk_gcc
#    cp -a $TESTDIR/CFS_TRUNK_NIGHTLY/share/matlab $DESTDIR/trunk_gcc/share
#    
#fi

# Change into CFS++ source directory and execute CTest for ICC.
source /opt/intel/Compiler/11.0/081/bin/iccvars.sh intel64
source /opt/intel/Compiler/11.0/081/bin/ifortvars.sh intel64

PerformTest "sedici_icc_nightly"

if [ -f $CFSBIN ] && [ -f $CFSTOOLBIN ] && [ -f $CPLREADER ]; then
    # Copy binaries to /opt/pckg/cfs_nightly
    ICC_SEDICI_DIR=$DESTDIR/trunk_icc_sedici
    rm -rf $ICC_SEDICI_DIR 
    mkdir $ICC_SEDICI_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/bin $ICC_SEDICI_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/lib64 $ICC_SEDICI_DIR
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/share $ICC_SEDICI_DIR
    cp -a $TESTDIR/CFS_TRUNK_NIGHTLY/share/matlab $ICC_SEDICI_DIR/share

#    rm -rf /opt/pckg/CFSDEPSCACHE/precompiled/*
#    cp -a $TESTDIR/CFSDEPSCACHE/precompiled/* /opt/pckg/CFSDEPSCACHE/precompiled
fi

# Change into CFS++ source directory and execute CTest for Schenk Pardiso
# using NETLIB, GOTOBLAS and ACML.
PerformTest "sedici_gcc_schenk_nightly"
PerformTest "sedici_gcc_acml_nightly"
PerformTest "sedici_gcc_gotoblas_nightly"
PerformTest "sedici_icc_gfortran_nightly"


# cd $TESTDIR
# rm -rf $TESTDIR/CFSDEPSCACHE/precompiled
# cd $TESTDIR/CFS_TRUNK_NIGHTLY
# ctest -V -S ctest_scripts/ctest_sedici_gcc44_nightly.cmake
# rm -rf $TESTDIR/bin
# mkdir $TESTDIR/bin
# cd $TESTDIR/bin
# ln -s /usr/bin/gcc-4.4 gcc
# ln -s /usr/bin/g++-4.4 g++
# ln -s /usr/bin/gfortran-4.4 gfortran
# export PATH=$TESTDIR/bin:$PATH
# cd $TESTDIR/CFS_TRUNK_NIGHTLY
# ctest -V -S ctest_scripts/ctest_sedici_icc_gcc44_nightly.cmake

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
