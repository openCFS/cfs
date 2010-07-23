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

function GetWorkingCopyRev {
    WCDIR=$1
    
    # Get previous revision for CFSDEPS
    WC_REV=$(svn st --xml --verbose $WCDIR | xsltproc $TESTDIR/CFSDEPS_NIGHTLY/utils/xslt/cfsdepsrev.xslt -)
    if [ ! "$?" = 0 ]; then WC_REV="-1"; fi

    change-svn-wc-format.py $WCDIR "1.5" --force # --verbose
}

function PerformTest {
    TESTNAME=$1;

    cd $TESTDIR/CFS_TRUNK_NIGHTLY

    # Let's build ParaView every Wednesday
    if [ "$DAYOFWEEK" = "3" -a "$TESTNAME" = "sedici_gcc_ifort_nightly" ]; then
        ADDITIONAL_SED_ARG="s/CPLREADER:BOOL\=ON/CPLREADER:BOOL\=ON\n   BUILD_PARAVIEW:BOOL=ON/g"
        cat ctest_scripts/ctest_$TESTNAME.cmake | sed $ADDITIONAL_SED_ARG > ctest_scripts/ctest_script.cmake
    else
        cat ctest_scripts/ctest_$TESTNAME.cmake > ctest_scripts/ctest_script.cmake
    fi

    rm -rf $TESTDIR/CFS_BUILD_NIGHTLY

    ctest -V -S ctest_scripts/ctest_script.cmake

    cd $TESTDIR
    tar cvzf cfs_build_$TESTNAME.tgz CFS_BUILD_NIGHTLY
    
    rm -rf CFS_BUILD_NIGHTLY
}

function SubmitErrorToCDash {
    BUILDSCRIPT=$(basename "$0");
    SITE=$HOSTNAME;
    BUILDNAME="ATTENTION: $BUILDSCRIPT";
    NOTENAME=$1;
    NOTE="$2";
    ERROR="$3";
    LOGLINE="$4";
    NOTETEMPFILE="$(mktemp).xml"
    BUILDFILE="$NOTETEMPFILE"
    BUILDSTAMP="$(date '+%Y%m%d-%H%M')-Nightly"
    
    echo $SITE
    echo $BUILDNAME
    echo $NOTENAME
    echo $NOTE
    echo $NOTETEMPFILE
    
    echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > $BUILDFILE
    echo "<Site BuildName=\"$BUILDNAME\"" >> $BUILDFILE
    echo "      BuildStamp=\"$BUILDSTAMP\"" >> $BUILDFILE
    echo "      Name=\"$SITE\"" >> $BUILDFILE
    echo "      Generator=\"$BUILDSCRIPT\"    OSName=\"Linux\"" >> $BUILDFILE
    echo "      Hostname=\"$SITE\">" >> $BUILDFILE
#    echo "        OSRelease=\"2.6.31.12-0.1-default\"" >> $BUILDFILE
#    echo "        OSVersion=\"#1 SMP 2010-01-27 08:20:11 +0100\"" >> $BUILDFILE
#    echo "        OSPlatform=\"x86_64\"" >> $BUILDFILE
#    echo "        Is64Bits=\"1\"" >> $BUILDFILE
#    echo "        VendorString=\"AuthenticAMD\"" >> $BUILDFILE
#    echo "        VendorID=\"Advanced Micro Devices\"" >> $BUILDFILE
#    echo "        FamilyID=\"16\"" >> $BUILDFILE
#    echo "        ModelID=\"4\"" >> $BUILDFILE
#    echo "        ProcessorCacheSize=\"512\"" >> $BUILDFILE
#    echo "        NumberOfLogicalCPU=\"16\"" >> $BUILDFILE
#    echo "        NumberOfPhysicalCPU=\"16\"" >> $BUILDFILE
#    echo "        TotalVirtualMemory=\"2047\"" >> $BUILDFILE
#    echo "        TotalPhysicalMemory=\"32239\"" >> $BUILDFILE
#    echo "        LogicalProcessorsPerPhysical=\"1\"" >> $BUILDFILE
#    echo "        ProcessorClockFrequency=\"2511.53\">" >> $BUILDFILE
    echo "  <Build>" >> $BUILDFILE
#    echo "    <StartDateTime>Mar 29 02:23 CEST</StartDateTime>" >> $BUILDFILE
#    echo "    <StartBuildTime>1269822221</StartBuildTime>" >> $BUILDFILE
    echo "    <StartDateTime>" >> $BUILDFILE
    date '+%b %d %R %Z' >> $BUILDFILE
    echo "    </StartDateTime>" >> $BUILDFILE
#    echo "    <StartBuildTime>" >> $BUILDFILE
#    date '+%N' >> $BUILDFILE
#    echo "    </StartBuildTime>" >> $BUILDFILE

    echo "    <BuildCommand>$BUILDSCRIPT</BuildCommand>" >> $BUILDFILE
    echo "    <Error>" >> $BUILDFILE
    echo "      <BuildLogLine>$LOGLINE</BuildLogLine>" >> $BUILDFILE
    echo "      <Text>" >> $BUILDFILE
    echo "      $ERROR" >> $BUILDFILE
    echo "      </Text>" >> $BUILDFILE
    echo "      <Precontext/>" >> $BUILDFILE
    echo "      <Postcontext/>" >> $BUILDFILE
    echo "      <RepeatCount>0</RepeatCount>" >> $BUILDFILE
    echo "    </Error>" >> $BUILDFILE
    echo "    <Log Encoding=\"base64\" Compression=\"/bin/gzip\"/>" >> $BUILDFILE
#    echo "    <EndDateTime>Mar 29 02:52 CEST</EndDateTime>" >> $BUILDFILE
#    echo "    <EndBuildTime>1269823955</EndBuildTime>" >> $BUILDFILE
    echo "    <EndDateTime>" >> $BUILDFILE
    date '+%b %d %R %Z' >> $BUILDFILE
    echo "    </EndDateTime>" >> $BUILDFILE
    echo "    <ElapsedMinutes>0</ElapsedMinutes>" >> $BUILDFILE
    echo "  </Build>" >> $BUILDFILE
    echo "</Site>" >> $BUILDFILE

    curl -T "{$BUILDFILE}" http://lse17.e-technik.uni-erlangen.de:2000/cdash/submit.php?project=CFS

    echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > $NOTETEMPFILE
    echo "<?xml-stylesheet type=\"text/xsl\" href=\"Dart/Source/Server/XSL/Build.xsl <file:///Dart/Source/Server/XSL/Build.xsl> \"?>" >> $NOTETEMPFILE
    echo "<Site BuildName=\"$BUILDNAME\" BuildStamp=\"$BUILDSTAMP\" Name=\"$SITE\" Generator=\"$BUILDSCRIPT\">" >> $NOTETEMPFILE
    echo "  <Notes>" >> $NOTETEMPFILE
    echo "    <Note Name=\"$NOTENAME\">" >> $NOTETEMPFILE
    echo -n "      <DateTime>" >> $NOTETEMPFILE
    date '+%b %d %R %Z' >> $NOTETEMPFILE
    echo "</DateTime>" >> $NOTETEMPFILE
    echo "      <Text>" >> $NOTETEMPFILE

    echo "$NOTE" >> $NOTETEMPFILE
    echo "      </Text>" >> $NOTETEMPFILE
    echo "    </Note>" >> $NOTETEMPFILE

    echo "  </Notes>" >> $NOTETEMPFILE
    echo "</Site>" >> $NOTETEMPFILE

    # Submit machine specific version of Notes.xml to dash board
    curl -T "{$NOTETEMPFILE}" http://lse17.e-technik.uni-erlangen.de:2000/cdash/submit.php?project=CFS
    
    rm -f $NOTETEMPFILE
    
    exit 1;
}

cd $TESTDIR

# Remove previous CFSDEPSCACHE directory and make sure that ACML sources
# get copied to the newly created one.
if [ "$DAYOFWEEK" = "7" ]; then
    mkdir -p update
    cp CFS_TRUNK_NIGHTLY/ctest_scripts/ctest_update* update
    cp CFS_TRUNK_NIGHTLY/ctest_scripts/CTestConfig.cmake update
    tar cvzf $TESTDIR/last_weeks_cfs_trunk_nightly.tgz CFS_TRUNK_NIGHTLY
    rm -rf $TESTDIR/CFS_TRUNK_NIGHTLY
    tar cvzf $TESTDIR/last_weeks_cfsdeps_nightly.tgz CFSDEPS_NIGHTLY
    rm -rf $TESTDIR/CFSDEPS_NIGHTLY
    tar cvzf $TESTDIR/last_weeks_cfs_testsuite_nightly.tgz CFS_TESTSUITE_NIGHTLY
    rm -rf $TESTDIR/CFS_TESTSUITE_NIGHTLY
    rm -rf $TESTDIR/CFSDEPSCACHE
    mkdir -p $TESTDIR/CFSDEPSCACHE/sources/acml
    cp /opt/pckg/CFSDEPSCACHE/sources/acml/*.tgz $TESTDIR/CFSDEPSCACHE/sources/acml

    rm -rf $TESTDIR/cfs_build_*.tgz
fi

# Checkout or update CFSDEPS
cd $TESTDIR/update
GetWorkingCopyRev $TESTDIR/CFSDEPS_NIGHTLY
CFSDEPS_PREV_REV=$WC_REV;
# Due to the fact, that CTest generates non-zero return values
# even if no error occured we do not check for that condition (cf. 
# http://public.kitware.com/Bug/bug_view_page.php?bug_id=8277&history=1).
# Instead we delete all working copies on Sunday (DAYOFWEEK=7) to
# make sure we get a fresh start from time to time.
ctest -V -S ctest_update_cfsdeps_klu.cmake
#|| \
#  SubmitErrorToCDash "ctest_update_cfsdeps_klu.cmake" "Update of CFSDEPS_NIGHTLY failed." "ctest_update_cfsdeps_klu.cmake: Update of CFSDEPS failed."
GetWorkingCopyRev $TESTDIR/CFSDEPS_NIGHTLY
CFSDEPS_CURRENT_REV=$WC_REV;

# Checkout or update CFS++
cd $TESTDIR/update
GetWorkingCopyRev $TESTDIR/CFS_TRUNK_NIGHTLY
CFS_PREV_REV=$WC_REV;
ctest -V -S ctest_update_cfs_klu.cmake
mkdir -p $TESTDIR/update
cp $TESTDIR/CFS_TRUNK_NIGHTLY/ctest_scripts/ctest_update* $TESTDIR/update
cp $TESTDIR/CFS_TRUNK_NIGHTLY/ctest_scripts/CTestConfig.cmake $TESTDIR/update
GetWorkingCopyRev $TESTDIR/CFS_TRUNK_NIGHTLY
CFS_CURRENT_REV=$WC_REV;

# Checkout or update test suite
cd $TESTDIR/update
GetWorkingCopyRev $TESTDIR/CFS_TESTSUITE_NIGHTLY
TESTSUITE_PREV_REV=$WC_REV;
ctest -V -S ctest_update_testsuite_klu.cmake
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

# Copy sources of GotoBLAS to CFSDEPSCACHE/source/gotoblas
mkdir -p CFSDEPSCACHE/sources/gotoblas
cp $DEVDIR/LOCAL_LIBS/GotoBLAS-1.26.tar.gz CFSDEPSCACHE/sources/gotoblas

# Copy PARDISO libs to CFSDEPSCACHE/source/pardiso
#mkdir -p CFSDEPSCACHE/sources/pardiso
#cp $DEVDIR/LOCAL_LIBS/libpardiso_* CFSDEPSCACHE/sources/pardiso
#cp $DEVDIR/LOCAL_LIBS/license.pdf CFSDEPSCACHE/sources/pardiso

# Determine distribution string
DISTRO=$($TESTDIR/CFS_TRUNK_NIGHTLY/share/scripts/distro.sh -u)

# Source Intel scripts to get environment variables
source /opt/intel/Compiler/11.0/081/bin/iccvars.sh intel64
source /opt/intel/Compiler/11.0/081/bin/ifortvars.sh intel64

# Change into CFS++ source directory and execute CTest for GCC.
PerformTest "sedici_gcc_ifort_nightly"

# Let's copy the ParaView binaries to /opt/pckg
if [ "$DAYOFWEEK" = "3" -a -d "$TESTDIR/CFS_BUILD_NIGHTLY/paraview"]; then
    cp -a $TESTDIR/CFS_BUILD_NIGHTLY/paraview/* /opt/pckg/paraview-weekly
fi

# Change into CFS++ source directory and execute CTest for ICC.
PerformTest "sedici_icc_nightly"

# Copy the nightly Intel binaries for sedici to /opt/pckg
CFSBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfsbin"
CFSTOOLBIN="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cfstoolbin"
CPLREADER="$TESTDIR/CFS_BUILD_NIGHTLY/bin/$DISTRO/cplreader"
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


echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
