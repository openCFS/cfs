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

    change-svn-wc-format.py $WCDIR "1.5" --force # --verbose
}

function PerformTest {
    TESTNAME=$1;

    cd $TESTDIR/CFS_TRUNK_NIGHTLY

    cat ctest_scripts/ctest_$TESTNAME.cmake > ctest_scripts/ctest_script.cmake

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
if [ "$DAYOFWEEK" = "2" ]; then
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
    rm -rf $TESTDIR/cfs_build_*.tgz
fi

# Copy current version of update scripts from server to update directory
mkdir -p $TESTDIR/update
cd $TESTDIR/update
FILES="ctest_update.cmake ctest_update_cfs_klu.cmake ctest_update_cfsdeps_klu.cmake ctest_update_testsuite_klu.cmake CTestConfig.cmake"
for FILE in $FILES; do
  wget --http-user=testuser-klu --http-password=$CFS_TESTUSER_PW \
       --no-check-certificate \
       https://lse17.e-technik.uni-erlangen.de:2001/svn/CFS++/trunk/ctest_scripts/$FILE
done


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
# cd "$TESTDIR/CFS_BUILD_NIGHTLY"      && \
# make doc-devel > /dev/null 2>&1      && \
# make doc-user > /dev/null 2>&1       && \
# make doc-user-html > /dev/null 2>&1;

# cp -a "$TESTDIR/CFS_TRUNK_NIGHTLY/share/doc/developer/html" "/srv/www/htdocs/cfsDocu/devel"
# cp -a "$TESTDIR/CFS_TRUNK_NIGHTLY/share/doc/user/xmlFile/html" "/srv/www/htdocs/cfsDocu/user"

# cd /srv/www/htdocs/cfsDocu/xmlManual/                  && \
# svn up > /dev/null 2>&1                                && \
# pdflatex -halt-on-error cfsManual.tex > /dev/null 2>&1 && \
# pdflatex -halt-on-error cfsManual.tex > /dev/null 2>&1 && \
# ./buildhtml > /dev/null 2>&1


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
