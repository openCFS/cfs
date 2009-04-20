#!/bin/sh

#mount | grep sedici/strieben &> /dev/null

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

if [ ! -f "/Volumes/strieben/Documents/dev/CFS_TRUNK_NIGHTLY.tgz" ]; then
    echo "Share sedici/strieben is not mounted!"
    exit 1
fi

. /Users/simon/.bashrc

echo "Updating CFS++..."
cd /Users/simon/Documents/dev && \
if [ ! -d "CFS_TRUNK_NIGHTLY" ]; then
    tar xvzf /Volumes/strieben/Documents/dev/CFS_TRUNK_NIGHTLY.tgz
fi
cd CFS_TRUNK_NIGHTLY && \
svn_lse up
if [ "$?" = "0" ]; then
    cd ..
    tar cvzf CFS_TRUNK_NIGHTLY.tgz CFS_TRUNK_NIGHTLY && \
    cp CFS_TRUNK_NIGHTLY.tgz /Volumes/strieben/Documents/dev/
    if [ "$?" = "0" ]; then
	echo "CFS=\"OK\"" > /Volumes/strieben/Documents/dev/CFS_NIGHTLY_BUILD.sh
	rm -rf CFS_TRUNK_NIGHTLY
    fi
fi
echo

echo "Updating CFSDEPS..."
cd /Users/simon/Documents/dev
if [ ! -d "CFSDEPS_NIGHTLY" ]; then
    tar xvzf /Volumes/strieben/Documents/dev/CFSDEPS_NIGHTLY.tgz
fi
cd CFSDEPS_NIGHTLY && \
svn_lse up
if [ "$?" = "0" ]; then
    cd ..
    tar cvzf CFSDEPS_NIGHTLY.tgz CFSDEPS_NIGHTLY && \
    cp CFSDEPS_NIGHTLY.tgz /Volumes/strieben/Documents/dev/
    if [ "$?" = "0" ]; then
	echo "CFSDEPS=\"OK\"" >> /Volumes/strieben/Documents/dev/CFS_NIGHTLY_BUILD.sh
	rm -rf CFSDEPS_NIGHTLY
    fi
fi
echo

echo "Updating CFS++ Testsuite..."
cd /Users/simon/Documents/dev
if [ ! -d "CFS_TESTSUITE_NIGHTLY" ]; then
    tar xvzf /Volumes/strieben/Documents/dev/CFS_TESTSUITE_NIGHTLY.tgz
fi
cd CFS_TESTSUITE_NIGHTLY && \
svn_lse up
if [ "$?" = "0" ]; then
    cd ..
    tar cvzf CFS_TESTSUITE_NIGHTLY.tgz CFS_TESTSUITE_NIGHTLY && \
    cp CFS_TESTSUITE_NIGHTLY.tgz /Volumes/strieben/Documents/dev/
    if [ "$?" = "0" ]; then
	echo "TESTSUITE=\"OK\"" >> /Volumes/strieben/Documents/dev/CFS_NIGHTLY_BUILD.sh
	rm -rf CFS_TESTSUITE_NIGHTLY
    fi
fi
echo

echo "-----------------------------------------------------------------------------"
echo "$0 Done!"
echo
echo
