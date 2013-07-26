#!/bin/sh

echo "`basename $0` started on `date`"
echo "-----------------------------------------------------------------------------"
echo

export HOME=/Users/simon
DEVDIR="$HOME/Documents/dev"
export PARDISO_LIC_PATH="$HOME/Documents/dev/pardiso/keys/mac"

. $HOME/.bashrc

cd $DEVDIR/CFS_TRUNK
ctest -V -S ctest_scripts/ctest_mac_gcc_experimental.cmake

source /opt/intel/cc/10.1.015/bin/iccvars.sh
source /opt/intel/fc/10.1.021/bin/ifortvars.sh
ctest -V -S ctest_scripts/ctest_mac_icc_experimental.cmake

echo "-----------------------------------------------------------------------------"
echo "`basename $0` finished on `date`"
echo
