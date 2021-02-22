#!/bin/sh

# Set some common variables for CFS start scripts
#
# strieben, june 2007
# ftoth, 2021

# Set language
LANG=en_US.UTF-8
export LANG

# Get name of CFS/bin directory
CFS_BIN_DIR=`dirname $0`

# the root dir of the installation
CFS_ROOT_DIR=$CFS_BIN_DIR/..

# set legacy arch-string and lib-dir 
# or in case we want to support different ones in the future (again)
ARCH_STR=LINUX_X86_64
LIB=lib64

# set LD_LIBRARY_PATH such that shared libs of the installation are found
if [ -z "$LD_LIBRARY_PATH" ]
then 
    LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR";
else
    LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR:$LD_LIBRARY_PATH";
fi
export LD_LIBRARY_PATH;
