#!/bin/sh

# Set some common variables for CFS start script and soapserve.sh
#
# strieben, june 2007

# Set language
LANG=en_US.UTF-8
export LANG

# Get name of CFS/bin directory
CFS_BIN_DIR=`dirname $0`

CFS_ROOT_DIR=$CFS_BIN_DIR/..

# Check if we are in BINARY tree or in DIST tree.
if [ ! -d $CFS_ROOT_DIR/source ]; then
    # We are in DIST tree
#    echo "Script has been started from DIST tree."
    BINARY_TREE=0
else
    # We are in BINARY tree
#    echo "Script has been started from BINARY tree."
    BINARY_TREE=1
fi

DISTRO_SH=$CFS_ROOT_DIR/share/scripts/distro.sh
if [ ! -x "$DISTRO_SH" ]; then
  echo "Could not execute $DISTRO_SH. Exiting..."
  exit 1
fi

# Get operating system
OS=$($DISTRO_SH -h | cut -d' ' -f1)

# Get architecture and distribution
ARCH_STR=$($DISTRO_SH -a)

DISTRO=`echo $ARCH_STR | cut -d' ' -f 1`
DISTRO_VER=`echo $ARCH_STR | cut -d' ' -f 2`
ARCH=`echo $ARCH_STR | cut -d' ' -f 3`
ARCH_STR="${DISTRO}_${DISTRO_VER}_${ARCH}"

# Set lib path according to architecture
case "$ARCH" in
    I386)
	LIB="lib" ;;
    X86_64)
	LIB="lib64" ;;
    IA64)
	LIB="lib" ;;
    *)
	LIB="lib" ;;
esac

case "$OS" in
    Linux)
	OS_SUPPORTED=1 ;;
    Mac)
	OS_SUPPORTED=1 ;;
    *)
	OS_SUPPORTED=0 ;;
esac

# At the moment we only support Linux
if [ ! "$OS_SUPPORTED" = "1" ]
then
    echo "Your operating system $OS is not supported."
    echo "CFS supports only Linux and Mac at the moment."
    exit 1
fi

# Set library path for current architecture
if [ "$BINARY_TREE" = "0" ]; then
    if [ -z "$LD_LIBRARY_PATH" ]
    then 
        LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR";
    else
        LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR:$LD_LIBRARY_PATH";
    fi
else
    if [ -z "$LD_LIBRARY_PATH" ]
    then 
        LD_LIBRARY_PATH="$CFSDEPS_ROOT/$LIB:$CFS_ROOT_DIR/$LIB/$ARCH_STR";
    else
        LD_LIBRARY_PATH="$CFSDEPS_ROOT/$LIB:$CFS_ROOT_DIR/$LIB/$ARCH_STR:$LD_LIBRARY_PATH";
    fi
fi
export LD_LIBRARY_PATH;

if [ "$DISTRO" = "MACOSX" ]; then
    LDD="otool -L"
else
    LDD="ldd"
fi
