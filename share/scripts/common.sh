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
eval $($DISTRO_SH -s)
OS_DUMMY=$(echo $OS | sed 's/ //g')
OS=$OS_DUMMY

# Get architecture and distribution
if [ ! "$DIST_FAMILY" = "" ]; then
  DISTRO=$DIST_FAMILY
  DISTRO_VER=$MAJOR_REV
else
  DISTRO=$DIST
  DISTRO_VER=$REV
fi

POSSIBLE_ARCHS="I386"
if [ "$ARCH" = "X86_64" ]; then
    POSSIBLE_ARCHS="$ARCH $POSSIBLE_ARCHS"
fi

ARCH_STR="${DISTRO}_${DISTRO_VER}_${ARCH}"

case "$OS" in
    LINUX)
	OS_SUPPORTED=1
	if [ ! -d "$CFS_ROOT_DIR/bin/$ARCH_STR" ]; then
	    ARCH_STR_FOUND=0
	    
	    for a in $POSSIBLE_ARCHS; do
		LINUX_BIN="$CFS_ROOT_DIR/bin/${DISTRO}_${DISTRO_VER}_${a}/cfsbin"
		if [ "$CFS_SCRIPT_DEBUG" = "1" ]; then
		    echo "Trying $LINUX_BIN..."
		fi
                if [ -f "$LINUX_BIN" ]; then
		    if [ "$CFS_SCRIPT_DEBUG" = "1" ]; then
			echo "Found $LINUX_BIN!"
			fi
		    ARCH=$a
		    ARCH_STR="${DISTRO}_${DISTRO_VER}_${ARCH}"
		    ARCH_STR_FOUND=1
		    break;
                fi
	    done
	fi
	;;
    MACOSX)
	LD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
	OS_SUPPORTED=1 
	if [ ! -d "$CFS_ROOT_DIR/bin/$ARCH_STR" ]; then
	    ARCH_STR_FOUND=0
	    VERSIONS="10.9 10.8 10.7 10.6 10.5"
	    
	    for v in $VERSIONS; do
		for a in $POSSIBLE_ARCHS; do
		    OSX_BIN="$CFS_ROOT_DIR/bin/${DISTRO}_${v}_${a}/cfsbin"
		    if [ "$CFS_SCRIPT_DEBUG" = "1" ]; then
			echo "Trying $OSX_BIN..."
		    fi
                    if [ -f "$OSX_BIN" ]; then
			if [ "$CFS_SCRIPT_DEBUG" = "1" ]; then
			    echo "Found $OSX_BIN!"
			fi
			DISTRO_VER=$v
			ARCH=$a
			ARCH_STR="${DISTRO}_${DISTRO_VER}_${ARCH}"
			ARCH_STR_FOUND=1
			break;
                    fi
		done
		if [ $ARCH_STR_FOUND = 1 ]; then
		    break
		fi
	    done
	fi
	;;
    *)
	OS_SUPPORTED=0 ;;
esac

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
        LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR";
    else
        LD_LIBRARY_PATH="$CFS_ROOT_DIR/$LIB/$ARCH_STR:$LD_LIBRARY_PATH";
    fi
fi

case "$OS" in
    LINUX)
	export LD_LIBRARY_PATH ;;
    MACOSX)
	DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
	export DYLD_LIBRARY_PATH ;;
    *)
	;;
esac


if [ "$DISTRO" = "MACOSX" ]; then
    LDD="otool -L"
else
    LDD="ldd"
fi
