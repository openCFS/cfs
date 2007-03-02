#!/bin/sh
# Detects which OS and if it is Linux then it will detect which Linux Distribution.


OS=`uname -s`
REV=`uname -r`
MACH=`uname -m`
ARCH=`echo $MACH | sed "s/i./i3/"`

LOWER='abcdefghijklmnopqrstuvwxyz'
UPPER='ABCDEFGHIJKLMNOPQRSTUVWXYZ'

GetVersionFromFile()
{
	VERSION=`cat $1 | tr "\n" ' ' | sed s/.*VERSION.*=\ // `
}

Help()
{
	echo "Usage: $0 [-h | -a | -u]"
	echo
	echo "-h      print distribution info human readable"
	echo "-a      print distribution info awk parsable"
	echo "-u      print distribution info unique"
}

if [ -z "$1" ]; then
	Help
	exit 1
fi

if [ "${OS}" = "SunOS" ] ; then
	OS=Solaris
	ARCH=`uname -p`	
	OSSTR="${OS} ${REV}(${ARCH} `uname -v`)"
elif [ "${OS}" = "AIX" ] ; then
	OSSTR="${OS} `oslevel` (`oslevel -r`)"
elif [ "${OS}" = "Linux" ] ; then
	KERNEL=`uname -r`
	if [ -f /etc/redhat-release ] ; then
		DIST='RedHat'
		PSEUDONAME=`cat /etc/redhat-release | sed s/.*\(// | sed s/\)//`
		REV=`cat /etc/redhat-release | sed s/.*release\ // | sed s/\ .*//`
	elif [ -f /etc/SuSE-release ] ; then
		DIST=`cat /etc/SuSE-release | tr "\n" ' '| sed s/VERSION.*// | awk '{print $1}'`
		REV=`cat /etc/SuSE-release | tr "\n" ' ' | sed s/.*=\ // | awk '{print $1}'`
	elif [ -f /etc/mandrake-release ] ; then
		DIST='Mandrake'
		PSEUDONAME=`cat /etc/mandrake-release | sed s/.*\(// | sed s/\)//`
		REV=`cat /etc/mandrake-release | sed s/.*release\ // | sed s/\ .*//`
	elif [ -f /etc/debian_version ] ; then
		DIST="Debian `cat /etc/debian_version`"
		REV=""
		. /etc/lsb-release
		DIST=$DISTRIB_ID
		REV=$DISTRIB_RELEASE
		PSEUDONAME=$DISTRIB_CODENAME
	fi
	if [ -f /etc/UnitedLinux-release ] ; then
		DIST="${DIST}[`cat /etc/UnitedLinux-release | tr "\n" ' ' | sed s/VERSION.*//`]"
	fi
	
	OSSTR="${OS} ${DIST} ${REV} (${PSEUDONAME} ${KERNEL} ${MACH})"

fi

      while :
      do
          case "$1" in
              -h) echo ${OSSTR} ;;
              -a) echo "${DIST} ${REV} ${ARCH}" | sed 'y/'$LOWER'/'$UPPER'/';;
              -u) echo "${DIST}_${REV}_${ARCH}" | sed 'y/'$LOWER'/'$UPPER'/' ;;
              -c) echo "${DIST};${REV};${ARCH}" | sed 'y/'$LOWER'/'$UPPER'/' ;;
              *) break ;;
          esac
          shift
      done

