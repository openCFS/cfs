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
                # On Mandrake/Mandriva/Fedora there exist also
                # /etc/redhat-release, /etc/mandrake-release,
                # /etc/mandriva-release, /etc/fedora-release files.
                # They all contain the same infos. I.e.
                # Mandriva Linux release 2007.0 (Official) for i586
                # Fedora Core release 6 (Zod)
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
		DIST="Debian"
                BASE_VERSION=`dpkg -p base-files | grep Version | cut -d' ' -f2`
		REV=`echo $BASE_VERSION | sed -e 's/.[^.]*$//'`
		case "$REV" in
		    "1.2") PSEUDONAME="rex";;
		    "1.3") PSEUDONAME="bo";;
		    "2.0") PSEUDONAME="hamm";;
		    "2.1") PSEUDONAME="slink";;
		    "2.2") PSEUDONAME="potato";;
		    "3.0") PSEUDONAME="woody";;
		    "3.1") PSEUDONAME="sid";;
                esac
                PSEUDONAME="$PSEUDONAME `cat /etc/debian_version`"

                SPINOFF=`echo $BASE_VERSION | cut -d'.' -f3 | sed -e 's/[0-9]*//g'`
		case "$SPINOFF" in
		    "ubuntu")
			. /etc/lsb-release;
			DIST=$DISTRIB_ID;
			REV=$DISTRIB_RELEASE;
			PSEUDONAME=$DISTRIB_CODENAME;;
                esac

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

