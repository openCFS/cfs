#!/bin/sh
# Detects which OS and if it is Linux then it will detect which Linux Distribution.
# Based upon: http://linuxmafia.com/faq/Admin/release-files.html

OS=`uname -s`
DIST=$OS
REV=`uname -r`
MACH=`uname -m`
ARCH=`echo $MACH | sed "s/i[0-9]/i3/"`

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
#	exit 1
fi

# Set default sub-architecture
SUBARCH="SUBARCHUNKNOWN"

if [ "${OS}" = "SunOS" ] ; then
	ARCH=`uname -p`
        DIST=`uname -n`
	OSSTR="${OS} ${REV} (${DIST} ${ARCH} `uname -v`)"
elif [ "${OS}" = "AIX" ] ; then
	OSSTR="${OS} `oslevel` (`oslevel -r`)"
elif [ "${OS}" = "Linux" ] ; then

        # If using 32-bit CMake on a 64-bit platform the uname command
        # returns i386 instead of ia64 or x86_64. The arch command does
        # the right thing.

        ARCH_AFFIRM=$(arch 2>&1)
	if [ $? = 0 ]; then
	    if [ "$ARCH_AFFIRM" = "ia64" ] ; then
		ARCH="IA64";
            fi
            if [ "$ARCH_AFFIRM" = "x86_64" ] ; then
		ARCH="X86_64";
            fi
	    ARCH_AFFIRM=`echo $ARCH_AFFIRM | sed "s/i[0-9]/i3/"`
            if [ "$ARCH_AFFIRM" = "i386" ] ; then
		ARCH="I386";
            fi
	fi

	KERNEL=`uname -r`

        # Now let's determine the sub-architecture. This can be EM64T or
	# OPTERON for X86_64 or SGI for IA64
        if [ "$ARCH" = "X86_64" ] ; then
            SUBARCH=$(grep -i amd /proc/cpuinfo)
	    
	    if [ "$SUBARCH" = "" ]; then
	         SUBARCH="EM64T"
            else
	         SUBARCH="OPTERON"
	    fi
        fi
        if [ "$ARCH" = "IA64" ] ; then
	    if [ -f /dev/sgi_fetchop ]; then
	         SUBARCH="SGI"
	    fi
        fi
	

	if [ -f /etc/redhat-release ] ; then
                # On Mandrake/Mandriva/Fedora there exist also
                # /etc/redhat-release, /etc/mandrake-release,
                # /etc/mandriva-release, /etc/fedora-release files.
                # They all con tain the same infos. I.e.
                # Mandriva Linux release 2007.0 (Official) for i586
                # Fedora Core release 6 (Zod)
		# DIST='RedHat'
                #
                # http://fedoraproject.org/wiki/History_of_Red_Hat_Linux
                # http://fedoraproject.org/wiki/Releases/HistoricalSchedules
	        DIST=`cat /etc/redhat-release | cut -d' ' -f1`
		PSEUDONAME=`cat /etc/redhat-release | sed s/.*\(// | sed s/\)//`
		REV=`cat /etc/redhat-release | sed s/.*release\ // | sed s/\ .*//`
	elif [ -f /etc/SuSE-release ] ; then
            SUSEREL="/etc/SuSE-release"
            FIRSTLINE=`head -1 $SUSEREL | sed 'y/'$LOWER'/'$UPPER'/'`
            ENTERPRISE=`echo $FIRSTLINE | cut -f3 -d' '`
            if [ "$ENTERPRISE" = "ENTERPRISE" ]; then
                DIST="SLE"
                REV=`echo $FIRSTLINE | cut -f5 -d' '`
                PSEUDONAME=`head -3 $SUSEREL | tail -1 | sed 's/ = //'`
            else
                DIST=`cat $SUSEREL | tr "\n" ' '| sed s/VERSION.*// | awk '{print $1}'`
                REV=`cat $SUSEREL | tr "\n" ' ' | awk '{print $2}'`
                PSEUDONAME=$(cat /etc/issue | tr ' ' '\n' | grep '\".*' | sed 's/\"//g')
            fi
	elif [ -f /etc/mandrake-release ] ; then
		DIST='Mandrake'
		PSEUDONAME=`cat /etc/mandrake-release | sed s/.*\(// | sed s/\)//`
		REV=`cat /etc/mandrake-release | sed s/.*release\ // | sed s/\ .*//`
	elif [ -f /etc/debian_version -o -f /etc/debian-version ]; then
		DIST="Debian"
                BASE_VERSION=`dpkg -p base-files 2> /dev/null | grep Version`
		if [ ! $? -eq 0 ]; then
                	BASE_VERSION=`apt-cache show base-files 2> /dev/null | grep Version`
		fi
                # echo $BASE_VERSION
		REV=$(echo $BASE_VERSION | cut -d' ' -f2 | grep -o '[0-9]*\.[0-9]*' )
		REV_VERSION=$(echo $BASE_VERSION | cut -d' ' -f2 | grep -o '[0-9]*$')
                if [ "$REV_VERSION" = "" ]; then REV_VERSION=0; fi
                # echo $REV_VERSION
		# echo $REV
		# See http://www.us.debian.org/doc/FAQ/ch-ftparchives.de.html
		case "$REV" in
		    "1.1") PSEUDONAME="buzz";;
		    "1.2") PSEUDONAME="rex";;
		    "1.3") PSEUDONAME="bo";;
		    "2.0") PSEUDONAME="hamm";;
		    "2.1") PSEUDONAME="slink";;
		    "2.2") PSEUDONAME="potato";;
		    "3.0") PSEUDONAME="woody";;
		    "3.1") PSEUDONAME="sarge";;
		    "4.0") PSEUDONAME="etch";;
                    "5.0") PSEUDONAME="lenny";;
                    "6.0") PSEUDONAME="squeeze";;
                    "7.0") PSEUDONAME="wheezy";;
                    "8.0") PSEUDONAME="sid";;
                esac
#                PSEUDONAME="$PSEUDONAME `cat /etc/debian_version`"

		if [ -f /etc/knoppix-version ]; then
			SPINOFF=knoppix;
		else
                	SPINOFF=`echo $BASE_VERSION | cut -d'.' -f3 | sed -e 's/[0-9]*//g'`;
                        if [ "$SPINOFF" = "" ]; then
                            SPINOFF=`echo $BASE_VERSION | sed -e 's/[0-9\.]*//g'`;
                        fi
		fi

		case "$SPINOFF" in
		    "ubuntu")
			. /etc/lsb-release;
			DIST=$DISTRIB_ID;
			REV=$DISTRIB_RELEASE;
			PSEUDONAME=$DISTRIB_CODENAME;
                        # https://wiki.ubuntu.com/DevelopmentCodeNames
			# http://en.wikipedia.org/wiki/History_of_Ubuntu
			case "$DISTRIB_CODENAME" in
			    "warty") PSEUDONAME="Warty Warthog";; # 4.10
			    "hoary") PSEUDONAME="Hoary Hedgehog";; # 5.04
			    "breezy") PSEUDONAME="Breezy Badger";; # 5.10
			    "dapper") PSEUDONAME="Dapper Drake";; # 6.06 LTS
			    "edgy") PSEUDONAME="Edgy Eft";; # 6.10
			    "feisty") PSEUDONAME="Feisty Fawn";; # 7.04
			    "gutsy") PSEUDONAME="Gutsy Gibbon";; # 7.10
			    "hardy") PSEUDONAME="Hardy Heron";; # 8.04 LTS
			    "intrepid") PSEUDONAME="Intrepid Ibex";; # 8.10
			    "jaunty") PSEUDONAME="Jaunty Jackalope";; # 9.04
			    "karmic") PSEUDONAME="Karmic Koala";; # 9.10
                            "lucid") PSEUDONAME="Lucid Lynx";; # 10.04 LTS
                            "maverick") PSEUDONAME="Maverick Meerkat";; # 10.10
                            "natty") PSEUDONAME="Natty Narwhal";; # 11.04
                            "oneiric") PSEUDONAME="Oneiric Ocelot";; # 11.10
                            "precise") PSEUDONAME="Precise Pangolin";; # 12.04
	                 esac;;
		    "knoppix")
			DIST=Knoppix;
			REV=`cat /etc/knoppix-version | cut -d' ' -f1`
		 	PSEUDONAME="Knoppix";;
                esac

        elif [ -f /etc/lsb-release ]; then
            . /etc/lsb-release;
            DIST=$DISTRIB_ID;
            REV=$DISTRIB_RELEASE;
            PSEUDONAME=$DISTRIB_CODENAME;
        fi
	if [ -f /etc/UnitedLinux-release ] ; then
		DIST="${DIST}[`cat /etc/UnitedLinux-release | tr "\n" ' ' | sed s/VERSION.*//`]"
	fi
	
	OSSTR="${OS} ${DIST} ${REV} (${PSEUDONAME} ${KERNEL} ${MACH})"

elif [ ${OS} = "Darwin" ]; then
        MACOSINFO=$(system_profiler SPSoftwareDataType | grep 'System Version')
	if [ $? -eq 0 ]; then
	    MACOSVER=$(echo $MACOSINFO | grep 'System Version' | cut -d':' -f2 | cut -d'X' -f2 | cut -d'(' -f1)
	    OS="Mac OS X"
            DIST="MACOSX"
	    FULL_REV=$MACOSVER
            MAJOR_REV=`echo $FULL_REV | sed s/\.[0-9]$//`

	    case "$MAJOR_REV" in
                "10.0") PSEUDONAME="Cheetah";;
                "10.1") PSEUDONAME="Puma";;
                "10.2") PSEUDONAME="Jaguar";;
                "10.3") PSEUDONAME="Panther";;
		"10.4") PSEUDONAME="Tiger";;
		"10.5") PSEUDONAME="Leopard";;
                "10.6") PSEUDONAME="Snow Leopard";;
                "10.7") PSEUDONAME="Lion";;
	    esac

            REV=$MAJOR_REV
            OSSTR="$OS $DIST $MAJOR_REV ($FULL_REV $PSEUDONAME ${MACH})"
	else
	    DIST="DARWIN"
            OSSTR="$OS $DIST $REV ($PSEUDONAME ${MACH})"
	fi
fi

      while :
      do
          case "$1" in
              -h) echo ${OSSTR} ;;
              -a) echo "${DIST} ${REV} ${ARCH} ${SUBARCH}" | sed 'y/'$LOWER'/'$UPPER'/';;
              -u) echo "${DIST}_${REV}_${ARCH}" | sed 'y/'$LOWER'/'$UPPER'/' ;;
              -c) echo "${DIST};${REV};${ARCH};${SUBARCH}" | sed 'y/'$LOWER'/'$UPPER'/' ;;
              *) break ;;
          esac
          shift
      done

