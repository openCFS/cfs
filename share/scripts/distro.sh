#!/bin/sh
# -*- mode: sh; coding: utf-8; indent-tabs-mode: nil; -*-
# vim:fenc=utf-8:ft=sh:et:sw=2:ts=2:sts=2
# kate: space-indent on; indent-width 2; encoding utf-8;
# kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#
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
    echo "-c      print distribution info for cmake"
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
    
    LSB_REL=$(which lsb_release 2> /dev/null)
    if [ "$?" = "0" -a ! "$LSB_REL" = "" -a "-f $LSB_REL" ]; then	
        ID=$($LSB_REL -i -s)
        DESC=$($LSB_REL -d -s)
        if [ "$ID" = "LinuxMint" ]; then
            LMDE=$(echo $DESC | cut -d' ' -f1)
            if [ "$LMDE" = "LMDE" ]; then
                DESC=$LMDE
            else
                DESC=$ID
            fi
        fi
        # Check for RHEL
        IS_RHEL=$(echo $ID | sed -n '/RedHatEnterprise/{s|\(RedHatEnterprise\)\(.*\)|\1|};p')
        DIST=$(echo $DESC | sed 's/"//g' | cut -d' ' -f1)
        REV=$($LSB_REL -r -s)
        PSEUDONAME=$($LSB_REL -c -s)
                
        case "$DIST" in
            "openSUSE") if grep -q Tumbleweed /etc/os-release; then
                # for the roling release tumbleweed the revision is something like 20150727 and no pesudoname is set 
                REV="Tumbleweed"
              fi ;;
            "SUSE") DIST="SLE" ;;
            "Debian") REV=$(echo $REV | sed 's/\.[0-9]*$//') ;;
            "Enterprise") DIST="ORACLE" ;;
            "Red") if [ "$IS_RHEL" = "RedHatEnterprise" ]; then                       
                       DIST="RHEL"
                       REV=$(echo $REV | sed 's/\.[0-9]*$//') 
                   fi ;;
        esac
    elif [ -f /etc/lsb-release ]; then
        . /etc/lsb-release;
        DIST=$DISTRIB_ID;
        REV=$DISTRIB_RELEASE;
        PSEUDONAME=$DISTRIB_CODENAME;
    elif [ -f /etc/redhat-release ] ; then
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
        DIST=$(rpm -qf /etc/redhat-release | cut -d'-' -f1)
        case "$DIST" in
            "enterprise")
               # https://blogs.oracle.com/VDIpier/entry/how_to_check_if_the
               DIST="ORACLE"
               PSEUDONAME=$(cat /etc/enterprise-release | cut -d'(' -f2 | cut -d ')' -f1)
               ;;
            *) DIST=$(cat /etc/redhat-release | cut -d' ' -f1)
               if [ "$DIST" = "Red" ]; then
                   DIST="RHEL";
               fi
               PSEUDONAME=`cat /etc/redhat-release | sed s/.*\(// | sed s/\)//`
               ;;
        esac
        
        REV=`cat /etc/redhat-release | sed s/.*release\ // | sed s/\ .*//`

    elif [ -f /etc/SuSE-release ] ; then
    	# as OpenSuse has lsb_release we don't come here in this elif case any more
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
    elif [ -f /etc/arch-release ] ; then
        DIST='Arch'
        REV='rolling'
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
            "8.0") PSEUDONAME="jessie";;
            "9.0") PSEUDONAME="sid";;
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
                    "warty")    PSEUDONAME="Warty Warthog";;    # 4.10
                    "hoary")    PSEUDONAME="Hoary Hedgehog";;   # 5.04
                    "breezy")   PSEUDONAME="Breezy Badger";;    # 5.10
                    "dapper")   PSEUDONAME="Dapper Drake";;     # 6.06 LTS
                    "edgy")     PSEUDONAME="Edgy Eft";;         # 6.10
                    "feisty")   PSEUDONAME="Feisty Fawn";;      # 7.04
                    "gutsy")    PSEUDONAME="Gutsy Gibbon";;     # 7.10
                    "hardy")    PSEUDONAME="Hardy Heron";;      # 8.04 LTS
                    "intrepid") PSEUDONAME="Intrepid Ibex";;    # 8.10
                    "jaunty")   PSEUDONAME="Jaunty Jackalope";; # 9.04
                    "karmic")   PSEUDONAME="Karmic Koala";;     # 9.10
                    "lucid")    PSEUDONAME="Lucid Lynx";;       # 10.04 LTS
                    "maverick") PSEUDONAME="Maverick Meerkat";; # 10.10
                    "natty")    PSEUDONAME="Natty Narwhal";;    # 11.04
                    "oneiric")  PSEUDONAME="Oneiric Ocelot";;   # 11.10
                    "precise")  PSEUDONAME="Precise Pangolin";; # 12.04 LTS
                    "quantal")  PSEUDONAME="Quantal Quetzal";;  # 12.10
                    "raring")   PSEUDONAME="Raring Ringtail";;  # 13.04
                    "saucy")    PSEUDONAME="Saucy Salamander";; # 13.10
                    "trusty")   PSEUDONAME="Trusty Tahr";;      # 14.04 LTS
                esac;;
            "knoppix")
                DIST=Knoppix;
                REV=`cat /etc/knoppix-version | cut -d' ' -f1`
                PSEUDONAME="Knoppix";;
        esac
    fi
       
    if [ -f /etc/UnitedLinux-release ] ; then
        DIST="${DIST}[`cat /etc/UnitedLinux-release | tr "\n" ' ' | sed s/VERSION.*//`]"
    fi

    OSSTR="${OS} ${DIST} ${REV} (${PSEUDONAME} ${KERNEL} ${MACH})"
elif [ ${OS} = "Darwin" ]; then
    MACOSINFO=$(system_profiler SPSoftwareDataType | grep 'System Version')
    if [ $? -eq 0 ]; then
        # up to 10.11 "OS X 10.11.6 (15G1004)" from 10.12 "macOS 10.12 (16A323"
        MACOSVER=$(echo $MACOSINFO | grep 'System Version' | cut -d':' -f2 | cut -d'X' -f2 | cut -d'S' -f2 | cut -d'(' -f1)
        OS="Mac OS X"
        DIST="MACOSX"
        DIST_FAMILY="MACOSX"
        FULL_REV=$(echo $MACOSVER | sed -e 's/^[ \t]*//')
        MAJOR_REV2=`echo $FULL_REV | sed -e 's/\(\.[0-9]\)\(\.[0-9]\)$/\1/' -e 's/Server //'`
        MAJOR_REV="$(echo $MAJOR_REV2 | cut -d. -f1).$(echo $MAJOR_REV2 | cut -d. -f2)"
        REV=`echo $FULL_REV | sed -e 's/Server //' -e 's/ //g'`

        # http://www.imore.com/os-x-version-code-names
        case "$MAJOR_REV" in
            "10.0") PSEUDONAME="Cheetah";;
            "10.1") PSEUDONAME="Puma";;
            "10.2") PSEUDONAME="Jaguar";;
            "10.3") PSEUDONAME="Panther (Pinot)";;
            "10.4") PSEUDONAME="Tiger (Merlot, Intel: Chardonay)";;
            "10.5") PSEUDONAME="Leopard (Chablis)";;
            "10.6") PSEUDONAME="Snow Leopard"; ARCH="X86_64";;
            "10.7") PSEUDONAME="Lion (Barolo)"; ARCH="X86_64";;
            "10.8") PSEUDONAME="Mountain Lion (Zinfandel)"; ARCH="X86_64";;
            "10.9") PSEUDONAME="Mavericks (Cabernet)"; ARCH="X86_64";;
            "10.10") PSEUDONAME="Yosemite (Sirah)"; ARCH="X86_64";;
            "10.11") PSEUDONAME="El Capitan"; ARCH="X86_64";;
            "10.12") PSEUDONAME="Sierra"; ARCH="X86_64";;
        esac

        OSSTR="$OS $DIST $MAJOR_REV ($FULL_REV $PSEUDONAME ${MACH})"
    else
        DIST="DARWIN"
        OSSTR="$OS $DIST $REV ($PSEUDONAME ${MACH})"
    fi
fi

case "$(echo $DIST | sed 'y/'$LOWER'/'$UPPER'/')" in
    "SCIENTIFIC") DIST_FAMILY="RHEL"; MAJOR_REV=$(echo ${REV} | sed -e 's/\.[0-9.]*$//') ;;
    "CENTOS") DIST_FAMILY="RHEL"; MAJOR_REV=$(echo ${REV} | sed -e 's/\.[0-9.]*$//');;
    "ORACLE") DIST_FAMILY="RHEL"; MAJOR_REV=$(echo ${REV} | sed -e 's/\.[0-9.]*$//') ;;
    "ROCKY") DIST_FAMILY="RHEL"; MAJOR_REV=$(echo ${REV} | sed -e 's/\.[0-9.]*$//') ;;
    "RHEL") DIST_FAMILY="RHEL"; MAJOR_REV=$(echo ${REV} | sed -e 's/\.[0-9.]*$//') ;;
    "SLE") DIST_FAMILY="SLE"; MAJOR_REV=${REV} ;;
    *) break ;;
esac

while :
do
    case "$1" in
        -h) # Human readable
echo
echo ${OSSTR}
echo $(cat << DELIM
#
OS: ${OS}#
DIST: ${DIST}#
DIST_FAMILY: ${DIST_FAMILY}#
REV: ${REV}#
MAJOR_REV: ${MAJOR_REV}#
ARCH: ${ARCH}#
SUBARCH: ${SUBARCH}#
DELIM
) | sed 'y/'$LOWER'/'$UPPER'/' | tr '#' '\n' ;;
        -a) echo "${DIST} ${REV} ${ARCH} ${SUBARCH}" | sed 'y/'$LOWER'/'$UPPER'/';;
        -u) echo "${DIST}_${REV}_${ARCH}" | sed 'y/'$LOWER'/'$UPPER'/' ;;
        -c) # CMake syntax
echo $(cat << DELIM
SET(OS "${OS}")
SET(DIST "${DIST}")
SET(DIST_FAMILY "${DIST_FAMILY}")
SET(REV "${REV}")
SET(MAJOR_REV "${MAJOR_REV}")
SET(ARCH "${ARCH}")
SET(SUBARCH "${SUBARCH}")
DELIM
) | sed 'y/'$LOWER'/'$UPPER'/;s/) /)#/g' | tr '#' '\n' ;;
        -s) # Shell syntax
echo $(cat << DELIM
OS="${OS}"
DIST="${DIST}"
DIST_FAMILY="${DIST_FAMILY}"
REV="${REV}"
MAJOR_REV="${MAJOR_REV}"
ARCH="${ARCH}"
SUBARCH="${SUBARCH}"
DELIM
) | sed 'y/'$LOWER'/'$UPPER'/;s/" /"#/g' | tr '#' '\n' ;;
        -p) # Perl syntax
echo $(cat << DELIM
\$OS='${OS}';
\$DIST='${DIST}';
\$DIST_FAMILY='${DIST_FAMILY}';
\$REV='${REV}';
\$MAJOR_REV='${MAJOR_REV}';
\$ARCH='${ARCH}';
\$SUBARCH='${SUBARCH}';
DELIM
) | sed 'y/'$LOWER'/'$UPPER'/;s/" /"#/g' | tr '#' '\n' ;;
        *) break ;;
    esac
    shift
done

