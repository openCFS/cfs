#!/bin/sh

MODE="-h"

. `dirname $0`/distro.sh -h

# This script is not primarily meant to be used on its own. It will get merged
# into    one   script    bootstrap_devel_machine.txt   with    distro.sh   by
# share/doc/devel/makedoc.sh before  the actual build of  the developer manual
# starts.   bootstrap_devel_machine.txt  will then  be  embedded  into to  the
# resulting PDF  file. If you  wish to make  changes to the  bootstrap process
# make your changes to share/scripts/boostrap_devel_machine.sh and commit your
# changes back to the repository.

ARCH=$(echo $ARCH | sed 'y/'$LOWER'/'$UPPER'/')
DIST=$(echo $DIST | sed 'y/'$LOWER'/'$UPPER'/')
REV=$(echo $REV | sed 'y/'$LOWER'/'$UPPER'/')

ENV_CFS=/etc/env_cfs.sh

uid=$(/usr/bin/id -u) && [ "$uid" = "0" ] ||
{ echo "$0 must be run as root!"; exit 1; }

echo "DIST: $DIST"
echo "REV: $REV"
echo "ARCH: $ARCH"

ExitFail() {
    echo "A problem has occured. Please try to fix the reason and try again."
    exit 1;
}

SetupDebian() {
    # Setup Debian, Ubuntu, Linux Mint, etc...
    PCKGS="subversion gcc g++ gfortran automake autoconf cmake \
        perl-base graphviz texlive-latex-base tex4ht xsltproc \
        python-pygments doxygen tcl-dev python-dev git-svn \
        cmake-curses-gui cmake-qt-gui gmsh default-jre openjdk-6-jdk \
        patch diff diffutils zip libxt-dev libxp6 tk-dev xsltproc \
        libgl1-mesa-dev libglu1-mesa-dev libxmuu-dev libncurses5-dev \
        util-linux gcc-multilib libxmu-dev"

    for pckg in $PCKGS; do
        apt-get install -y --force-yes -f $pckg
    done
    rm -f /usr/lib/libXext.so /usr/lib/libXmu.so
    ln -s /usr/lib/libXext.so.6.4.0 /usr/lib/libXext.so || ExitFail
    ln -s /usr/lib/libXmu.so.6.2.0 /usr/lib/libXmu.so || ExitFail
}

SetupSuse() {
    MAJOR=$(echo $REV | cut -d'.' -f1)
    MINOR=$(echo $REV | cut -d'.' -f2)

    if [ "$ARCH" = "I386" ]; then SUSEARCH=x86; fi
    if [ "$ARCH" = "X86_64" ]; then SUSEARCH=x64; fi

    PCKGS="subversion gcc gcc-c++ gcc-fortran automake autoconf \
        cmake perl-base perl-Switch graphviz texlive-latex texlive-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        java-1_6_0-openjdk-devel cmake-gui xorg-x11-libXt-devel \
        diffutils patch zip xorg-x11-libXp tk-devel Mesa-devel \
        ncurses-devel perl util-linux glibc-devel-32bit libXmu-devel\
        gcc-32bit gcc-c++-32bit gcc-fortran-32bit glu-devel"

    
    # We need to add the SDK DVDs as repos in case of SLE
    if [ "$DIST" = "SLE" ]; then 
        SLETYPE=$(cat /etc/SuSE-release | grep SUSE | cut -d' ' -f4)

        case "$SLETYPE" in
            Server) SLE="SLES";;
            Desktop) SLE="SLED";;
        esac

        unset MINOR

        if [ "$MAJOR" = "10" ]; then
            zypper sa http://demeter.uni-regensburg.de/${SLE}10SP3-$SUSEARCH/DVD1 DVD1-Regensburg
            # zypper sa http://demeter.uni-regensburg.de/${SLE}10SP3-$SUSEARCH/DVD2 DVD2-Regensburg
            zypper sa http://demeter.uni-regensburg.de/SLES10SP3-SDK-$SUSEARCH/DVD1 SDK-DVD1
            zypper sa http://demeter.uni-regensburg.de/SLES10SP3-SDK-$SUSEARCH/DVD2 SDK-DVD2

	    PCKGS="$PCKGS java-1_5_0-ibm java-1_5_0-ibm-devel"
        fi

        if [ "$MAJOR" = "11" ]; then
            zypper ar http://demeter.uni-regensburg.de/${SLE}11SP1-$SUSEARCH/DVD1 DVD1-Regensburg
            zypper ar http://demeter.uni-regensburg.de/${SLE}11SP1-$SUSEARCH/DVD2 DVD2-Regensburg
            zypper ar http://demeter.uni-regensburg.de/SLE11SP1-SDK-$SUSEARCH/DVD1 SDK-DVD1
            zypper ar http://demeter.uni-regensburg.de/SLE11SP1-SDK-$SUSEARCH/DVD2 SDK-DVD2

	    PCKGS="$PCKGS java-1_6_0-ibm java-1_6_0-ibm-devel"
        fi
    else
        PCKGS="$PCKGS java-1_7_0-openjdk-devel"
        PCKGS="$PCKGS java-1_6_0-openjdk-devel"
    fi

    for pckg in $PCKGS
    do
        zypper --non-interactive install $pckg
    done

    if [ "$MAJOR" = "11" ] && [ $MINOR -ge 4 ]; then
        REPO="http://download.opensuse.org/repositories/science/openSUSE_$REV/science.repo"
        for repo in $REPO
        do
          zypper addrepo $repo
          if [ $? -eq 0 ]; then
            break;
          fi
        done
        if [ $? -ne 0 ]; then
          ExitFail
        fi
        zypper install gmsh || ExitFail
    fi

    if [ "$MAJOR" = "12" ]; then
        REPO="http://download.opensuse.org/repositories/science/openSUSE_$REV/science.repo"
        for repo in $REPO
        do
          zypper addrepo $repo
          if [ $? -eq 0 ]; then
            break;
          fi
        done
        if [ $? -ne 0 ]; then
          ExitFail
        fi
        zypper install gmsh || ExitFail
    fi

}

SetupFedora() {
    # /usr/bin/audit from package audit conflicts with binutils.
    yum remove audit

    yum install subversion gcc gcc-c++ gcc-gfortran automake autoconf cmake \
        perl graphviz texlive-latex tetex-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        cmake-gui java-1.6.0-openjdk-devel java-1.7.0-openjdk-devel tk-devel \
        patch diffutils zip libXt-devel libXp ncurses-devel \
        mesa-libGL-devel mesa-libGLU-devel libXmu-devel mesa-libglapi || ExitFail
}

SetupRHEL() {
    # Setup Red Hat Enterprise Linux, CentOS, Oracle, Scientific Linux, etc.
    RHEL_REL=$(echo $REV | cut -d'.' -f1)

    # At the moment we support RHEL 5-7
    SUPPORTED=0
    # For RHEL 7 we only tested with the public Beta version for which no
    # additional repositories were available yet (2014-02-27).
    ADD_ADDITIONAL_REPOS=1
    case "${RHEL_REL}" in
	5) SUPPORTED=1 ;;
	6) SUPPORTED=1 ;;
	7) SUPPORTED=1; ADD_ADDITIONAL_REPOS=0; ;;
	*)
            echo "RHEL release ${RHEL_REL} is NOT supported!"
	    exit 1
            ;;
    esac
    if [ "$SUPPORTED" = 1 ]; then
        echo "Fine! RHEL release ${RHEL_REL} is supported!";
    fi

    if [ "$ADD_ADDITIONAL_REPOS" = 1 ]; then
	cd /etc/yum.repos.d && \
	    rm -f graphviz-rhel.repo || ExitFail 
	wget http://www.graphviz.org/graphviz-rhel.repo || ExitFail
	
	YC=atrpms.repo
	echo "[atrpms]" > $YC && \
	    echo "name=Redhat Enterprise Linux RHEL\$releasever - \$basearch - ATrpms" >> $YC && \
	    echo "baseurl=http://dl.atrpms.net/el${RHEL_REL}-\$basearch/atrpms/stable/" >> $YC && \
	    echo "gpgkey=http://ATrpms.net/RPM-GPG-KEY.atrpms" >> $YC && \
	    echo "gpgcheck=1" >> $YC || ExitFail
	rpm --import http://ATrpms.net/RPM-GPG-KEY.atrpms
	
	YC=epel.repo
	EPEL_MIRROR="http://ftp.uni-bayreuth.de/linux/fedora-epel"
	echo "[epel]" > $YC && \
	    echo "name=EPEL RHEL\$releasever - \$basearch" >> $YC && \
	    echo "baseurl=${EPEL_MIRROR}/${RHEL_REL}/\$basearch" >> $YC || ExitFail
	rm -f RPM-GPG-KEY-EPEL-${RHEL_REL} || ExitFail
	wget ${EPEL_MIRROR}/RPM-GPG-KEY-EPEL-${RHEL_REL} || ExitFail
	rpm --import RPM-GPG-KEY-EPEL-${RHEL_REL}
	
	ARCH=$(uname -m | sed 's/i[0-9]86/i386/') || ExitFail
	BASE=http://pkgs.repoforge.org/rpmforge-release
	case "${RHEL_REL}" in
	    5) RPM=rpmforge-release-0.5.3-1.el5.rf.$ARCH.rpm
		;;
	    6) RPM=rpmforge-release-0.5.3-1.el6.rf.$ARCH.rpm
		;;
	    *)
		echo "RHEL release ${RHEL_REL} is not supported!"
		;;
	esac
	wget $BASE/$RPM || ExitFail
	rpm -Uhv --force $RPM || ExitFail
    fi

    # http://public-yum.oracle.com/
    # http://linux.oracle.com/switch/
    #
    # https://linux.oracle.com/switch/centos/
    # curl -O https://linux.oracle.com/switch/centos2ol.sh
    # sh centos2ol.sh
    if [ "$DIST" = "ORACLE" ]; then
      case "${RHEL_REL}" in
	5) wget http://public-yum.oracle.com/public-yum-el${RHEL_REL}.repo
	    ;;
	6) wget http://public-yum.oracle.com/public-yum-ol${RHEL_REL}.repo
	    ;;
	*)
            echo "RHEL release ${RHEL_REL} is not supported!"
            ;;
      esac
    fi

    yum makecache || ExitFail

    cd /opt && \
    rm -f org.tmatesoft.svn_1.3.5.standalone.zip || ExitFail
    wget http://www.svnkit.com/org.tmatesoft.svn_1.3.5.standalone.zip && \
    unzip org.tmatesoft.svn_1.3.5.standalone.zip || ExitFail

    if [ "$DIST" = "CENTOS" ]; then
	ENABLE_REPO="--enablerepo=centosplus"
    fi

    yum $ENABLE_REPO install fuse-sshfs subversion gcc gcc-c++ \
                perl graphviz.$(uname -m) tetex-latex tetex-tex4ht \
                automake autoconf cmake gcc-gfortran ncurses-devel \
                java-1.6.0-openjdk-devel tk-devel python-pygments doxygen \
                tcl-devel python-devel git-svn patch diffutils zip \
                libXt-devel libXp mesa-libGLU-devel libXmu-devel make \
                glibc-devel.i386 glibc-devel.i686 util-linux-ng util-linux || ExitFail
           
    if [ "$ARCH" = "X86_64" ]; then
	LIB="lib64"
    else
	LIB="lib"
    fi 
    rm -f /usr/$LIB/libXext.so
    ln -s /usr/$LIB/libXext.so.6.4.0 /usr/$LIB/libXext.so || ExitFail

    printf "JAVA_HOME=/usr\n" >> $ENV_CFS
    printf "PATH=/opt/svnkit-1.3.5.7406:\$PATH\n" >> $ENV_CFS
    printf "export JAVA_HOME PATH\n" >> $ENV_CFS

}

SetupMacOS() {
    ISOK=1

    # Install Xcode from the MacOS X installation DVD (optional packages -> Xcode)
    if [ ! -f /Developer/Applications/Xcode.app/Contents/MacOS/Xcode ]; then
	echo "Xcode is  not installed. Please install it from your  MacOS X DVD or"
	echo "download it from Apple."
	ISOK=0
    fi

    if [ ! -f /usr/bin/gcc ]; then
	echo "UNIX  Dev  Support  is  not  installed.  Please  install  it  while"
	echo "installing Xcode."
	ISOK=0
    fi

    # Check for installed MacPorts
    if [ ! -f /opt/local/bin/port ]; then
	echo "MacPorts is not  installed. Please download it from www.macports.org"
	echo "and install it."
	ISOK=0
    fi

    # Check for gfortran.
    if [ ! -f /usr/bin/gfortran-4.2 ]; then
	echo "Could not find /usr/bin/gfortran-4.2!"
	echo "No  proper gfortran is installed.  You can download  and install "
	echo "gfortran from http://r.research.att.com/tools/"
        echo
	echo "Different binaries are also available from:"
	echo "http://www.macresearch.org/files/gfortran/gfortran-4.3-Nov.mpkg.zip"
	echo "http://www.macresearch.org/gfortran-leopard"
	echo "http://hpc.sourceforge.net."
	ISOK=0
    fi

    # Find CMake
    TMPFILE=$(mktemp -t bootstrap) || exit 1
    find /Applications -name 'CMake*.app' | while read dir; do 
	CMAKE_VERSION=$("$dir/Contents/bin/cmake" --version | cut -d' ' -f3 | sed 's/\(2.8\)\(.*\)/\1/')
	CMAKE_MAJOR_VERSION=$(echo $CMAKE_VERSION | cut -d'.' -f1)
	CMAKE_MINOR_VERSION=$(echo $CMAKE_VERSION | cut -d'.' -f2)
	if [ $CMAKE_MAJOR_VERSION -ge 2 -a $CMAKE_MINOR_VERSION -ge 8 ]; then
	    echo "CMAKEDIR=\"$dir/Contents/bin\"" > $TMPFILE;
	fi
    done

    . $TMPFILE
    rm $TMPFILE
    if [ "$CMAKEDIR" = "" ]; then
	echo "CMake 2.8 not found! Please  go to www.cmake.org and download the latest"
	echo "CMake package for Mac.";
	ISOK=0
    fi

    if [ ! $ISOK = 1 ]; then
	echo "Some dependencies are missing. Please fix this first!"
	exit 1;
    fi

    # Get latest updates and package lists
    /opt/local/bin/port selfupdate || ExitFail

    # Install required packages
    /opt/local/bin/port install doxygen graphviz teTeX texlive_texmf-minimal wget  || ExitFail
    /opt/local/bin/port install git-core +svn || ExitFail

    # Make sure CMake 2.8 is on PATH
    PATH=$CMAKEDIR:$PATH
    export PATH

    printf "PATH=\"$CMAKEDIR\":\$PATH\n" >> $ENV_CFS
    printf "export PATH\n" >> $ENV_CFS
}

SetupNetBSD() {
    BSDARCH=$(echo $ARCH | tr [A-Z] [a-z])
    cd /tmp
    ftp ftp://ftp.NetBSD.org/pub/pkgsrc/current/pkgsrc.tar.gz
    ftp ftp://ftp7.de.netbsd.org/pub/ftp.netbsd.org/pub/NetBSD/NetBSD-$REV/${BSDARCH}/binary/sets/comp.tgz
    ftp ftp://ftp7.de.netbsd.org/pub/ftp.netbsd.org/pub/NetBSD/NetBSD-$REV/${BSDARCH}/binary/sets/xcomp.tgz
    tar -xzf pkgsrc.tar.gz -C /usr
    tar -xzf comp.tgz -C /
    tar -xzf xcomp.tgz -C /
    rm pkgsrc.tar.gz comp.tgz xcomp.tgz

    PKG_PATH="http://ftp.NetBSD.org/pub/pkgsrc/packages/NetBSD/${BSDARCH}/$REV/All"
    export PKG_PATH

    pkg_add cmake gcc44 perl openjdk7 python27 doxygen subversion \
	    binutils automake autoconf gmake scmgit valgrind patch \
	    libxml2 libxslt unzip wget

    echo "PATH=/usr/pkg/gcc44/bin:/usr/pkg/${BSDARCH}--netbsdelf/bin:$PATH" >> $ENV_CFS
}

SetupCMake() {
    CMAKE_VERSION=$(cmake --version | cut -d' ' -f3 | cut -d'-' -f1)
    CMAKE_MAJOR_VERSION=$(echo $CMAKE_VERSION | cut -d'.' -f1)
    CMAKE_MINOR_VERSION=$(echo $CMAKE_VERSION | cut -d'.' -f2)
    CMAKE_PATCH_LEVEL=$(echo $CMAKE_VERSION | cut -d'.' -f3)

    if [ $CMAKE_MAJOR_VERSION -ge 2 ] && [ $CMAKE_MINOR_VERSION -ge 8 ] && [ $CMAKE_PATCH_LEVEL -ge 10 ]; then
        return 1
    fi

    CMAKE_MAJOR_VERSION=2
    CMAKE_MINOR_VERSION=8
    CMAKE_PATCH_LEVEL=12

    PCKG_BASE_NAME="cmake-$CMAKE_MAJOR_VERSION.$CMAKE_MINOR_VERSION.$CMAKE_PATCH_LEVEL";
    MYTMPDIR="$TMPDIR/$(basename $0).$$"
    echo "$MYTMPDIR"

    (umask 077 && mkdir "$MYTMPDIR") || exit 1

    cd "$MYTMPDIR"

    # Define list of mirrors
    mirrors="http://www.cmake.org/files/v$CMAKE_MAJOR_VERSION.$CMAKE_MINOR_VERSION/$PCKG_BASE_NAME.tar.gz
             ftp://lse17.e-technik.uni-erlangen.de:40065/cfsdeps/sources/cmake/$PCKG_BASE_NAME.tar.gz"

    MD5SUM="17c6513483d23590cbce6957ec6d1e66"

    # Download source
    for mirror in $mirrors; do
        if [ -f $PCKG_BASE_NAME.tar.gz ]; then
            rm -f $PCKG_BASE_NAME.tar.gz
        fi
        wget --timeout=30 $mirror
        if [ $? -eq 0 ]; then break; fi
    done

    MD5SUM_ACTUAL=""
    # Check MD5 sum on CMake >= 2.8
    MD5SUM_ACTUAL="$MD5SUM_ACTUAL $(cmake -E md5sum $PCKG_BASE_NAME.tar.gz 2> /dev/null | cut -d' ' -f1)"
    # Check MD5 sum on Mac OS X
    MD5SUM_ACTUAL="$MD5SUM_ACTUAL $(md5 2> /dev/null | cut -d'=' -f2 | cut -d' ' -f2)"
    # Check MD5 sum on Linux
    MD5SUM_ACTUAL="$MD5SUM_ACTUAL $(md5sum $PCKG_BASE_NAME.tar.gz 2> /dev/null | sed -e 's/ .*//')"

    isokay=0
    for sum in $MD5SUM_ACTUAL; do
        if [ "$sum" = "$MD5SUM" ]; then isokay=1; break; fi
    done

    if [ $isokay -lt 1 ]; then
        echo "MD5 sums for $PCKG_BASE_NAME do not match!";
        exit 1
    fi

    tar xzf $PCKG_BASE_NAME.tar.gz

    cd $PCKG_BASE_NAME

    CC="gcc"; export CC
    CXX="g++"; export CXX

    sh ./configure --prefix=/opt/$PCKG_BASE_NAME --no-qt-gui

    make

    make install

    cd ..
    rm -rf "$MYTMPDIR"

    printf "\n# Setting PATH to CMake.\n" >> $ENV_CFS
    printf "PATH=/opt/$PCKG_BASE_NAME/bin:\$PATH\n" >> $ENV_CFS
}

HINTSTR="Please add the command\n\n"
HINTSTR="${HINTSTR}. $ENV_CFS\n\n"
HINTSTR="${HINTSTR}to one of the following files\n\n"
HINTSTR="${HINTSTR}\$HOME/.bashrc (Linux)\n"
HINTSTR="${HINTSTR}\$HOME/.profile (Mac)\n"
HINTSTR="${HINTSTR}/etc/profile.local\n"
HINTSTR="${HINTSTR}/etc/profile\n\n"
HINTSTR="${HINTSTR}and open a new shell for all environment settings to become active.\n\n"

echo "# Configuration for CFS++ shell environment." > $ENV_CFS
printf "$HINTSTR" | sed 's/^/# /' >> $ENV_CFS

case "$DIST" in
     MACOSX) SetupMacOS;;
     OPENSUSE) SetupSuse ;;
     SLE) SetupSuse ;;
     DEBIAN) SetupDebian ;;
     UBUNTU) SetupDebian ;;
     LINUXMINT) SetupDebian ;;
     LMDE) SetupDebian ;;
     FEDORA) SetupFedora ;;
     RHEL) SetupRHEL ;;
     CENTOS) SetupRHEL ;;
     SCIENTIFIC) SetupRHEL ;;
     ORACLE) SetupRHEL ;;
     NETBSD) SetupNetBSD ;;
     *)
        echo "Your distribution $DIST is currently not supported by this script."
        echo "You are encouraged to contribute a new boostrap routine by taking"
        echo "the Setup* functions in $0 as a reference for implementing a function"
        echo "for $DIST. When finished please place your modified $0 into"
        echo "CFS_SOURCE_DIR/share/scripts and commit it to the Subversion repo"
        echo "or send it to one of the CFS++ developers."
        exit 1 ;;
esac

SetupCMake

echo 
echo "The configuration to set up an environment for CFS++ development has been"
echo "written to $ENV_CFS."
printf "$HINTSTR"
