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
        cmake-curses-gui cmake-qt-gui gmsh openjdk-6-jdk \
        patch diff diffutils zip libxt-dev libxp6 tk-dev \
        libgl1-mesa-dev libglu1-mesa-dev libxmuu-dev libncurses5-dev"

    for pckg in $PCKGS; do
        apt-get install -y -f $pckg
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
    
    # We need to add the SDK DVDs as repos in case of SLE
    if [ "$DIST" = "SLE" ]; then 
        SLETYPE=$(cat /etc/SuSE-release | grep SUSE | cut -d' ' -f4)

        if [ ! "$MAJOR" = "11" ]; then
            echo "Only SLE 11 is supported at the moment!"
        fi 
       
        unset MINOR
 
        case "$SLETYPE" in
            Server) SLE="SLES";;
            Desktop) SLE="SLED";;
        esac
        
        zypper ar http://demeter.uni-regensburg.de/${SLE}11SP1-$SUSEARCH/DVD1 DVD1-Regensburg
        zypper ar http://demeter.uni-regensburg.de/${SLE}11SP1-$SUSEARCH/DVD2 DVD2-Regensburg
        zypper ar http://demeter.uni-regensburg.de/SLE11SP1-SDK-$SUSEARCH/DVD1 SDK-DVD1
        zypper ar http://demeter.uni-regensburg.de/SLE11SP1-SDK-$SUSEARCH/DVD2 SDK-DVD2
    fi

 

    zypper install subversion gcc gcc-c++ gcc-fortran automake autoconf \
        cmake perl-base perl-Switch graphviz texlive-latex texlive-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        java-1_6_0-openjdk-devel cmake-gui xorg-x11-libXt-devel \
        diffutils patch zip xorg-x11-libXp tk-devel Mesa-devel \
        ncurses-devel


    if [ "$MAJOR" = "11" ] && [ $MINOR -ge 3 ]; then
	REPO="http://download.opensuse.org/repositories/home:/tsokar/openSUSE_$REV \
              http://download.opensuse.org/repositories/home:/scorot/openSUSE_$REV"
        for repo in $REPO
        do
          zypper addrepo --check $repo Gmsh
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
    yum install subversion gcc gcc-c++ gcc-gfortran automake autoconf cmake \
        perl graphviz texlive-latex tetex-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        cmake-gui java-1.6.0-openjdk-devel tk-devel \
        patch diffutils zip libXt-devel libXp ncurses-devel || ExitFail
}

SetupRHEL() {
    # Setup Red Hat Enterprise Linux, CentOS, Oracle, Scientific Linux, etc.
    RHEL_REL=$(echo $REV | cut -d'.' -f1)

    case "${RHEL_REL}" in
	5) echo "Fine! RHEL release ${RHEL_REL} is supported!";;
	6) echo "Fine! RHEL release ${RHEL_REL} is supported!";;
	*)
            echo "RHEL release ${RHEL_REL} is NOT supported!"
            ;;
    esac

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
    echo "[epel]" > $YC && \
    echo "name=EPEL RHEL\$releasever - \$basearch" >> $YC && \
    echo "baseurl=http://ftp.ucr.ac.cr/epel/${RHEL_REL}/\$basearch" >> $YC || ExitFail
    rm -f RPM-GPG-KEY-EPEL-${RHEL_REL} || ExitFail
    wget http://ftp.ucr.ac.cr/epel/RPM-GPG-KEY-EPEL-${RHEL_REL} || ExitFail
    rpm --import RPM-GPG-KEY-EPEL-${RHEL_REL}

    ARCH=$(uname -m | sed 's/i[0-9]86/i386/') || ExitFail
    BASE=http://apt.sw.be/redhat/el${RHEL_REL}/en
    case "${RHEL_REL}" in
	5) RPM=$ARCH/rpmforge/RPMS/rpmforge-release-0.3.6-1.el5.rf.$ARCH.rpm
	    ;;
	6) RPM=$ARCH/rpmforge/RPMS/rpmforge-release-0.5.2-2.el6.rf.$ARCH.rpm
	    ;;
	*)
            echo "RHEL release ${RHEL_REL} is not supported!"
            ;;
    esac
    rpm -Uhv --force $BASE/$RPM || ExitFail

    yum makecache || ExitFail

    cd /opt && \
    rm -f org.tmatesoft.svn_1.3.5.standalone.zip || ExitFail
    wget http://www.svnkit.com/org.tmatesoft.svn_1.3.5.standalone.zip && \
    unzip org.tmatesoft.svn_1.3.5.standalone.zip || ExitFail

    yum --enablerepo=centosplus install fuse-sshfs subversion gcc gcc-c++ \
                perl graphviz.$(uname -m) tetex-latex tetex-tex4ht \
                automake autoconf cmake gcc-gfortran ncurses-devel \
                java-1.6.0-openjdk-devel tk-devel python-pygments doxygen \
                tcl-devel python-devel git-svn patch diffutils zip \
                libXt-devel libXp mesa-libGLU-devel libXmu-devel || ExitFail
           
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

    if [ $CMAKE_MAJOR_VERSION -ge 2 ] && [ $CMAKE_MINOR_VERSION -ge 8 ]; then
        return 1
    fi

    PCKG_BASE_NAME="cmake-2.8.6";
    MYTMPDIR="$TMPDIR/$(basename $0).$$"
    echo "$MYTMPDIR"

    (umask 077 && mkdir "$MYTMPDIR") || exit 1

    cd "$MYTMPDIR"

    # Define list of mirrors
    mirrors="http://distfiles.macports.org/cmake/$PCKG_BASE_NAME.tar.gz
             http://gentoo.tups.lv/source/distfiles/$PCKG_BASE_NAME.tar.gz
             http://130.230.54.100/gentoo/distfiles/$PCKG_BASE_NAME.tar.gz
             http://www.cmake.org/files/v2.8/$PCKG_BASE_NAME.tar.gz"

    MD5SUM="2147da452fd9212bb9b4542a9eee9d5b"

    # Download source
    for mirror in $mirrors; do
        if [ -f $PCKG_BASE_NAME.tar.gz ]; then
            rm -f $PCKG_BASE_NAME.tar.gz
        fi
        wget $mirror
        if [ $? -eq 0 ]; then break; fi
    done

    MD5SUM_ACTUAL=""
    # Check MD5 sum on CMake >= 2.8
    MD5SUM_ACTUAL="$MD5SUM_ACTUAL $(cmake -E md5sum cmake-2.8.6.tar.gz 2> /dev/null | cut -d' ' -f1)"
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
     FEDORA) SetupFedora ;;
     RHEL) SetupRHEL ;;
     CENTOS) SetupRHEL ;;
     SCIENTIFIC) SetupRHEL ;;
     NETBSD) SetupNetBSD ;;
     *)
        echo "Your distribution $DIST is currently not supported by this script."
        echo "You are encouraged to contribute a new boostrap routine by taking"
        echo "the Setup* functions in $0 as a reference for implementing a function"
        echo "for $DIST. When finished please place your modified $0 into"
	echo "CFS_SOURCE_DIR/share/scripts and commit it to the Subversion repo"
	echo "or send it to one of the CFS++ developers."
        exit 1
        ;;
esac

SetupCMake

echo 
echo "The configuration to set up an environment for CFS++ development has been"
echo "written to $ENV_CFS."
printf "$HINTSTR"
