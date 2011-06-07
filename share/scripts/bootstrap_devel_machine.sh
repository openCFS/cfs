
ARCH=$(echo $ARCH | sed 'y/'$LOWER'/'$UPPER'/')
DIST=$(echo $DIST | sed 'y/'$LOWER'/'$UPPER'/')
REV=$(echo $REV | sed 'y/'$LOWER'/'$UPPER'/')

uid=$(/usr/bin/id -u) && [ "$uid" = "0" ] ||
{ echo "$0 must be run as root!"; exit 1; }

echo "DIST: $DIST"
echo "REV: $REV"
echo "ARCH: $ARCH"

function ExitFail {
    echo "A problem has occured. Please try to fix the reason and try again."
    exit 1;
}

function SetupDebian {
    apt-get install subversion gcc g++ gfortran automake autoconf cmake \
        perl-base graphviz texlive-latex-base tex4ht \
        python-pygments doxygen tcl-dev python-dev git-svn \
        cmake-curses-gui cmake-qt-gui gmsh openjdk-6-jdk \
        patch diffutils zip libxt-dev libxp6 tk-dev \
        libgl1-mesa-dev libglu1-mesa-dev libxmuu-dev || ExitFail

    ln -s /usr/lib/libXext.so.6.4.0 /usr/lib/libXext.so || ExitFail
    ln -s /usr/lib/libXmu.so.6.2.0 /usr/lib/libXmu.so || ExitFail
}

function SetupOpenSuse {
    zypper install subversion gcc gcc-c++ gcc-fortran automake autoconf \
        cmake perl-base graphviz texlive-latex texlive-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        java-1_6_0-openjdk-devel cmake-gui xorg-x11-libXt-devel \
        diffutils patch zip xorg-x11-libXp tk-devel Mesa-devel || ExitFail

    MAJOR=$(echo $REV | cut -d'.' -f1)
    MINOR=$(echo $REV | cut -d'.' -f2)

    if [ "$MAJOR" = "11" ] && [ $MINOR -ge 3 ]; then
	REPO="http://download.opensuse.org/repositories/home:/tsokar/openSUSE_$REV"
	zypper addrepo $REPO Gmsh || ExitFail
	zypper install gmsh || ExitFail
    fi
}

function SetupFedora {
    yum install subversion gcc gcc-c++ gcc-gfortran automake autoconf cmake \
        perl graphviz texlive-latex tetex-tex4ht \
        python-pygments doxygen tcl-devel python-devel git-svn \
        cmake-gui java-1.6.0-openjdk-devel tk-devel \
        patch diffutils zip libXt-devel libXp || ExitFail
}

function SetupCentOS {
    cd /etc/yum.repos.d && \
    wget http://www.graphviz.org/graphviz-rhel.repo || ExitFail

    YC=/etc/yum.conf

    echo "[atrpms]" >> $YC && \
    echo "name=Redhat Enterprise Linux RHEL$releasever - $basearch - ATrpms" >> $YC && \
    echo "baseurl=http://dl.atrpms.net/el5-$basearch/atrpms/stable/" >> $YC && \
    echo "gpgkey=http://ATrpms.net/RPM-GPG-KEY.atrpms" >> $YC && \
    echo "gpgcheck=1" >> $YC || ExitFail

    echo "[epel]" >> $YC && \
    echo "name=EPEL RHEL$releasever - $basearch" >> $YC && \
    echo "baseurl=http://ftp.ucr.ac.cr/epel/5/$basearch" >> $YC || ExitFail

    ARCH=$(uname -m | sed 's/i[0-9]86/i386/') || ExitFail
    BASE=http://apt.sw.be/redhat/el5/en
    RPM=$ARCH/rpmforge/RPMS/rpmforge-release-0.3.6-1.el5.rf.$ARCH.rpm
    rpm -Uhv $BASE/$RPM || ExitFail

    yum makecache || ExitFail

    cd /opt && \
    wget http://www.svnkit.com/org.tmatesoft.svn_1.3.5.standalone.zip && \
    unzip org.tmatesoft.svn_1.3.5.standalone.zip || ExitFail

    yum install subversion gcc gcc-c++ gcc-gfortran automake autoconf cmake \
                perl graphviz tetex-latex tetex-tex4ht \
                python-pygments doxygen tcl-devel python-devel git-svn \
                cmake-gui java-1.6.0-openjdk-devel tk-devel\
                patch diffutils zip libXt-devel libXp mesa-libGLU-devel libXmu-devel || ExitFail
            
    ln -s /usr/lib64/libXext.so.6.4.0 /usr/lib64/libXext.so || ExitFail

    echo "export JAVA_HOME=/usr" >> $HOME/.bashrc || ExitFail
    echo "export PATH=/opt/svnkit-1.3.5.7406:$PATH" >> $HOME/.bashrc || ExitFail

    echo "Please open a new shell for all environment settings to become active."
}

function SetupMacOS {
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
    if [ ! -f /usr/local/bin/gfortran ]; then
	echo "No  proper gfortran is installed.  You can download  and install the"
	echo "one from"
	echo "http://www.macresearch.org/files/gfortran/gfortran-4.3-Nov.mpkg.zip"
	echo "Also be sure to visit"
	echo "http://www.macresearch.org/gfortran-leopard and"
	echo "http://hpc.sourceforge.net."
	ISOK=0
    fi

    # Find CMake
    for dir in "`find /Applications -name 'CMake*'`"; do
	if [ -f "$dir/Contents/bin/cmake" ]; then
	    echo "Using cmake in $dir";
	    CMAKEDIR="$dir/Contents/bin";
            break;
	fi
    done

    if [ "$CMAKEDIR" == "" ]; then
	echo "CMake not found! Please  go to www.cmake.org and download the latest"
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

    # Adjust environment of user with a small perl script
    cat > /tmp/adjustenv.pl << EOF
open(IN, "<\$ENV{'HOME'}/.profile");
@profile=<IN>;
close(IN);

open(OUT, ">\$ENV{'HOME'}/.profile2");
foreach \$line (@profile) {
    if(\$line =~ /export PATH/) {
	@tokens = split(/=/, \$line);
	@tokens = split(/:/, \$tokens[1]);

	my \$last = undef; 
	@pathcontents = grep { (\$last ne \$_) && (\$last = \$_) } sort(@tokens, @pathcontents); 
	
	print OUT "# \$line"
    } else {
	print OUT \$line;
    }
}

\$newpath="export PATH=";
foreach \$token (@pathcontents) {
    if(\$token !~ /\\\$PATH/) {
	chomp(\$token);
#	print "\$token\n";
	\$newpath="\$newpath\$token:";
    }

    if(\$token =~ /\/opt\/local\/bin/) {
	\$haveoptlocalbin = 1;
    }

    if(\$token =~ /\/opt\/local\/sbin/) {
	\$haveoptlocalsbin = 1;
    }
}

if( ! \$haveoptlocalbin ) {
    \$newpath="\${newpath}/opt/local/bin:";
}

if( ! \$haveoptlocalsbin ) {
    \$newpath="\${newpath}/opt/local/sbin:";
}

\$newpath="\${newpath}'$CMAKEDIR':";
\$newpath="\${newpath}\\\$PATH";

print OUT "\$newpath\n";
close(OUT);

print "A new configuration for your environment has been written to \$ENV{'HOME'}/.profile2.\n";
print "Either move it to \$ENV{'HOME'}/.profile or add the following line to your \$ENV{'HOME'}/.profile:\n\n";
print "\$newpath\n";

EOF

    perl /tmp/adjustenv.pl || ExitFail

    rm -f /tmp/adjustenv.pl || ExitFail

    echo "Please open a new shell for all environment settings to become active."
}

case "$DIST" in
     MACOSX) SetupMacOS;;
     OPENSUSE) SetupOpenSuse ;;
     DEBIAN) SetupDebian ;;
     UBUNTU) SetupDebian ;;
     FEDORA) SetupFedora ;;
     CENTOS) SetupCentOS ;;
     *)
        echo "Your distribution $DIST is currently not supported by this script."
        echo "You are encouraged to contribute a new boostrap routine by taking"
        echo "the Setup* functions in $0 as a reference for implementing a function"
        echo "for $DIST. When finished please place your modified $0 into"
	echo "CFS_SOURCE_DIR/share/scripts and commit it to the Subversion repo"
	echo "or send it to one of the CFS++ developers."
        ;;
esac
