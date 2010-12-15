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
sudo /opt/local/bin/port selfupdate

# Install required packages
sudo /opt/local/bin/port install doxygen graphviz teTeX texlive_texmf-minimal wget 
sudo /opt/local/bin/port install git-core +svn

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

perl /tmp/adjustenv.pl

rm -f /tmp/adjustenv.pl