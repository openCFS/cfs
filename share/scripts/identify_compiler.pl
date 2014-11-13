#!/usr/bin/perl

use feature qw(switch);

$COMPILER=$ARGV[0];
$SOURCE=$ARGV[1];
$OUTPUT=$ARGV[2];

# print "CC $CC\nSOURCE $SOURCE\nOUTPUT $OUTPUT\n";

foreach $line (`$COMPILER -E $SOURCE 2>&1`) {

    if( $line =~ /CC_ID/ ) {
	$line =~ s/CC_ID //;
	chomp($line);
	$CC_ID=$line;
    }

    if( $line =~ /^CC_VERSION/ ) {
	$line =~ s/CC_VERSION //;
	chomp($line);
	$CC_VERSION=$line;
    }

    if( $line =~ /^CC_GCC_VERSION/ ) {
	$line =~ s/CC_GCC_VERSION //;
	chomp($line);
	$CC_GCC_VERSION=$line;
    }

    if( $line =~ /CXX_ID/ ) {
	$line =~ s/CXX_ID //;
	chomp($line);
	$CXX_ID=$line;
    }

    if( $line =~ /^CXX_VERSION/ ) {
	$line =~ s/CXX_VERSION //;
	chomp($line);
	$CXX_VERSION=$line;
    }

    if( $line =~ /^CXX_GCC_VERSION/ ) {
	$line =~ s/CXX_GCC_VERSION //;
	chomp($line);
	$CXX_GCC_VERSION=$line;
    }


    if( $line =~ /FC_ID/ ) {
	$line =~ s/FC_ID //;	
	chomp($line);
	$FC_ID=$line;
    }

    if( $line =~ /^FC_VERSION/ ) {
	$line =~ s/FC_VERSION //;
	chomp($line);
	$FC_VERSION=$line;
    }
}   

if ($CC_ID ne "") {
    given ($CC_ID) {
	when ("GCC") { 
	    $CC_VERSION =~ s/ //g;
	    $CC_GCC_VERSION =~ s/ //g; 
	}
	when ("ICC") { 
	    $ICC_VER=$CC_VERSION;
	    $ICC_VER =~ s/ [0-9].*$//;
	    {
		use integer;
		$MAJOR_VER = $ICC_VER / 100;
		$MINOR_VER = $ICC_VER % 100;
		$MINOR_VER = $MINOR_VER / 10;
	    }
	    $ICC_DATE=$CC_VERSION;
	    $ICC_DATE =~ s/^[0-9].* //;
	    $CC_VERSION = "$MAJOR_VER.$MINOR_VER $ICC_DATE";
	    $CC_GCC_VERSION =~ s/ //g; 
	}
	when ("OPEN64") {
	    $CC_VERSION =~ s/ //g;
	    $CC_GCC_VERSION =~ s/ //g; 
	}
	when ("MSVC") {
	    $CC_VERSION =~ s/ //g;
	}
	default {
	    print "C compiler not supported!\n";
	    exit 1;
	}
    }

    given ($OUTPUT) {
	when ("cmake") { 
	    print "SET(CC_ID \"$CC_ID\")\n"; 
	    print "SET(CC_VERSION \"$CC_VERSION\")\n"; 
	    print "SET(CC_GCC_VERSION \"$CC_GCC_VERSION\")\n";
	    exit 0;
	}
	when ("perl") { 
	    print "\$CC_ID=\"$CC_ID\";\n"; 
	    print "\$CC_VERSION=\"$CC_VERSION\";\n"; 
	    print "\$CC_GCC_VERSION=\"$CC_GCC_VERSION\";\n";
	    exit 0;
	}
    }
} 

if ($CXX_ID ne "") {
    given ($CXX_ID) {
	when ("GCC") { 
	    $CXX_VERSION =~ s/ //g;
	    $CXX_GCC_VERSION =~ s/ //g; 
	}
	when ("ICC") {
	    # Intel 12.x define __ICC and __INTEL_COMPILER to be 9999 therefore
            # we read the output of $COMPILER --version
	    ($version, $build_date) = split(/ /, $CXX_VERSION, 2);
	    if($version eq "9999") {
		foreach $line (`$COMPILER --version 2>&1`) {
		    chomp($line);
		    ($icc, $ICC_UPPER, $version, $date) = split(/ /, $line, 4);
		    if($icc eq "icpc") {
			$version =~ s/\.//g;
			$CXX_VERSION="$version $date";
		    }
		}
	    }


	    $ICC_VER=$CXX_VERSION;
	    $ICC_VER =~ s/ [0-9].*$//;
	    {
		use integer;
		$MAJOR_VER = $ICC_VER / 100;
		$MINOR_VER = $ICC_VER % 100;
		$MINOR_VER = $MINOR_VER / 10;
	    }
	    $ICC_DATE=$CXX_VERSION;
	    $ICC_DATE =~ s/^[0-9].* //;
	    $CXX_VERSION = "$MAJOR_VER.$MINOR_VER $ICC_DATE";
	    $CXX_GCC_VERSION =~ s/ //g; 
	}
	when ("OPEN64") {
	    $CXX_VERSION =~ s/ //g;
	    $CXX_GCC_VERSION =~ s/ //g; 
	}
	when ("MSVC") {
	    $CC_VERSION =~ s/ //g;
	}
	default {
	    print "C++ compiler not supported!\n";
	    exit 1;
	}
    }

    given ($OUTPUT) {
	when ("cmake") { 
	    print "SET(CXX_ID \"$CXX_ID\")\n"; 
	    print "SET(CXX_VERSION \"$CXX_VERSION\")\n"; 
	    print "SET(CXX_GCC_VERSION \"$CXX_GCC_VERSION\")\n";
	    exit 0;
	}
	when ("perl") { 
	    print "\$CXX_ID=\"$CXX_ID\";\n"; 
	    print "\$CXX_VERSION=\"$CXX_VERSION\";\n"; 
	    print "\$CXX_GCC_VERSION=\"$CXX_GCC_VERSION\";\n";
	    exit 0;
	}
    }
}

if ($FC_ID ne "") {
    given ($FC_ID) {
	when ("GNU") { 
	    $FC_VERSION =~ s/ //g;
	}
	when ("IFORT") { 
	    # Intel 12.x define __ICC and __INTEL_COMPILER to be 9999 therefore
            # we read the output of $COMPILER --version
	    ($version, $build_date) = split(/ /, $FC_VERSION, 2);
	    if($version eq "9999") {
		foreach $line (`$COMPILER --version 2>&1`) {
		    chomp($line);
		    ($icc, $ICC_UPPER, $version, $date) = split(/ /, $line, 4);
		    if($icc eq "ifort") {
			$version =~ s/\.//g;
			$FC_VERSION="$version $date";
		    }
		}
	    }

	    $ICC_VER=$FC_VERSION;
	    $ICC_VER =~ s/ [0-9].*$//;
	    {
		use integer;
		$MAJOR_VER = $ICC_VER / 100;
		$MINOR_VER = $ICC_VER % 100;
		$MINOR_VER = $MINOR_VER / 10;
	    }
	    $ICC_DATE=$FC_VERSION;
	    $ICC_DATE =~ s/^[0-9].* //;
	    $FC_VERSION = "$MAJOR_VER.$MINOR_VER $ICC_DATE";
	}
	when ("OPEN64") {
	    $FC_VERSION =~ s/ //g;
	}
	default {
	    print "Fortran compiler not supported!\n";
	    exit 1;
	}
    }

    given ($OUTPUT) {
	when ("cmake") { 
	    print "SET(FC_ID \"$FC_ID\")\n"; 
	    print "SET(FC_VERSION \"$FC_VERSION\")\n"; 
	    exit 0;
	}
	when ("perl") { 
	    print "\$FC_ID=\"$FC_ID\";\n"; 
	    print "\$FC_VERSION=\"$FC_VERSION\";\n"; 
	    exit 0;
	}
    }
}

print "Compiler not supported!\n";
exit 1;

