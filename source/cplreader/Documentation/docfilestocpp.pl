#!/usr/bin/perl

# Hex dump program to generate C++ from arbitrary files
#
# Written by Peter N Lewis a long time ago
# Released in to the Public Domain
# modified by Rudif c/o Perlmonks.org, to handle CRLF conversion
#
# Used to generate cplreader documentation
# Usage: ./hexdump.pl test.cc GetHTMLDocumentation '*.html *.png'

use strict;
use warnings;

# usage() if $ARGV[0] and $ARGV[0] =~ m!^-[^-]!;

our $filepos = 0;
our $linechars = '';
our $filesize = 0;
our $actpos = 1;

#    my @files = glob "*.html *.png";
my $filename;
my $i = 1;
my $outfile="$ARGV[0]";
my $fncname="$ARGV[1]";
my @files = glob "$ARGV[2]";

open OUTFILE, "> $outfile";

print( OUTFILE "#include \"$fncname.hh\"\n\n");

print( OUTFILE "namespace CoupledField {\n\n");

print( OUTFILE "  void $fncname(std::map<std::string, std::vector<unsigned char> >& doc) {\n\n");

foreach $filename (@files) {
    $filepos = 0;
    $linechars = '';
    $actpos = 1;

    open FILE, '<:raw', $filename or die "no such file $_";
    seek FILE, 0, 2;
    $filesize = tell FILE;
    seek FILE, 0, 0;

    print( OUTFILE "    static unsigned char doc_file${i}[] =\n    {\n" );
    while (<FILE>) {
	dump_char($_, 0, \*OUTFILE) foreach split(//);
    }
    print( OUTFILE "    };\n\n" );
    print( OUTFILE "    std::copy( doc_file$i, doc_file$i + sizeof(doc_file$i), std::back_inserter(doc[\"$filename\"]) );\n\n" );
#    dump_char( ' ', 1 ) while length($linechars) != 0;
    close FILE;

    $i++;
}

print( OUTFILE "  } // end of namespace\n\n");
print( OUTFILE "}\n");

close OUTFILE;

sub dump_char {
    my ( $char, $blank, $OUTFILE ) = @_;

#    if ( length( $linechars ) == 0 ) {
#	printf( "%06X: ", $filepos );
#    }
    if ( length( $linechars ) == 0 ) {
	print( $OUTFILE "      " );
    }

    $linechars .= ( $char =~ m#[!-~ ]# ) ? $char : '.';
    if ( $blank ) {
	print $OUTFILE '   ';
    } else {
#	print( "$filesize $actpos " );
	if($actpos != $filesize)
	{
	    printf( $OUTFILE "0x%02X, ", ord($char) );
	}
	else
	{
	    printf( $OUTFILE "0x%02X\n", ord($char) );
	}
	$actpos++;
    }
    print $OUTFILE ' ' if length( $linechars ) % 4 == 0;
    if ( length( $linechars ) == 12 ) {
#	print( $linechars, "\n" );
	print( $OUTFILE "\n" );
	$linechars = '';
	$filepos += 12;
    }
}

#sub usage {
#    print STDERR <<EOM;
#Usage: hdump.pl [file]...
#    Example `hdump.pl .cshrc' or `ls -l | hdump.pl'
#EOM
#  exit( 0 );
#}
