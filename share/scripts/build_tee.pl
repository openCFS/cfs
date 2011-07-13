#!/usr/bin/perl
#
# This script executes a command with two arguments and sends its output to
# stdout and a logfile. This is what the shell command "tee" does. However,
# we cannot use "tee", because we need the exit code of the actual command,
# not the exit code of "tee". So here is the Perl version...

# First 3 arguments form the command line
$cmd="$ARGV[0] $ARGV[1] $ARGV[2]";
# print command line (goes into CMake's logfile)
print "Command: $cmd\n";

# open our own logfile
open(my $logfile, ">", "$ARGV[3]") or die "Cannot open logfile '$ARGV[3]'.";

# execute command and open pipe to its stdout stream
open(my $child, "$cmd 2>&1 |") or die "Cannot execute '$ARGV[0]'.";

# print output to stdout and our logfile ("tee")
while (<$child>) {
  print STDOUT $_;    # $_ contains the current line read from pipe
  print $logfile $_;
}

# close pipe and wait until command terminates
close($child);
# if close sets $!=1 then there was something wrong with the pipe
if ($!) { warn "Error closing pipe."; }
# get rid of Perl's error codes (8 least significant bits)
# and obtain the exit code of our command
my $RC = $? >> 8;
# print the error code if it isn't "success"
print "Exit code: $RC" if $RC != 0;

close($logfile);

# give the correct exit code to CMake
exit $RC;
