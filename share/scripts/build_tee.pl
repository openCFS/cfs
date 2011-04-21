#!/usr/bin/perl

print "$ARGV[0]\n";
print "$ARGV[1]\n";
print "$ARGV[2]\n";
print "$ARGV[3]\n";

`"$ARGV[0]" "$ARGV[1]" "$ARGV[2]" | tee "$ARGV[3]"`

