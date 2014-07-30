#!/usr/bin/perl

if(@ARGV != 3){
 die "Syntax: voidmesh.pl <density file> <mesh file> <threshold>\n\n";
}

open(DENSITY, '<', $ARGV[0]) or die $!;
while(<DENSITY>){
 if(/^[\t ]*<element[^>]* nr="([0-9]+)"[^>]* design="([0-9.eE+-]+)"[^>]*\/>[\t ]*$/){
  $value[$1] = $2 > $ARGV[2];
 }
}
close(DENSITY);

$elements3d = false;

open(MESH, '<', $ARGV[1]) or die $!;
while(<MESH>){
 $elements3d = true if /^[3D Elements]$/;
 $elements3d = false if /^$/;
 if($elements3d && /^([0-9]+)( [0-9]+ [0-9]+ )([a-zA-Z_].*)$/){
  print $1.$2.($value[$1]?$3:"void")."\n";
 }else{
  print;
 }
}
close(MESH);

