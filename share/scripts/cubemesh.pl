#!/usr/bin/perl

$usage = "
Usage: cubemesh.pl [2d] <sizes>
 Options:
 2d    generate 2D mesh instead of 3D
Sizes:
 3D: n | n s | n l u | n sx sy sz | nx ny nz sx sy sz | n lx ly lz ux uy uz | nx ny nz lx ly lz ux uy uz
 2D: n | n s | nx ny sx sy | n lx ly ux uy | nx ny lx ly ux uy
 n, nx, ny, nz: elements in corresponding direction
 s, sx, sy, sz: size (default=1)
 l, lx, ly, lz: left-lower coordinate (default=0)
 u, ux, uy, uz: upper-right coordinate (default=s)

";

use warnings;

$dim = 3;
if(@ARGV > 0 && $ARGV[0] =~ /2d/i){
 $dim = 2;
 shift(@ARGV);
}

$lx = $ly = $lz = 0.0;
$ux = $uy = $uz = 1.0;

if($dim == 2){
 $nz = 0;
 if(@ARGV == 1 || @ARGV == 2 || @ARGV == 5){
  $nx = $ny = shift(@ARGV);
 }elsif(@ARGV == 4 || @ARGV == 6){
  $nx = shift(@ARGV);
  $ny = shift(@ARGV);
 }
 if(@ARGV == 1){
  $ux = $uy = shift(@ARGV);
 }elsif(@ARGV == 2){
  $ux = shift(@ARGV);
  $uy = shift(@ARGV);
 }elsif(@ARGV == 4){
  $lx = shift(@ARGV);
  $ly = shift(@ARGV);
  $ux = shift(@ARGV);
  $uy = shift(@ARGV);
 }
}else{ #3d
 if(@ARGV == 1 || @ARGV == 2 || @ARGV == 3 || @ARGV == 4 || @ARGV == 7){
  $nx = $ny = $nz = shift(@ARGV);
 }else{
  $nx = shift(@ARGV);
  $ny = shift(@ARGV);
  $nz = shift(@ARGV);
 }
 if(@ARGV == 1){
  $ux = $uy = $uz = shift(@ARGV);
 }elsif(@ARGV == 2){
  $lx = $ly = $lz = shift(@ARGV);
  $ux = $uy = $uz = shift(@ARGV);
 }elsif(@ARGV == 3){
  $ux = shift(@ARGV);
  $uy = shift(@ARGV);
  $uz = shift(@ARGV);
 }elsif(@ARGV == 6){
  $lx = shift(@ARGV);
  $ly = shift(@ARGV);
  $lz = shift(@ARGV);
  $ux = shift(@ARGV);
  $uy = shift(@ARGV);
  $uz = shift(@ARGV);
 }
}

if(!defined($nx) || @ARGV != 0){ # at least n was given and no parameters are left
 print($usage);
 exit(1);
}

$nnx = $nx + 1;
$nny = $ny + 1;
$nnz = $nz + 1;

print(STDERR $nx.",".$ny.",".$nz);
if($dim == 3){
 $elem3d = $nx*$ny*$nz;
 $elem2d = 2 * ($nx*$ny + $nx*$nz + $ny*$nz);
}else{
 $elem3d = 0;
 $elem2d = $nx*$ny;
}

# generate header
print("[Info]\n");
printf("Version %d\n", 1);
printf("Dimension %d\n", $dim);
printf("NumNodes %d\n", $nnx*$nny*$nnz);
printf("Num3DElements %d\n", $elem3d);
printf("Num2DElements %d\n", $elem2d);
printf("Num1DElements %d\n", 0);
printf("NumNodeBC %d\n", $dim == 3 ? 2 * ($nnx*$nny + $nnx*$nnz + $nny*$nnz) : 2*($nnx+$nny));
printf("NumSaveNodes %d\n", 0);
printf("NumSaveElements %d\n", 0);
printf("Num 2d-line      : %d\n", 0);
printf("Num 2d-line,quad : %d\n", 0);
printf("Num 3d-line      : %d\n", 0);
printf("Num 3d-line,quad : %d\n", 0);
printf("Num triangle     : %d\n", 0);
printf("Num triangle,quad: %d\n", 0);
printf("Num quadr        : %d\n", $elem2d);
printf("Num quadr,quad   : %d\n", 0);
printf("Num tetra        : %d\n", 0);
printf("Num tetra,quad   : %d\n", 0);
printf("Num brick        : %d\n", $elem3d);
printf("Num brick,quad   : %d\n", 0);
printf("Num pyramid      : %d\n", 0);
printf("Num pyramid,quad : %d\n", 0);
printf("Num wedge        : %d\n", 0);
printf("Num wedge,quad   : %d\n", 0);
print("\n");

print("[Nodes]\n");
print("#NodeNr x-coord y-coord z-coord\n");
# generate nodes
$c = 0;
$hx = ($ux - $lx) / $nx;
$hy = ($uy - $ly) / $ny;
$hz = $nz > 0.0 ? ($uz - $lz) / $nz : 0.0;
for($z = 0; $z < $nnz; ++$z){
 for($y = 0; $y < $nny; ++$y){
  for($x = 0; $x < $nnx; ++$x){
   printf("%8d %20.13e %20.13e %20.13e\n", ++$c, $lx+$hx*$x, $ly+$hy*$y, $lz+$hz*$z);
  }
 }
}
print("\n");

# return the node number
# note that the node number is 1+$x+$nnx*($y+$nny*$z) by construction
# note this is slow and therefore not used
#sub nodenr{
# my($x,$y,$z) = @_;
# return 1+$x+$nnx*($y+$nny*$z);
#}

# generate elements
$c = 0;
print("[1D Elements]\n");
print("#ElemNr  ElemType  NrOfNodes  Level\n");
print("#Node1 Node2 ... NodeNrOfNodes\n");
print("\n");

print("[3D Elements]\n");
print("#ElemNr  ElemType  NrOfNodes  Level\n");
print("#Node1 Node2 ... NodeNrOfNodes\n");
if($dim == 3){
 for($z = 0; $z < $nz; ++$z){
  for($y = 0; $y < $ny; ++$y){
   for($x = 0; $x < $nx; ++$x){
    printf("%d %d %d %s\n", ++$c, 10, 8, "mech");
#   printf("%d %d %d %d %d %d %d %d\n", nodenr($x,$y,$z+1), nodenr($x,$y+1,$z+1), nodenr($x+1,$y+1,$z+1), nodenr($x+1,$y,$z+1), nodenr($x,$y,$z), nodenr($x,$y+1,$z), nodenr($x+1,$y+1,$z), nodenr($x+1,$y,$z));
    printf("%d %d %d %d %d %d %d %d\n", (($z+1)*$nny+$y)*$nnx+$x+1, (($z+1)*$nny+$y+1)*$nnx+$x+1, (($z+1)*$nny+$y+1)*$nnx+$x+2, (($z+1)*$nny+$y)*$nnx+$x+2, ($z*$nny+$y)*$nnx+$x+1, ($z*$nny+$y+1)*$nnx+$x+1, ($z*$nny+$y+1)*$nnx+$x+2, ($z*$nny+$y)*$nnx+$x+2);
   }
  }
 }
}
print("\n");

print("[2D Elements]\n");
print("#ElemNr  ElemType  NrOfNodes  Level\n");
print("#Node1 Node2 ... NodeNrOfNodes\n");
# we have surface elements at all sides
# left (x=0)
for($z = 0; $z < $nz; ++$z){
 for($y = 0; $y < $ny; ++$y){
  printf("%d %d %d %s\n", ++$c, 6, 4, "leftsurface");
#  printf("%d %d %d %d\n", nodenr(0,$y,$z), nodenr(0,$y,$z+1), nodenr(0,$y+1,$z+1), nodenr(0,$y+1,$z));
  printf("%d %d %d %d\n", ($z*$nny+$y)*$nnx+1, (($z+1)*$nny+$y)*$nnx+1, (($z+1)*$nny+$y+1)*$nnx+1, ($z*$nny+$y+1)*$nnx+1);
 }
}
# right (x=1)
for($z = 0; $z < $nz; ++$z){
 for($y = 0; $y < $ny; ++$y){
  printf("%d %d %d %s\n", ++$c, 6, 4, "rightsurface");
#  printf("%d %d %d %d\n", nodenr($n,$y,$z), nodenr($n,$y+1,$z), nodenr($n,$y+1,$z+1), nodenr($n,$y,$z+1));
  printf("%d %d %d %d\n", ($z*$nny+$y)*$nnx+$nnx, ($z*$nny+$y+1)*$nnx+$nnx, (($z+1)*$nny+$y+1)*$nnx+$nnx, (($z+1)*$nny+$y)*$nnx+$nnx);
 }
}
# bottom (y=0)
for($z = 0; $z < $nz; ++$z){
 for($x = 0; $x < $nx; ++$x){
  printf("%d %d %d %s\n", ++$c, 6, 4, "bottomsurface");
#  printf("%d %d %d %d\n", nodenr($x,0,$z), nodenr($x,0,$z+1), nodenr($x+1,0,$z+1), nodenr($x+1,0,$z));
  printf("%d %d %d %d\n", $z*$nny*$nnx+$x+1, ($z+1)*$nny*$nnx+$x+1, ($z+1)*$nny*$nnx+$x+2, $z*$nny*$nnx+$x+2);
 }
}
# top (y=1)
for($z = 0; $z < $nz; ++$z){
 for($x = 0; $x < $nx; ++$x){
  printf("%d %d %d %s\n", ++$c, 6, 4, "topsurface");
#  printf("%d %d %d %d\n", nodenr($x,$n,$z), nodenr($x+1,$n,$z), nodenr($x+1,$n,$z+1), nodenr($x,$n,$z+1));
  printf("%d %d %d %d\n", ($z*$nny+$ny)*$nnx+$x+1, ($z*$nny+$ny)*$nnx+$x+2, (($z+1)*$nny+$ny)*$nnx+$x+2, (($z+1)*$nny+$ny)*$nnx+$x+1);
 }
}
# back (z=0)
if($dim == 3){
 $region = "backsurface";
}else{
 $region = "mech";
}
for($y = 0; $y < $ny; ++$y){
 for($x = 0; $x < $nx; ++$x){
  printf("%d %d %d %s\n", ++$c, 6, 4, $region);
#  printf("%d %d %d %d\n", nodenr($x,$y,0), nodenr($x+1,$y,0), nodenr($x+1,$y+1,0), nodenr($x,$y+1,0));
  printf("%d %d %d %d\n", $y*$nnx+$x+1, $y*$nnx+$x+2, ($y+1)*$nnx+$x+2, ($y+1)*$nnx+$x+1);
 }
}
if($dim == 3){ # skip this for 2D
 # front (z=1)
 for($y = 0; $y < $ny; ++$y){
  for($x = 0; $x < $nx; ++$x){
   printf("%d %d %d %s\n", ++$c, 6, 4, "frontsurface");
 #  printf("%d %d %d %d\n", nodenr($x,$y,$n), nodenr($x,$y+1,$n), nodenr($x+1,$y+1,$n), nodenr($x+1,$y,$n));
   printf("%d %d %d %d\n", ($nz*$nny+$y)*$nnx+$x+1, ($nz*$nny+$y+1)*$nnx+$x+1, ($nz*$nny+$y+1)*$nnx+$x+2, ($nz*$nny+$y)*$nnx+$x+2);
  }
 }
}
print("\n");

print("[Node BC]\n");
print("#NodeNr Level\n");
# we name the nodes at all sides, we have (n+1)^2 nodes at every side, nodes are numbered z,y,x (from slowest to fastest increasing)
# note that the node number is $x+($n+1)*($y+($n+1)*$z)
# left (x=0)
for($z = 0; $z < $nnz; ++$z){
 for($y = 0; $y < $nny; ++$y){
#  printf("%8d %s\n", nodenr(0,$y,$z), "left");
  printf("%8d %s\n", ($z*$nny+$y)*$nnx+1, "left");
 }
}
# right (x=1)
for($z = 0; $z < $nnz; ++$z){
 for($y = 0; $y < $nny; ++$y){
#  printf("%8d %s\n", nodenr($n,$y,$z), "right");
  printf("%8d %s\n", ($z*$nny+$y)*$nnx+$nx+1, "right");
 }
}
# bottom (y=0)
for($z = 0; $z < $nnz; ++$z){
 for($x = 0; $x < $nnx; ++$x){
#  printf("%8d %s\n", nodenr($x,0,$z), "bottom");
  printf("%8d %s\n", ($z*$nny)*$nnx+$x+1, "bottom");
 }
}
# top (y=1)
for($z = 0; $z < $nnz; ++$z){
 for($x = 0; $x < $nnx; ++$x){
#  printf("%8d %s\n", nodenr($x,$n,$z), "top");
  printf("%8d %s\n", ($z*$nny+$ny)*$nnx+$x+1, "top");
 }
}
if($dim == 3){ # skip back/front in 3d
 # back (z=0)
 for($y = 0; $y < $nny; ++$y){
  for($x = 0; $x < $nnx; ++$x){
 #  printf("%8d %s\n", nodenr($x,$y,0), "back");
   printf("%8d %s\n", $y*$nnx+$x+1, "back");
  }
 }
 # front (z=1)
 for($y = 0; $y < $nny; ++$y){
  for($x = 0; $x < $nnx; ++$x){
 #  printf("%8d %s\n", nodenr($x,$y,$n), "front");
   printf("%8d %s\n", ($nz*$nny+$y)*$nnx+$x+1, "front");
  }
 }
}
print("\n");

print("[Save Nodes]\n");
print("#NodeNr Level\n");
print("\n");

print("[Save Elements]\n");
print("#ElemNr Level\n");
print("\n");

