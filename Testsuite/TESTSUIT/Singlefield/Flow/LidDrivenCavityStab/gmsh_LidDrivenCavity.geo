// Gmsh (http://geuz.org/gmsh) channel model  for Stokes flow according to PhD
// thesis of  Gerhard Link: A Finite Element  Scheme for Fluid–Solid–Acoustics
// Interactions and its Application to Human Phonation. pages 61 ff.

// Lid Driven Cavity
W = 1;
H = 1;

Point(1) = {-W, -H, 0};
Point(2) = { W, -H, 0};
Point(3) = { W,  H, 0};
Point(4) = {-W,  H, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface{6} = {};

Transfinite Line {2} = 11 Using Progression 1.0;
Transfinite Line {-4} = 11 Using Progression 1.0;
Transfinite Line {-1} = 11 Using Progression 1.0;
Transfinite Line {3} = 11 Using Progression 1.0;

Recombine Surface{6};

Physical Surface("cavity") = {6};

Physical Line("bottom") = 1;
Physical Line("top") = 3;
Physical Line("right") = 2;
Physical Line("left") = 4;
