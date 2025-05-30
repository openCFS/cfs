// Gmsh  (http://geuz.org/gmsh)  model  for  a  simple model  of  first  order
// absorbing boundary conditions
//
// To create the mesh run
// gmsh gmsh_Abc2d.geo -2 -v 0 -format msh -order 1 -o gmsh_Abc2d.msh

W = 15e-2;

esize = 1e-2;

div_edge = 31;

Point(1) = {-W, -W, 0, esize};
Point(2) = { W, -W, 0, esize};
Point(3) = { W,  W, 0, esize};
Point(4) = {-W,  W, 0, esize};

Line(11) = {1, 2};
Line(12) = {2, 3};
Line(13) = {3, 4};
Line(14) = {4, 1};
Line Loop(30) = {11, 12, 13, 14};
Plane Surface(31) = {30};

Transfinite Line {11} = div_edge Using Progression 1.0;
Transfinite Line {12} = div_edge Using Progression 1.0;
Transfinite Line {13} = div_edge Using Progression 1.0;
Transfinite Line {14} = div_edge Using Progression 1.0;

Transfinite Surface{31} = {};

Recombine Surface{31};

Physical Surface("fluid") = {31};
Physical Line("dombnd") = {11, 12, 13, 14};

  
