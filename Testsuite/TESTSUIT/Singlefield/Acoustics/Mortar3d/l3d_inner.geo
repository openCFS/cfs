/*
Input file for Gmsh for generating inner part of mesh.
*/

Point(1) = { 0.0, 0.0, 0.0 };
Point(2) = { 0.5, 0.0, 0.0 };
Point(3) = { 0.5, 0.5, 0.0 };
Point(4) = { 0.0, 0.5, 0.0 };

Line(5) = { 1, 2 };
Line(6) = { 2, 3 };
Line(7) = { 3, 4 };
Line(8) = { 4, 1 };

Line Loop(9) = {5, 6, 7, 8}; Plane Surface(10) = {9};

Characteristic Length { 1, 2, 3, 4 } = 0.40;

Extrude{0, 0, 0.5} { Surface{10}; }

Physical Surface ( "inner_iface1" ) = { 23 };
Physical Surface ( "inner_iface2" ) = { 27 };
Physical Surface ( "inner_iface3" ) = { 32 };

Physical Volume ( "inner" ) = { 1 };

