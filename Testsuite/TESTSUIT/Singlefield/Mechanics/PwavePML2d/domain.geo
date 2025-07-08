// GMSH geometry definition
// vim: set filetype=c : //

// specify the dimensions & discretisation
a = 1.0;
h = 1.0/10.0;
l = 0.5;
N_pml = 5;//$PML;
d = h*N_pml;

// set view options
Geometry.LineNumbers = 1;
Geometry.PointNumbers = 1;
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 1; // show surface numbers
Geometry.VolumeNumbers = 1; // show volume numbers

Point(1) = { 0,  0,  0, h};
Point(2) = { a,  0,  0, h};
Point(3) = { a,  a, 0, h};
Point(4) = { 0,  a, 0, h};
Point(5) = { l,  0, 0, h};
Point(6) = { a+d, 0, 0, h};
Point(7) = { a+d, a, 0, h};
Point(8) = { a+d, a+d, 0, h};
Point(9) = { a, a+d, 0, h};
Point(10) = { 0, a+d, 0, h};
Point(11) = {0, l, 0, h};
Line(1) = {1, 5};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 11};
Line(5) = {5, 2};
Line(6) = {2, 6};
Line(7) = {3, 7};
Line(8) = {3, 9};
Line(9) = {4, 10};
Line(10) = {6, 7};
Line(11) = {7, 8};
Line(12) = {8, 9};
Line(13) = {9, 10};
Circle(14) = { 5, 1, 11};
Line(15) = {1,11};

Line Loop(1) = {5, 2, 3, 4,-14};
Plane Surface(1) = {1};
Line Loop(2) = {-2, -7, 10, 6};
Ruled Surface(2) = {2};
Line Loop(3) = {7, 11, 12, -8};
Ruled Surface(3) = {3};
Line Loop(4) = {13, -9, -3, 8};
Ruled Surface(4) = {4};
// inside circle
Line Loop(5) = {-15, 14, 1};
Plane Surface(5) = {5};

// regular mesh in PML domain
Transfinite Line{6,7,8,9,11,12} = N_pml+1;
Surf_pml[] = {2,3,4};
Transfinite Surface{Surf_pml[]};
Recombine Surface{Surf_pml[]};

// define regions
Physical Surface("pml") = {Surf_pml[]}; // top of the beam
Physical Surface("outside") = {1};
Physical Surface("inside") = {5};
Physical Line("abc") = {2,3};
Physical Line("x0") = {15,4,9};
Physical Line("y0") = {1,5,6};
Physical Line("excitation_line") = {1};
Physical Line("excitation_circle") = {14};
Physical Line("x") = {5,6};
Physical Line("y") = {4,9};
Physical Point("origin") = {1};
Physical Point("circle_point") = {5};


// specify mesh parameters
//Mesh.RecombineAll = 1; // recombine triangles to quads
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.ElementOrder = 2; // use second order elements
Mesh.SecondOrderIncomplete = 1; // use 20-noded hexas

