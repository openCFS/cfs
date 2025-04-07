// GMSH geometry definition
// vim: set filetype=c : //

// specify the dimensions & discretisation

R_i = 20; // inner radius of pipe
R_a = 22; // outer radius of pipe
N_c = 50; // number of elements on radius
N_r = 5;  // number of elements in radial direction

Geometry.LineNumbers = 1; // show line numbers
Geometry.PointNumbers = 1; // show point numbers
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 1; // show surface numbers

Point(1) = {R_i, 0.0, 0.0, 5}; //
Point(2) = {R_a, 0.0, 0.0, 5}; //
Point(3) = {0.0, R_a, 0.0, 5};
Point(4) = {0.0, R_i, 0.0, 5};
Point(5) = {0.0, 0.0, 0.0, 5};
 
Line(1) = {1 , 2}; 
Circle(2) = {2, 5, 3};
Line(3) = {3, 4};
Circle(4) = {4, 5, 1};

Line Loop(10) = {1, 2, 3, 4};
Plane Surface(1) = {10};

Transfinite Line{1, 3} = N_r+1; // division in radial direction
Transfinite Line{2, 4} = N_c+1; // division in circular direction
Transfinite Surface {1}; // mesh the surface with transfinite algorithm

Recombine Surface {1};

Mesh.RecombineAll = 1; // recombine triangles to quads
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.ElementOrder = 2; // use second order elements
Mesh.SecondOrderIncomplete = 1; // use 8-noded quads

Physical Line("inner") = {4}; //inner side of pipe
Physical Line("outer") = {2}; //outer side of pipe
Physical Line("x-axis") = {1}; //line symmetric to the x-axis
Physical Line("y-axis") = {3}; //line symmetric to the y-axis
Physical Surface("pipe") = {1}; //pipe surface


