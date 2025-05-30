// GMSH geometry definition
// vim: set filetype=c : //

// specify the dimensions & discretisation


// set view options
Geometry.LineNumbers = 1;
Geometry.PointNumbers = 1;
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 1; // show surface numbers
Geometry.VolumeNumbers = 1; // show volume numbers

Lx = 0.3;
Ly = Lx/10;
Lz = Lx/3;
a = 0.001;

Point(1) = {0,0,0};
Point(2) = {0,a,0};
Point(3) = {0,a,a};
Point(4) = {0,0,a};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};
Surface(1) = {5};
Transfinite Line{1,2,3,4} = 1;
out[] = Extrude {Lx,Ly,Lz} {Surface{1}; Layers{200}; Recombine;};
Physical Surface("S_in") = {1};
Physical Surface("S_out") = {out[0]};
Physical Surface("S_side") = {out[2],out[3],out[4],out[5]};
Physical Volume("V") = {out[1]};


// specify mesh parameters
Mesh.RecombineAll = 1; // recombine triangles to quads
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.ElementOrder = 2; // use second order elements
Mesh.SecondOrderIncomplete = 1; // use 20-noded hexas
