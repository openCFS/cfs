// GMSH geometry definition
// vim: set filetype=c : //

// specify the dimensions & discretisation
L = 1.2; // Length of channel (both directions)
h = 0.02; // mesh size
a = h; // height of channel

// set view options
Geometry.LineNumbers = 1;
Geometry.PointNumbers = 1;
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 1; // show surface numbers
Geometry.VolumeNumbers = 1; // show volume numbers

Point(1) = { 0,  0,  0, h};
Point(2) = { 0,  a,  0, h};
Line(1) = {1, 2};

// regular mesh in PML domain
Transfinite Line{1} = 1; // one element in thickness direction

// extrude to left and right
outR[] = Extrude{ L,0,0 } { Line{1}; Layers{L/h}; Recombine; };
//outL[] = Extrude{ -L,0,0 } { Line{1}; Layers{L/h}; Recombine; };

// define regions
Physical Surface("S_flow") = { outR[1] }; // top of the beam
Physical Line("L_left") = { 1 };
Physical Line("L_right") = {outR[0]};
Physical Line("L_top") = {outR[2] };
Physical Line("L_bottom") = {-outR[3]};
