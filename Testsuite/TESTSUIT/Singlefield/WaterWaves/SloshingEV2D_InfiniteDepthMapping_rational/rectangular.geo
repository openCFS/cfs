// gmsh -3 SloshingEV2D_InfiniteDepth.geo
d_x = 3; // mesh size in x-direction
h = 20; // hight
l = 150.0; // length
N_y = 5; // number of elements in depth direction

// mesh options
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Point(1) = { l/2, 0, 0};
Point(2) = {-l/2, 0, 0};
Point(3) = {-l/2,-h, 0};
Point(4) = { l/2,-h, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(10) = {1, 2, 3, 4};
Plane Surface(10) = {10};

Transfinite Line {1,3} = l/d_x + 1;
Transfinite Line {2,4} = N_y + 1;
Transfinite Surface { 10 };

Physical Line('SURFACE') = {1};
Physical Line('LEFT') = {2};
Physical Line('BOTTOM') = {3};
Physical Surface('FLUID') = {10};
