//define constants
a = 3;
dx = 0.1;
dy = 0.1;

//define points
Point(1) = {-a, -a, 0};
Point(2) = { a, -a, 0};
Point(3) = { a,  a, 0};
Point(4) = {-a,  a, 0};

//define lines
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

//define surfaces
Line Loop(5) = {5, 6, 7, 8};
Plane Surface(6) = {5};

//define lines growth ratios
Transfinite Line{5,7} = a/dx+1;
Transfinite Line{6,8} = a/dy+1;

//Define structured meshing and recombine triangles to quadrilaterals
Transfinite Surface{6} = {1,2,3,4};
Recombine Surface{6};

//Generate regions to address later in CFS++
Physical Surface("fluid") = {6};

Physical Line("bottom") = {5};
Physical Line("top") = {7};
Physical Line("inlet") = {8};
Physical Line("outlet") = {6};

//perform 3D meshing with linear elements
Mesh.ElementOrder = 1;
Mesh 2;
