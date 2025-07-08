//define constants
a = 1;
c = 0.1;
dx = 0.02;
dy = 0.01;

//define points
Point(1) = {0, -c, 0};
Point(2) = { a, -c, 0};
Point(3) = { a,  c, 0};
Point(4) = {0,  c, 0};


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
Transfinite Line{6,8} = c/dy+1;

//Define structured meshing and recombine triangles to quadrilaterals
Transfinite Surface{6} = {1,2,3,4};
Recombine Surface{6};


//Generate regions to address later in CFS++
Physical Surface("S_inner") = {6};

Physical Line("bottom") = {5};
Physical Line("top") = {7};
Physical Line("left") = {8};
Physical Line("right") = {6};


