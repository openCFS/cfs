//define constants
a = 10;
c = 0.5;
b = 2;
dx = 0.2;
dy = 0.06;
dxB = 0.2;

//define points
Point(1) = {0, -c, 0};
Point(2) = { a, -c, 0};
Point(3) = { a,  c, 0};
Point(4) = {0,  c, 0};
Point(5) = { a+b,  c, 0};
Point(6) = { a+b,  -c, 0};


//define lines
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {3, 5};
Line(10) = {5, 6};
Line(11) = {6, 2};


//define surfaces
Line Loop(5) = {5, 6, 7, 8};
Plane Surface(6) = {5};

Line Loop(7) = {9, 10, 11, 6};
Plane Surface(8) = {7};

//define lines growth ratios
Transfinite Line{5,7} = a/dx+1;
Transfinite Line{6,8,10} = c/dy+1;
Transfinite Line{9,11} = b/dxB+1;

//Define structured meshing and recombine triangles to quadrilaterals
Transfinite Surface{6} = {1,2,3,4};
Recombine Surface{6};
Transfinite Surface{8} = {3,5,6,2};
Recombine Surface{8};

//Generate regions to address later in CFS++
Physical Surface("S_inner") = {6};
Physical Surface("S_Inf") = {8};

Physical Line("bottom") = {5};
Physical Line("bottom_Inf") = {11};
Physical Line("top") = {7};
Physical Line("top_Inf") = {9};
Physical Line("left") = {8};
Physical Line("right") = {10};
Physical Line("rightInner") = {6};

