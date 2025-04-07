//define constants
a = 5;
r = 0.1;
//Number of elements
N = 50;

//define points Rectangle large
Point(1) = {-a, -a, 0};
Point(2) = { a, -a, 0};
Point(3) = { a,  a, 0};
Point(4) = {-a,  a, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//define points hole small + origin
Point(5) = {0, 0, 0};
Point(6) = {-r, -r, 0};
Point(7) = { r, -r, 0};
Point(8) = { r,  r, 0};
Point(9) = {-r,  r, 0};
Circle(5) = {6, 5, 7};
Circle(6) = {7, 5, 8};
Circle(7) = {8, 5, 9};
Circle(8) = {9, 5, 6};

//define block lines
Line(9) = {1, 6};
Line(10) = {2, 7};
Line(11) = {3, 8};
Line(12) = {4, 9};

//define lines growth ratios
Transfinite Line{1,2,3,4,5,6,7,8} = N+1;
Transfinite Line{9,10,11,12} = N+1 Using Progression 0.9;

//Connect lines to surfaces
Line Loop(13) = {9, 5, -10, -1};
Plane Surface(14) = {13};
Line Loop(15) = {10, 6, -11, -2};
Plane Surface(16) = {15};
Line Loop(17) = {11, 7, -12, -3};
Plane Surface(18) = {17};
Line Loop(19) = {12, 8, -9, -4};
Plane Surface(20) = {19};

//Define structured meshing and recombine triangles to quadrilaterals
Transfinite Surface {14};
Transfinite Surface {16};
Transfinite Surface {18};
Transfinite Surface {20} Right;
Recombine Surface {14};
Recombine Surface {16};
Recombine Surface {18};
Recombine Surface {20};


Physical Surface("internal") = {14,16,18,20};
Physical Line("bottom") = {1};
Physical Line("left") = {4};
Physical Line("top") = {3};
Physical Line("right") = {2};





