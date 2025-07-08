d = 5.0; // mesh size
h = 30.0; // hight
l = 150.0; // length
//t = 5.0; // absorbing thickness

// mesh options
Mesh.ElementOrder = 2;
Mesh.RecombineAll = 1;

Point(1) = { l/2, 0, 0, d};
Point(2) = {-l/2, 0, 0, d};
Point(3) = {-l/2,-h, 0, d};
Point(4) = { l/2,-h, 0, d};
//Point(5) = {-l/2,-h-t, 0, d};
//Point(6) = { l/2,-h-t, 0, d};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//Line(5) = {3, 5};
//Line(6) = {5, 6};
//Line(7) = {6, 4};

Line Loop(10) = {1, 2, 3, 4};
Plane Surface(10) = {10};
//Line Loop(11) = {-3, 5, 6, 7};
//Plane Surface(11) = {11};

Transfinite Line {1,3} = l/d + 1;
Transfinite Line {2,4} = h/d + 1;
//Transfinite Line {7,5} = t/d + 1;
//Transfinite Surface { 11,10 };
Transfinite Surface { 10 };

Physical Line('SURFACE') = {1};
//Physical Surface('INF') = {11};
Physical Surface('FLUID') = {10};
