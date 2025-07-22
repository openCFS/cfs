l = 1.0;
a = 1/3;
lx = l;
ly = a*l;
n = 12;
d = 1.0;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Point(1) = {-0.5*lx, -0.5*ly, 0, d};
Point(2) = { 0.5*lx, -0.5*ly, 0, d};
Point(3) = { 0.5*lx,  0.5*ly, 0, d};
Point(4) = {-0.5*lx,  0.5*ly, 0, d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = n+1;
Transfinite Line{2,4} = n*a;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};

//Extrude {0, 0, l} {Surface{1};Layers{1};Recombine;}

// corner points
Physical Point('NO') = {3};
Physical Point('NW') = {4};
Physical Point('SO') = {2};
Physical Point('SW') = {1};
// coner edges
Physical Line('N') = {3};
Physical Line('O') = {2};
Physical Line('S') = {1};
Physical Line('W') = {4};
// faces
Physical Surface('surf') = {1};
