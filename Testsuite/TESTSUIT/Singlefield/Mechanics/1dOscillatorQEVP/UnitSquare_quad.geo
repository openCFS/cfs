l = 1.0;
n = 1;
d = l/n;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Point(1) = {-0.5*l, -0.5*l, 0, d};
Point(2) = { 0.5*l, -0.5*l, 0, d};
Point(3) = { 0.5*l,  0.5*l, 0, d};
Point(4) = {-0.5*l,  0.5*l, 0, d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,2,3,4} = n+1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};

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
// volume
Physical Surface('surf') = {1};
