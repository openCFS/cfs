l = 1.0;
n = 6;
d = l/n;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Point(1) = {-0.5*l, -0.5*l, -0.5*l, d};
Point(2) = { 0.5*l, -0.5*l, -0.5*l, d};
Point(3) = { 0.5*l,  0.5*l, -0.5*l, d};
Point(4) = {-0.5*l,  0.5*l, -0.5*l, d};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,2,3,4} = n+1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};

Extrude {0, 0, l} {Surface{1};Layers{n};Recombine;}

// corner points
Physical Point('NOB') = {3};
Physical Point('NWB') = {4};
Physical Point('SOB') = {2};
Physical Point('SWB') = {1};
Physical Point('NOT') = {10};
Physical Point('NWT') = {14};
Physical Point('SOT') = {6};
Physical Point('SWT') = {5};
// coner edges
Physical Line('NW') = {20};
Physical Line('SW') = {11};
Physical Line('SE') = {12};
Physical Line('NO') = {16};
Physical Line('NB') = {3};
Physical Line('OB') = {2};
Physical Line('SB') = {1};
Physical Line('WB') = {4};
Physical Line('NT') = {8};
Physical Line('OT') = {7};
Physical Line('ST') = {6};
Physical Line('WT') = {9};
// faces
Physical Surface('N') = {21};
Physical Surface('O') = {17};
Physical Surface('S') = {13};
Physical Surface('W') = {25};
Physical Surface('B') = {1};
Physical Surface('T') = {26};
// volume
Physical Volume('SOLID') = {1};

// rotate
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/6} { Volume{1}; }
