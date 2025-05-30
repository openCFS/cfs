l = 1.0;
d = 0.5;

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

Transfinite Line{1,3} = 4;
Transfinite Line{2,4} = 4;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface{1};

// corner points
Physical Point('N_NO') = {3};
Physical Point('N_NW') = {4};
Physical Point('N_SO') = {2};
Physical Point('N_spring') = {1};

// corner edges
Physical Line('S_N') = {3};
Physical Line('S_O') = {2};
Physical Line('S_fixedy') = {1};
Physical Line('S_fixedx') = {4};

// domain
Physical Surface('V_square') = {1};

// export mesh 
Save "UnitSquareRefined.msh";

