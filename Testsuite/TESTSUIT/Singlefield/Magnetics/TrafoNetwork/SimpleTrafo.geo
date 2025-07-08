// simple trafo
a = 1.0;
b1 = 1.0;
b2 = 1.0;
b3 = b2;
b4 = b1;
t = 1.0;
ny1 = 1;
ny2 = 1;
ny3 = 1;
ny4 = 1;
nx = 1;
nz = 1;


//+
Point(1) = {-a/2, -b1/2, -t/2, a/nx};
//+
Point(2) = {a/2, -b1/2, -t/2, a/nx};
//+
Line(1) = {1, 2};
//+
Extrude {0, b1, 0} {
  Line{1}; Layers{ny1}; Recombine;
}
//+
Extrude {0, b2, 0} {
  Line{2}; Layers{ny2}; Recombine;
}
//+
Extrude {0, b3, 0} {
  Line{6}; Layers{ny3}; Recombine;
}
Extrude {0, b4, 0} {
  Line{10}; Layers{ny4}; Recombine;
}
//+
Extrude {0, 0, t} {
  Surface{13}; Surface{9}; Surface{5}; Surface{17};  Layers{nz}; Recombine;
}
//+
Physical Volume("V1") = {3};
//+
Physical Volume("V2") = {2};
//+
Physical Volume("V3") = {1};
Physical Volume("V4") = {4};
// x+ = O, x- = W
// y+ = N, y- = S
// z+ = T, z- = B
Physical Surface("S_T") = {83, 61, 39, 105};
Physical Surface("S_B") = {5, 9, 13, 17};
Physical Surface("S_N") = {100};
Physical Surface("S_S") = {70};
Physical Surface("S1_W") = {82};
