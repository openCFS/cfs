// simple trafo
a = 1.0;
b1 = 1.0;
b2 = 1.0;
b3 = b1;
t = 1.0;
ny1 = 10;
ny2 = 10;
ny3 = 10;
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
//+
Extrude {0, 0, t} {
  Surface{13}; Surface{9}; Surface{5}; Layers{nz}; Recombine;
}
//+
Physical Volume("V1") = {3};
//+
Physical Volume("V2") = {2};
//+
Physical Volume("V3") = {1};
//+
Physical Surface("S1S") = {66};
//+
Physical Surface("S1O") = {70};
//+
Physical Surface("S1W") = {78};
//+
Physical Surface("S2W") = {56};
//+
Physical Surface("S3W") = {34};
//+
Physical Surface("S3O") = {26};
//+
Physical Surface("S2O") = {48};
//+
Physical Surface("S12") = {70, 44};
//+
Physical Surface("S23") = {22};
//+
Physical Surface("S3N") = {30};
//+
Physical Surface("S1T") = {5};
//+
Physical Surface("S2T") = {9};
//+
Physical Surface("S3T") = {13};
//+
Physical Surface("S3B") = {35};
//+
Physical Surface("S2B") = {57};
//+
Physical Surface("S1B") = {79};
