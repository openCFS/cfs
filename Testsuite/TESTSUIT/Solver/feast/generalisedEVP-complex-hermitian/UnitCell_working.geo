l = 1.0;
a_x = 1/3;
a_y = 1/3;
n = 1;
d = l/n;
Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;
//Mesh.Algorithm3D = 4; // Frontal

Point(1) = {  0.5*l,  0.5*l, 0.0, d};
Point(2) = { -0.5*l,  0.5*l, 0.0, d};
Point(3) = { -0.5*l, -0.5*l, 0.0, d};
Point(4) = {  0.5*l, -0.5*l, 0.0, d};
Point(5)  = {  0.5*l*a_x,  0.5*l, 0.0, d};
Point(6)  = { -0.5*l*a_x,  0.5*l, 0.0, d};
Point(7)  = {  0.5*l*a_x,  0.5*l*a_y, 0.0, d};
Point(8)  = { -0.5*l*a_x,  0.5*l*a_y, 0.0, d};
Point(9)  = {  0.5*l*a_x,  -0.5*l*a_y, 0.0, d};
Point(10) = { -0.5*l*a_x,  -0.5*l*a_y, 0.0, d};
Point(11) = {  0.5*l*a_x,  -0.5*l, 0.0, d};
Point(12) = { -0.5*l*a_x,  -0.5*l, 0.0, d};
Point(13) = { 0.5*l,  0.5*l*a_y, 0.0, d};
Point(14) = { 0.5*l, -0.5*l*a_y, 0.0, d};
Point(15) = {-0.5*l,  0.5*l*a_y, 0.0, d};
Point(16) = {-0.5*l, -0.5*l*a_y, 0.0, d};

Line(1) = {1, 5};
Line(2) = {5, 7};
Line(3) = {7, 9};
Line(4) = {9, 11};
Line(5) = {11, 4};
Line(7) = {5, 6};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {8, 10};
Line(11) = {10, 9};
Line(12) = {10, 12};
Line(13) = {12, 11};
Line(14) = {12, 3};
Line(16) = {2, 6};
Line(17) = {2, 15};
Line(18) = {15, 8};
Line(19) = {15, 16};
Line(20) = {16, 10};
Line(21) = {16, 3};
Line(22) = {9, 14};
Line(23) = {14, 4};
Line(24) = {14, 13};
Line(25) = {7, 13};
Line(26) = {13, 1};

Line Loop(27) = {16, 8, -18, -17};
Plane Surface(28) = {27};
Line Loop(29) = {7, 8, 9, -2};
Plane Surface(30) = {29};
Line Loop(31) = {1, 2, 25, 26};
Plane Surface(32) = {31};
Line Loop(33) = {18, 10, -20, -19};
Plane Surface(34) = {33};
Line Loop(35) = {20, 12, 14, -21};
Plane Surface(36) = {35};
Line Loop(37) = {11, 4, -13, -12};
Plane Surface(38) = {37};
Line Loop(39) = {4, 5, -23, -22};
Plane Surface(40) = {39};
Line Loop(41) = {25, -24, -22, -3};
Plane Surface(42) = {41};
Line Loop(43) = {9, 3, -11, -10};
Plane Surface(44) = {43};

//Extrude {0, 0, 1} {
//  Surface{30, 28, 34, 36, 38, 44, 40, 42, 32};
//  Layers{n};
//  Recombine;
//}
//
Transfinite Line{7,9,11,13,136,134,48,46} = n*a_x/l + 1 ; // center x
Transfinite Line{24,3,10,19,93,91,157,201} = n*a_y/l + 1; // center y
Transfinite Line{225,49,47,71,26,2,8,17,180,135,113,115,23,4,12,21} = n*(1-a_y)/l/2 + 1; // outside y
Transfinite Line{68,70,92,114,16,18,20,14,5,22,25,1,222,200,181,179} = n*(1-a_x)/l/2 + 1; // outside x
//Transfinite Line{227,206,51,192,60,188,140,52,144,56,73,100,82,122,104,126} = n+1; // z
//
//Transfinite Surface{28,30,32,34,36,38,40,42,44}; // xy prallel (bottom)
//Transfinite Surface{88,110,132,154,176,66,242,220,198};// xy parallel (top)
//Transfinite Surface{75,83,105,127,53,61,141,149,229,207,197,189}; // xz prallel
//Transfinite Surface{87,109,131,57,101,123,65,167,145,241,211,193};// yz parallel
//
////surfaceVector[] = Extrude {0, 0, n} {
////  Surface{28,30,32,34,36,38,40,42,44};
////  Layers{1};
////  Recombine;
////};
//
//
//Physical Surface("T") = {88, 66, 242, 220, 176, 110, 132, 154, 198};
//Physical Surface("O") = {193, 211, 241};
//Physical Surface("N") = {229, 53, 75};
//Physical Surface("W") = {87, 109, 131};
//Physical Surface("B") = {28, 30, 32, 42, 44, 34, 36, 38, 40};
//Physical Surface("S") = {189, 149, 127};
//Physical Volume("CORE") = {6};
//Physical Volume("MATRIX") = {9, 8, 7, 5, 4, 3, 2, 1};
////Physical Volume("SOLID") = {6, 9, 8, 7, 5, 4, 3, 2, 1};
//Physical Point("NOT") = {75};
//Physical Point("NWT") = {27};
//Physical Point("SWT") = {56};
//Physical Point("SOT") = {68};
//Physical Point("NOB") = {1};
//Physical Point("NWB") = {2};
//Physical Point("SWB") = {3};
//Physical Point("SOB") = {4};
//Physical Line("NW") = {73};
//Physical Line("SW") = {126};
//Physical Line("SO") = {188};
//Physical Line("NO") = {227};
//Physical Line("NT") = {68, 46, 222};
//Physical Line("OT") = {225, 201, 180};
//Physical Line("ST") = {179, 136, 114};
//Physical Line("WT") = {115, 93, 71};
//Physical Line("NB") = {1, 7, 16};
//Physical Line("WB") = {17, 19, 21};
//Physical Line("SB") = {14, 13, 5};
//Physical Line("OB") = {23, 24, 26};
//+
Physical Line("left") = {17, 19, 21};
//+
Physical Line("right") = {26, 24, 23};
//+
Physical Line("down") = {14, 13, 5};
//+
Physical Line("up") = {16, 7, 1};
//+
Physical Surface("outer") = {28, 34, 36, 38, 40, 42, 32, 30};
//+
Physical Surface("inner") = {44};
//+
Transfinite Surface{28,30,32,34,36,38,40,42,44};
