radius = 1;
box = 1;
box2 = 2;

//create outer rectangle
Point (7) = {-box2,-box2,0};
Point (8) = {box2,-box2,0};
Point (9) = {box2,box2,0};
Point (10) = {-box2,box2,0};

Point (11) = {-box,-box,0};
Point (12) = {box,-box,0};
Point (13) = {box,box,0};
Point (14) = {-box,box,0};

Point (15) = {-box,-box2,0};
Point (16) = {-box2,-box,0};

Point (17) = {box,-box2,0};
Point (18) = {box2,-box,0};

Point (19) = {box,box2,0};
Point (20) = {box2,box,0};

Point (21) = {-box,box2,0};
Point (22) = {-box2,box,0};

Line(6) = {13, 20};
Line(7) = {13, 19};
Line(8) = {14, 21};
Line(9) = {14, 22};
Line(10) = {11, 16};
Line(11) = {11, 15};
Line(12) = {12, 18};
Line(13) = {12, 17};
Line(14) = {17, 8};
Line(15) = {8, 18};
Line(16) = {18, 20};
Line(17) = {20, 9};
Line(18) = {9, 19};
Line(19) = {19, 21};
Line(20) = {21, 10};
Line(21) = {10, 22};
Line(22) = {22, 16};
Line(23) = {16, 7};
Line(24) = {7, 15};
Line(25) = {15, 17};
Line(27) = {13, 12};
Line(28) = {12, 11};
Line(29) = {11, 14};
Line(30) = {14, 13};

Coherence;

Line Loop(31) = {12, 16, -6, 27};
Plane Surface(32) = {31};
Line Loop(33) = {6, 17, 18, -7};
Plane Surface(34) = {33};
Line Loop(35) = {30, 7, 19, -8};
Plane Surface(36) = {35};
Line Loop(37) = {9, -21, -20, -8};
Plane Surface(38) = {37};
Line Loop(39) = {10, -22, -9, -29};
Plane Surface(40) = {39};
Line Loop(41) = {24, -11, 10, 23};
Plane Surface(42) = {41};
Line Loop(43) = {25, -13, 28, 11};
Plane Surface(44) = {43};
Line Loop(45) = {14, 15, -12, 13};
Plane Surface(46) = {45};

Physical Line ("outerR") = {27};
Physical Line ("outerL") = {29};
Physical Line ("outerT") = {30};
Physical Line ("outerB") = {28};
Physical Line ("abc") = {16,17,18,19,20,21,22,23,24,25,14,15};

Physical Surface ("outer") = {44,42,-40,-38,36,34,32,46};

Transfinite Line {14, 13, 15, 12, 6, 7, 17, 18, 8, 20, 21, 9, 24, 11, 10, 23} = 5 Using Progression 1;
Transfinite Line {27, 16, 29, 22, 28, 25, 30, 19} = 5 Using Progression 1;

Transfinite Surface '*';
Mesh.RecombineAll = 1;

//Mesh 2;

