cl1 = 8;
Point(1) = {7.5, 6, 0, cl1};
Point(2) = {0, 0, 0, cl1};
Point(3) = {15, 0, 0, cl1};
Point(4) = {15, 12, 0, cl1};
Point(5) = {0, 12, 0, cl1};
Point(6) = {40, 0, 0, cl1};
Point(7) = {40, 12, 0, cl1};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};
Circle(5) = {6, 1, 6};
Line Loop(8) = {2, 3, 4, 1, -5};
Plane Surface(8) = {8};

Point(20) = {3, 0, 0, cl1};
Point(21) = {4, 0, 0, cl1};
Point(22) = {7, 0, 0, cl1};
Point(23) = {8, 0, 0, cl1};
Point(24) = {11, 0, 0, cl1};
Point(25) = {12, 0, 0, cl1};

Line(30) = {2, 20};
Line(31) = {20, 21};
Line(32) = {21, 22};
Line(33) = {22, 23};
Line(34) = {23, 24};
Line(35) = {24, 25};
Line(36) = {25, 3};
Extrude {0, 3, 0} {
  Line{30, 31, 32, 33, 34, 35, 36};
}
Extrude {0, 1, 0} {
  Line{37, 41, 45, 49, 53, 57, 61};
}
Extrude {0, 5, 0} {
  Line{65, 69, 73, 77, 81, 85, 89};
}
Extrude {0, 3, 0} {
  Line{93, 97, 101, 105, 109, 113, 117};
}

Transfinite Line {31,33,35,41,49,57,69,77,85,97,105,113,125,133,141,66,67,71,75,79,83,87,91} = 2 Using Progression 1;

Transfinite Line {30,32,34,36,37,45,53,61,65,73,81,89,93,101,109,117,121,129,137,145,38,39,43,47,51,55,59,63,122,123,127,131,135,139,143,147} = 4 Using Progression 1;

Transfinite Line {94,95,99,103,107,111,115,119} = 6 Using Progression 1;

Transfinite Surface {40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148} = {};

Mesh.RecombineAll = 1;

Physical Line("lower_outer") = {1};
Physical Line("right_outer") = {2};
Physical Line("upper_outer") = {3};
Physical Line("left_outer") = {4};

Physical Surface("outer") = { 8 };

  
Physical Surface("coil") = {100, 72, 76, 80, 108, 84, 88, 116};
Physical Surface("heat_inner") = {40, 44, 48, 52, 56, 60, 64, 68, 92, 96, 104, 112, 120, 124, 128, 132, 136, 140, 144, 148};

Physical Line("left_inner") = {38,66,94,122};
Physical Line("right_inner") = {63,91,119,147};
Physical Line("lower_inner") = {30,31,32,33,34,35,36};
Physical Line("upper_inner") = {121,125,129,133,137,141,145};

// cfstool -m convert heat_outer.msh  out.h5
// cfstool --formatArgs '{"input":{"hdf5":{"readEntities":"outer left_outer right_outer lower_outer upper_outer","linearizeEntities":"__none__"}}}' -m convert results_hdf5/out.h5  outer.h5
// cfstool --formatArgs '{"input":{"hdf5":{"readEntities":"coil heat_inner left_inner right_inner lower_inner upper_inner","linearizeEntities":"__all__"}}}' -m convert results_hdf5/out.h5  inner.h5