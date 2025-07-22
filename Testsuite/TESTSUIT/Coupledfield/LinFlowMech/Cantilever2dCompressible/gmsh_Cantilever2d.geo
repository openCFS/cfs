// Channel model for  Stokes flow according to M.Sc.   thesis of Gerhard Link:
// Numerical Simulation of FSI. pages 41 ff.

// Thin channel
L_C = 12e-2;
H_C = 3e-3;
D = 10e-2;

Point(1) = {0, 0, 0};
Point(2) = {L_C, 0, 0};
Point(3) = {L_C+D, 0, 0};

Point(4) = {0, D, 0};
Point(5) = {L_C, D, 0};
Point(6) = {L_C+D, D, 0};

Point(7) = {0, D+H_C, 0};
Point(8) = {L_C, D+H_C, 0};
Point(9) = {L_C+D, D+H_C, 0};

Point(10) = {0, D+H_C+D, 0};
Point(11) = {L_C, D+H_C+D, 0};
Point(12) = {L_C+D, D+H_C+D, 0};

// Line(1) = {1, 2};
// Line(2) = {2, 3};
// Line(3) = {3, 4};
// Line(4) = {4, 1};
// Line Loop(6) = {4, 1, 2, 3};
// Plane Surface(6) = {6};
// Transfinite Surface{6} = {};

// Transfinite Line {2} = 2 Using Progression 1.0;
// Transfinite Line {-4} = 2 Using Progression 1.0;
// Transfinite Line {-1} = 3 Using Progression 1.0;
// Transfinite Line {3} = 3 Using Progression 1.0;

// In den  angegebenen Surfaces sollen  die strukturierten Dreiecke  wieder zu
// Vierecken zusammengefasst werden.
// Recombine Surface{6};

// Physical Surface("channel") = {6};
// Physical Line("inlet") = {4};
// Physical Line("outlet") = {2};
// Physical Line("noslip") = {1,3};

  
  
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {7, 8};
Line(6) = {8, 9};
Line(7) = {10, 11};
Line(8) = {11, 12};
Line(9) = {1, 4};
Line(10) = {4, 7};
Line(11) = {7, 10};
Line(12) = {2, 5};
Line(13) = {5, 8};
Line(14) = {8, 11};
Line(15) = {3, 6};
Line(16) = {6, 9};
Line(17) = {9, 12};
Line Loop(18) = {1, 12, -3, -9};
Plane Surface(19) = {18};
Transfinite Surface{19} = {};
Line Loop(20) = {2, 15, -4, -12};
Plane Surface(21) = {20};
Transfinite Surface{21} = {};
Line Loop(22) = {3, 13, -5, -10};
Plane Surface(23) = {22};
Transfinite Surface{23} = {};
Line Loop(24) = {4, 16, -6, -13};
Plane Surface(25) = {24};
Transfinite Surface{25} = {};
Line Loop(26) = {5, 14, -7, -11};
Plane Surface(27) = {26};
Transfinite Surface{27} = {};
Line Loop(28) = {6, 17, -8, -14};
Plane Surface(29) = {28};
Transfinite Surface{29} = {};

Transfinite Line {1} = 9 Using Progression 1.0;
Transfinite Line {3} = 9 Using Progression 1.0;
Transfinite Line {5} = 9 Using Progression 1.0;
Transfinite Line {7} = 9 Using Progression 1.0;

Transfinite Line {2} = 7 Using Progression 1.0;
Transfinite Line {4} = 7 Using Progression 1.0;
Transfinite Line {6} = 7 Using Progression 1.0;
Transfinite Line {8} = 7 Using Progression 1.0;

Transfinite Line {9} = 8 Using Progression 1.0;
Transfinite Line {12} = 8 Using Progression 1.0;
Transfinite Line {15} = 8 Using Progression 1.0;

Transfinite Line {10} = 2 Using Progression 1.0;
Transfinite Line {13} = 2 Using Progression 1.0;
Transfinite Line {16} = 2 Using Progression 1.0;

Transfinite Line {11} = 8 Using Progression 1.0;
Transfinite Line {14} = 8 Using Progression 1.0;
Transfinite Line {17} = 8 Using Progression 1.0;

Recombine Surface{19, 21, 23, 25, 27, 29};

Physical Surface("cantilever") = {23};
Physical Surface("water") = {19, 21, 25, 27, 29};
Physical Line("wall") = {9, 11};
Physical Line("free") = {1, 2, 15, 16, 17, 8, 7};
Physical Line("fix") = {10};
Physical Line("cplsurf") = {3,13,5};
