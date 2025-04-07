// Gmsh (http://geuz.org/gmsh) model for Couette  flow according to p.  2, Sec
// 2.1  in Techreport  (Kumar2006) Kumar,  V.  &  Huber, C.   Determination of
// Viscous Forces  on Oscillating  Measuring Tubes, Part  I :  Verification of
// Numerical Tools Endress+Hauser, 2006

R = 8.8e-3;
R_in = R / 10;

Point(1) = {0, 0, 0};
Point(2) = {R_in, 0, 0};
Point(3) = {R, 0, 0};
Point(4) = {0, R_in, 0};
Point(5) = {0, R, 0};
Point(6) = {-R_in, 0, 0};
Point(7) = {-R, 0, 0};
Point(8) = {0, -R_in, 0};
Point(9) = {0, -R, 0};

Line(1) = {2, 3};
Line(2) = {4, 5};
Circle(3) = {4, 1, 2};
Circle(4) = {5, 1, 3};

Line(5) = {6, 7};
Line(6) = {8, 9};

Line Loop(6) = {1, -4, -2, 3};
Plane Surface(6) = {6};
Circle(7) = {4, 1, 6};
Circle(8) = {6, 1, 8};
Circle(9) = {8, 1, 2};
Circle(10) = {7, 1, 5};
Circle(11) = {7, 1, 9};
Circle(12) = {9, 1, 3};
Line Loop(13) = {10, -2, 7, 5};
Plane Surface(14) = {13};
Line Loop(15) = {5, 11, -6, -8};
Plane Surface(16) = {15};
Line Loop(17) = {6, 12, -1, -9};
Plane Surface(18) = {17};

Transfinite Surface{6} = {};
Transfinite Surface{14} = {};
Transfinite Surface{16} = {};
Transfinite Surface{18} = {};

Transfinite Line {-1} = 70 Using Progression 1.1;
Transfinite Line {-2} = 70 Using Progression 1.1;
Transfinite Line {-5} = 70 Using Progression 1.1;
Transfinite Line {-6} = 70 Using Progression 1.1;
Transfinite Line {3} = 30 Using Progression 1.0;
Transfinite Line {4} = 30 Using Progression 1.0;
Transfinite Line {7} = 30 Using Progression 1.0;
Transfinite Line {8} = 30 Using Progression 1.0;
Transfinite Line {9} = 30 Using Progression 1.0;
Transfinite Line {10} = 30 Using Progression 1.0;
Transfinite Line {11} = 30 Using Progression 1.0;
Transfinite Line {12} = 30 Using Progression 1.0;

Recombine Surface{6};
Recombine Surface{14};
Recombine Surface{16};
Recombine Surface{18};

Physical Surface("pipe") = {6, 14, 16, 18};

//Physical Line("right") = {1};
//Physical Line("left") = {2};
  
Physical Line("inner") = {3, 7, 8, 9};
Physical Line("outer") = {4, 10, 11, 12};
