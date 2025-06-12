// Gmsh (http://geuz.org/gmsh) model  Couette flow according to p.  2, Sec 2.1
// in Techreport (Kumar2006)  Kumar, V.  & Huber, C.  Determination of Viscous
// Forces on Oscillating  Measuring Tubes, Part I :  Verification of Numerical
// Tools Endress+Hauser, 2006

W = 0.2;
H = 4;
Point(1) = {0, 0, 0};
Point(2) = {W, 0, 0};
Point(3) = {W, H, 0};
Point(4) = {0, H, 0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Transfinite Surface{6} = {};
Transfinite Line {2} = 70 Using Progression 1.1;
Transfinite Line {-4} = 70 Using Progression 1.1;
Transfinite Line {-1} = 3 Using Progression 1.0;
Transfinite Line {3} = 3 Using Progression 1.0;

Recombine Surface{6};

Physical Surface("channel") = {6};
  
