// Gmsh  (http://geuz.org/gmsh)  model  for  Stokes  flow  around  a  circular
// obstacle according  to p.  63 ff. in  PhD thesis of Gerhard  Link: A Finite
// Element Scheme  for Fluid–Solid–Acoustics Interactions  and its Application
// to Human Phonation.

R = 8e-4;
R2 = R/3;
L_in = 0.01;
L_out = 3*L_in;
H = L_in;

esize1 = 0.0005;
esize2 = H/4;
div_circum = 5;
div_radial = 3;

Point(1) = {0, 0, 0, esize1};

Point(2) = {R, 0, 0, esize1};
Point(3) = {0, R, 0, esize1};
Point(4) = {-R, 0, 0, esize1};
Point(5) = {0, -R, 0, esize1};

Point(6) = {R + R2, 0, 0, esize2};
Point(7) = {0, R + R2, 0, esize2};
Point(8) = {-R - R2, 0, 0, esize2};
Point(9) = {0, -R - R2, 0, esize2};

Point(10) = {-L_in, -H, 0, esize2};
Point(11) = {0, -H, 0, esize2};
Point(12) = {L_in, -H, 0, esize2};
Point(13) = {L_out, -H, 0, esize1};

Point(14) = {-L_in, H, 0, esize2};
Point(15) = {L_in, H, 0, esize2};
Point(16) = {0, H, 0, esize2};
Point(17) = {L_out, H, 0, esize1};

Point(18) = {-L_in, 0, 0, esize2};
Point(19) = {L_in, 0, 0, esize2};
Point(20) = {L_out, 0, 0, esize1};


Circle(1) = {3, 1, 4};
Circle(2) = {4, 1, 5};
Circle(3) = {5, 1, 2};
Circle(4) = {2, 1, 3};
Circle(5) = {7, 1, 8};
Circle(6) = {8, 1, 9};
Circle(7) = {9, 1, 6};
Circle(8) = {6, 1, 7};
Line(9) = {7, 16};
Line(10) = {8, 18};
Line(11) = {9, 11};
Line(12) = {6, 19};
Line(13) = {3, 7};
Line(14) = {4, 8};
Line(15) = {5, 9};
Line(16) = {2, 6};
Line(17) = {19, 20};
Line(18) = {10, 11};
Line(19) = {11, 12};
Line(20) = {12, 13};
Line(21) = {13, 20};
Line(22) = {20, 17};
Line(23) = {17, 15};
Line(24) = {15, 16};
Line(25) = {16, 14};
Line(26) = {14, 18};
Line(27) = {18, 10};
Line(28) = {19, 12};
Line(29) = {19, 15};
Line Loop(30) = {4, 13, -8, -16};
Plane Surface(31) = {30};
Line Loop(32) = {1, 14, -5, -13};
Plane Surface(33) = {32};
Line Loop(34) = {2, 15, -6, -14};
Plane Surface(35) = {34};
Line Loop(36) = {3, 16, -7, -15};
Plane Surface(37) = {36};
Line Loop(38) = {12, 29, 24, -9, -8};
Plane Surface(39) = {38};
Line Loop(40) = {5, 10, -26, -25, -9};
Plane Surface(41) = {40};
Line Loop(42) = {10, 27, 18, -11, -6};
Plane Surface(43) = {42};
Line Loop(44) = {11, 19, -28, -12, -7};
Plane Surface(45) = {44};
Line Loop(46) = {20, 21, -17, 28};
Plane Surface(47) = {46};
Line Loop(48) = {17, 22, 23, -29};
Plane Surface(49) = {48};

Transfinite Line {1} = div_circum Using Progression 1.0;
Transfinite Line {2} = div_circum Using Progression 1.0;
Transfinite Line {3} = div_circum Using Progression 1.0;
Transfinite Line {4} = div_circum Using Progression 1.0;

Transfinite Line {5} = div_circum Using Progression 1.0;
Transfinite Line {6} = div_circum Using Progression 1.0;
Transfinite Line {7} = div_circum Using Progression 1.0;
Transfinite Line {8} = div_circum Using Progression 1.0;

Transfinite Line {13} = div_radial Using Progression 1.0;
Transfinite Line {14} = div_radial Using Progression 1.0;
Transfinite Line {15} = div_radial Using Progression 1.0;
Transfinite Line {16} = div_radial Using Progression 1.0;

Transfinite Surface{31} = {};
Transfinite Surface{33} = {};
Transfinite Surface{35} = {};
Transfinite Surface{37} = {};

Recombine Surface{31};
Recombine Surface{33};
Recombine Surface{35};
Recombine Surface{37};

Transfinite Line {17} = div_circum*2 Using Progression 1.0;
Transfinite Line {20} = div_circum*2 Using Progression 1.0;
Transfinite Line {23} = div_circum*2 Using Progression 1.0;

Transfinite Line {28} = div_circum Using Progression 1.0;
Transfinite Line {29} = div_circum Using Progression 1.0;
Transfinite Line {21} = div_circum Using Progression 1.0;
Transfinite Line {22} = div_circum Using Progression 1.0;

Transfinite Surface{47} = {};
Transfinite Surface{49} = {};

Recombine Surface{47};
Recombine Surface{49};

Physical Surface("channel") = {41, 43, 39, 45, 49, 47, 31, 33, 35, 37};
Physical Line("wall") = {18, 19, 20, 23, 24, 25};
Physical Line("cyl") = {1, 2, 3, 4};
Physical Line("inflow") = {26, 27};
Physical Line("outflow") = {21, 22};

  
