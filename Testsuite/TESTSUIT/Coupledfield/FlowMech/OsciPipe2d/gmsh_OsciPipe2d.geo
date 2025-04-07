// kate: indent-width 2; mixedindent on; indent-mode cstyle; line-numbers on;
// kate: syntax c; replace-tabs on;
// kate: remove-trailing-space on; bracket-highlight-color #ff00ff;

nu = 1e-3;
f = 725;
pi = 4*Atan(1);
omega = 2*pi*f;
// Grenzschicht
G = 6*Sqrt(nu / omega);


Printf("Kin. Viskosität nu: %f", nu);
Printf("Frequenz f: %f", f);
Printf("pi: %f", pi);
Printf("omega: %f", omega);
Printf("Grenzschichtdicke G: %f", G);

// Radius
r = 8.8e-3;

// Radius für inneres Quadrat
R_q = r / 5;

// Rohraußenradius
R_a = r+12e-4;

// Rohrinnenradius
// R_i = R_a - W;
R_i = r;

// Grenzschichtradius
R_g = R_i - G;

N=Ceil((R_g-R_q) / G);
Printf("N: %f", N);

// Verfeinerung
feinheit = 1;

N_radial = (4*feinheit)+1;
N_grenz = (4*feinheit)+1;
N_rohr = 1;
prog = 1.0 + 0.75^feinheit;
// prog = 1;

Printf("Progression prog: %f", prog);

Point(1) = {0, 0, 0};

Point(2) = {-R_q, 0, 0};
Point(3) = {0, -R_q, 0};
Point(4) = {R_q, 0, 0};
Point(5) = {0, R_q, 0};

Point(6) = {-R_g, 0, 0};
Point(7) = {0, -R_g, 0};
Point(8) = {R_g, 0, 0};
Point(9) = {0, R_g, 0};

Point(10) = {-R_i, 0, 0};
Point(11) = {0, -R_i, 0};
Point(12) = {R_i, 0, 0};
Point(13) = {0, R_i, 0};

Point(14) = {-R_a, 0, 0};
Point(15) = {0, -R_a, 0};
Point(16) = {R_a, 0, 0};
Point(17) = {0, R_a, 0};

Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 2};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 10};
Circle(1009) = {14, 1, 15};
Circle(1010) = {15, 1, 16};
Circle(1011) = {16, 1, 17};
Circle(1012) = {17, 1, 14};

Line(13) = {10, 6};
Line(14) = {6, 2};
Line(15) = {2, 1};
Line(16) = {11, 7};
Line(17) = {7, 3};
Line(18) = {3, 1};
Line(19) = {12, 8};
Line(20) = {8, 4};
Line(21) = {4, 1};
Line(22) = {13, 9};
Line(23) = {9, 5};
Line(24) = {5, 1};

Line(1013) = {14, 10};
Line(1014) = {15, 11};
Line(1015) = {16, 12};
Line(1016) = {17, 13};

Line Loop(25) = {1, 2, 3, 4};
Plane Surface(26) = {25};
// Line Loop(27) = {18, -21, -2};
// Plane Surface(28) = {27};
// Line Loop(29) = {21, -24, -3};
// Plane Surface(30) = {29};
// Line Loop(31) = {24, -15, -4};
// Plane Surface(32) = {31};

Transfinite Surface{26} = {};
Recombine Surface{26};
// Recombine Surface{28};
// Recombine Surface{30};
// Recombine Surface{32};

Line Loop(33) = {14, 1, -17, -5};
Ruled Surface(34) = {33};
Line Loop(35) = {17, 2, -20, -6};
Ruled Surface(36) = {35};
Line Loop(37) = {20, 3, -23, -7};
Ruled Surface(38) = {37};
Line Loop(39) = {23, 4, -14, -8};
Ruled Surface(40) = {39};

Transfinite Surface{34} = {};
Transfinite Surface{36} = {};
Transfinite Surface{38} = {};
Transfinite Surface{40} = {};

Transfinite Line {14} = (N*feinheit)+1 Using Progression 1.0;
Transfinite Line {1} = N_radial Using Progression 1.0;
Transfinite Line {17} = (N*feinheit)+1 Using Progression 1.0;
Transfinite Line {5} = N_radial Using Progression 1.0;
Recombine Surface{34};

Transfinite Line {2} = N_radial Using Progression 1.0;
Transfinite Line {20} = (N*feinheit)+1 Using Progression 1.0;
Transfinite Line {6} = N_radial Using Progression 1.0;
Recombine Surface{36};

Transfinite Line {3} = N_radial Using Progression 1.0;
Transfinite Line {23} = (N*feinheit)+1 Using Progression 1.0;
Transfinite Line {7} = N_radial Using Progression 1.0;
Recombine Surface{38};

Transfinite Line {4} = N_radial Using Progression 1.0;
Transfinite Line {8} = N_radial Using Progression 1.0;
Recombine Surface{40};

Line Loop(41) = {13, 5, -16, -9};
Ruled Surface(42) = {41};
Line Loop(43) = {16, 6, -19, -10};
Ruled Surface(44) = {43};
Line Loop(45) = {19, 7, -22, -11};
Ruled Surface(46) = {45};
Line Loop(47) = {22, 8, -13, -12};
Ruled Surface(48) = {47};

Transfinite Surface{42} = {};
Transfinite Surface{44} = {};
Transfinite Surface{46} = {};
Transfinite Surface{48} = {};

Line Loop(1041) = {1013, 9, -1014, -1009};
Ruled Surface(1042) = {1041};
Line Loop(1043) = {1014, 10, -1015, -1010};
Ruled Surface(1044) = {1043};
Line Loop(1045) = {1015, 11, -1016, -1011};
Ruled Surface(1046) = {1045};
Line Loop(1047) = {1016, 12, -1013, -1012};
Ruled Surface(1048) = {1047};

Transfinite Surface{42} = {};
Transfinite Surface{44} = {};
Transfinite Surface{46} = {};
Transfinite Surface{48} = {};

Transfinite Surface{1042} = {};
Transfinite Surface{1044} = {};
Transfinite Surface{1046} = {};
Transfinite Surface{1048} = {};

Transfinite Line {13} = N_grenz Using Progression prog;
Transfinite Line {16} = N_grenz Using Progression prog;
Transfinite Line {9} = N_radial Using Progression 1.0;
Recombine Surface{42};

Transfinite Line {19} = N_grenz Using Progression prog;
Transfinite Line {10} = N_radial Using Progression 1.0;
Recombine Surface{44};

Transfinite Line {22} = N_grenz Using Progression prog;
Transfinite Line {11} = N_radial Using Progression 1.0;
Recombine Surface{46};

Transfinite Line {12} = N_radial Using Progression 1.0;
Recombine Surface{48};



Transfinite Line {1013} = N_rohr Using Progression 1.0;
Transfinite Line {1014} = N_rohr Using Progression 1.0;
Transfinite Line {1015} = N_rohr Using Progression 1.0;
Transfinite Line {1016} = N_rohr Using Progression 1.0;

Transfinite Line {1009} = N_radial Using Progression 1.0;
Transfinite Line {1010} = N_radial Using Progression 1.0;
Transfinite Line {1011} = N_radial Using Progression 1.0;
Transfinite Line {1012} = N_radial Using Progression 1.0;

Recombine Surface{1042};
Recombine Surface{1044};
Recombine Surface{1046};
Recombine Surface{1048};


Physical Surface("water") = {26,34,36,38,40,42,44,46,48};
Color Grey50 { Surface {26}; }
Color Grey85 { Surface {34,36,38,40}; }
Color Yellow { Surface {42,44,46,48}; }

Physical Line("pipe_inner") = {9,10,11,12};
// Color Green { Line {9,10,11,12} };

Physical Surface("pipe") = {1042,1044,1046,1048};
Color Green { Surface {1042,1044,1046,1048}; }

Physical Line("pipe_outer") = {1009,1010,1011,1012};
// Color Green { Line {9,10,11,12} };

