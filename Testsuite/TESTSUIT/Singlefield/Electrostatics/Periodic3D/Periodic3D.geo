// Gmsh project created on Thu Nov 28 09:41:57 2013
c = 1.0e-0;
a = 0.1;
k1 = 1;
k2 = 1;
Point(1) = {-k1*a, -k2*a, 0, k1*c*a};
Point(2) = {k1*a, -k2*a, 0, k1*c*a};
Point(3) = {k1*a, k2*a, 0, k1*c*a};
Point(4) = {-k1*a, k2*a, 0, k1*c*a};

d = -2*k1*a;
Extrude {0, 0, d} {
  Point{4, 3, 1, 2};
}

Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {8, 6};
Line(10) = {5, 7};
Line(11) = {8, 7};
Line(12) = {6, 5};
Line Loop(13) = {5, 6, 7, 8};
Plane Surface(14) = {13};
Line Loop(15) = {9, 12, 10, -11};
Plane Surface(16) = {-15};
Line Loop(17) = {4, 9, -2, -6};
Plane Surface(18) = {17};
Line Loop(19) = {1, 10, -3, -8};
Plane Surface(20) = {19};
Line Loop(21) = {2, 12, -1, -7};
Plane Surface(22) = {21};
Line Loop(23) = {3, -11, -4, -5};
Plane Surface(24) = {23};
Surface Loop(25) = {16, 24, 20, 22, 18, 14};
Volume(26) = {25};

Physical Volume("mat") = {26};
Physical Surface("left") = {20};
Physical Surface("right") = {18};
Physical Surface("top") = {22};
Physical Surface("bottom") = {24};
//Physical Surface("front") = {14};
//Physical Surface("rear") = {16};
