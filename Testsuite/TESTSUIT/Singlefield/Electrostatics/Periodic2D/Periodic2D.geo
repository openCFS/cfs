// Gmsh project created on Thu Nov 28 09:41:57 2013
c = 0.5e-0;
a = 1.0;
k1 = 0.5;
k2 = 1.0;
Point(1) = {-a, -a, 0.0, k1*c*a};
Point(2) = {a, -a, 0.0, k2*c*a};
Point(3) = {a, a, 0.0, k1*c*a};
Point(4) = {-a, a, 0.0, k2*c*a};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(7) = {1, 2, 3, 4};
Plane Surface(8) = {7};

Physical Line("left") = {1};
Physical Line("right") = {3};
Physical Line("bottom") = {2};
Physical Line("top") = {4};
Physical Surface("mat") = {8};
