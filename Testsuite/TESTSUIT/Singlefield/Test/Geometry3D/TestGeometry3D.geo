// create geometry

r = 1.0;
h = 0.25;//1.2;
Mesh.ElementOrder = 2;
//+
Point(1) = {0, 0, 0, h};
//+
Point(2) = {r, 0, 0, h};
//+
Point(3) = {0, r, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {1, 3};
//+
Ellipse(3) = {2, 1, 3, 3};
//+
Line Loop(1) = {1, 3, -2};
//+
Plane Surface(1) = {1};
//+
Extrude {{0, -1, 0}, {0, 0, 0}, Pi/2} {
  Surface{1}; 
}
//+
Physical Point("P_origin") = {1};
//+
Physical Line("L_x") = {1};
//+
//Physical Line("L_y") = {2};
//+
Physical Line("L_z") = {5};
//+
Physical Surface("S_8") = {13};
//+
Physical Surface("S_yz") = {15};
//+
Physical Surface("S_xz") = {10};
//+
Physical Surface("S_xy") = {1};
//+
Physical Volume("V_8") = {1};
