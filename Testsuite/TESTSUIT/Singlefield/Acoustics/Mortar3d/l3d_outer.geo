wprop = 0.5;
wpml = 0.15;
numelemsI = 5;
numelemsPML = 3;

Point(1) = {0.0, 0.0, 0.0, 1.0};
Point(2) = {wprop, 0.0, 0.0, 1.0};
Point(3) = {wprop, wprop, 0.0, 1.0};
Point(4) = {0.0, wprop, 0.0, 1.0};
Point(5) = {wprop+wpml, 0.0, 0.0, 1.0};
Point(6) = {wprop+wpml, wprop, 0.0, 1.0};
Point(7) = {wprop+wpml, wprop+wpml, 0.0, 1.0};
Point(8) = {0.0, wprop+wpml, 0.0, 1.0};
Point(9) = {wprop, wprop+wpml, 0.0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {1, 4};
Line(5) = {4, 8};
Line(6) = {8, 9};
Line(7) = {3, 9};
Line(8) = {9, 7};
Line(9) = {3, 6};
Line(10) = {6, 7};
Line(11) = {2, 5};
Line(12) = {5, 6};
Line Loop(13) = {1, 2, -3, -4};
Plane Surface(14) = {13};
Line Loop(15) = {3, 7, -6, -5};
Plane Surface(16) = {15};
Line Loop(17) = {9, 10, -8, -7};
Plane Surface(18) = {17};
Line Loop(19) = {11, 12, -9, -2};
Plane Surface(20) = {19};
Extrude {0, 0, wprop} {
  Surface{14, 16, 18, 20};
}
Extrude {0, 0, wpml} {
  Surface{42, 108, 64, 86};
}
Transfinite Line {110, 22, 1, 112, 24, 3, 156, 46, 6, 58, 54, 76, 36, 32, 72, 27, 28, 94, 12, 2, 4, 25, 23, 89, 133, 111, 113} = numelemsI Using Progression 1;
Transfinite Line {115, 124, 168, 164, 68, 178, 134, 186, 66, 142, 155, 45, 120, 177, 67, 132, 116, 88, 138, 11, 9, 8, 7, 10, 5, 47, 157} = numelemsPML Using Progression 1;

Transfinite Surface "*";
Transfinite Volume{5} = {38, 39, 11, 10, 47, 43, 15, 19};
Transfinite Volume{1} = {1, 2, 3, 4, 10, 11, 15, 19};
Transfinite Volume{7} = {19, 15, 25, 29, 47, 43, 59, 63};
Transfinite Volume{8} = {15, 31, 35, 25, 43, 53, 69, 59};
Transfinite Volume{6} = {11, 37, 31, 15, 39, 49, 53, 43};
Transfinite Volume{4} = {2, 5, 6, 3, 11, 37, 31, 15};
Transfinite Volume{3} = {3, 6, 7, 9, 15, 31, 35, 25};
Transfinite Volume{2} = {4, 3, 9, 8, 19, 15, 25, 29};

Mesh.RecombineAll = 1;

Physical Volume("pml") = {2,3,4,5,6,7,8};

Physical Surface("outer_iface1") = {33};
Physical Surface("outer_iface2") = {37};
Physical Surface("outer_iface3") = {42};
