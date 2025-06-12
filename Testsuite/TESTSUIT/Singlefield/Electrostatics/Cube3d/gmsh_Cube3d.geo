Point(1) = {-1, -1, -1, 2.0};
Point(2) = {1, -1, -1, 2.0};
Point(3) = {1, 1, -1, 2.0};
Point(4) = {-1, 1, -1, 2.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {2, 3, 4, 1};
Plane Surface(6) = {5};
// Extrude {0, 0, 2} {
//   Surface{6};
// }

Extrude {0, 0, 2} {
  Surface{6}; Layers{2};
  Recombine;
}

Physical Surface("lower") = {27};
Physical Surface("upper") = {19};
Physical Volume("elec3d") = {1};
