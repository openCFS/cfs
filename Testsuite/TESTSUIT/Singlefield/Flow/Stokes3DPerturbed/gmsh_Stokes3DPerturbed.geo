h=0.1;
phi=3.14159/6;

Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {Cos(phi), Sin(phi), 0, h};
Point(4) = {1.3*Cos(phi), 1.3*Sin(phi), 0, h};
Point(5) = {1.3, 0, 0, h};
Line(1) = {2, 5};
Circle(2) = {5, 1, 4};
Line(3) = {4, 3};
Circle(4) = {2, 1, 3};
Line Loop(5) = {2, 3, -4, 1};
Ruled Surface(6) = {5};
Extrude {0, 0, 0.3} {
  Surface{6};
}

Physical Volume("channel") = {1};
Physical Surface("inflow") = {27};
Physical Surface("outflow") = {19};
Physical Surface("wall") = {23, 6, 15, 28};
