/* Trench heater micro region, finger  in center of multiple std cells
 * S7
 *
 * Author: Sebastian Eiser
 * Date: Jan. 2013
 * Use: This file assembles building blocks of the s7 trench poly heater of Mischa (cf. Exmon5/ Polyheater04)
 */

General.Terminal=1;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.Algorithm=8; // delquad mesher (better Blossom results)

Point(1) = {0, 0, 0, 0.4};
Point(2) = {0, 1, 0, 0.6};
Point(3) = {1, 1, 0, 0.3};
Point(4) = {1, 0, 0, 0.7};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};
Extrude {0, 0, 20} {
  Surface{6};
  Layers{6};
  Recombine;
}
Extrude {0, 0, 10} {
  Surface{28};
  Layers{6};
  Recombine;
}
Physical Surface("hot") = {50};
Physical Surface("cold") = {6};
Physical Surface("meas") = {28};
Physical Volume("vWire1") = {1};
Physical Volume("vWire2") = {2};
