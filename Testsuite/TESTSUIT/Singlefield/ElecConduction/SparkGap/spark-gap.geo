/* Trench heater micro region, finger  in center of multiple std cells
 * S7
 *
 * Author: Sebastian Eiser
 * Date: Jan. 2013
 * Use: This file assembles building blocks of the s7 trench poly heater of Mischa (cf. Exmon5/ Polyheater04)
 */

elxy = 2;
elzWire = 4; // Layers in wire
elzCh = 3; // Layers in channel


General.Terminal=1;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.Algorithm=8; // delquad mesher (better Blossom results)

Point(1) = {0, 0, 0, 9.4};
Point(2) = {0, 1, 0, 9.6};
Point(3) = {1, 1, 0, 9.3};
Point(4) = {1, 0, 0, 9.7};
Line(1) = {2, 3};
Line(2) = {3, 4};
Line(3) = {4, 1};
Line(4) = {1, 2};
Transfinite Line {2, 1, 4, 3} = elxy Using Progression 1;
//Transfinite Line {1,2,3,4} = 2;
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

Extrude {0, 0, 5}{
  Surface{6};
  Layers{elzWire};
  Recombine;
}
Extrude {0,0,3} {
  Surface{28};
  Layers{elzCh};
  Recombine;
}
Extrude {0, 0, 5} {
  Surface{50};
  Layers{elzWire};
  Recombine;
}
Physical Surface("hot") = {72};
Physical Surface("cold") = {6};
Physical Surface("anode") = {28};
Physical Surface("cathode") = {50};
Physical Volume("vWire1") = {1};
Physical Volume("chan") = {2};
Physical Volume("vWire2") = {3};

