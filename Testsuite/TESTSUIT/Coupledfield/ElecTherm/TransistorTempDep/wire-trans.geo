/* 
 * 
 */

Mesh.ScalingFactor = 1e-3;

elxy = 1;
elzWire = 5; // Layers in wire
elzCh = 1; // Layers in channel


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
Transfinite Surface {6};

Extrude {0, 0, .2}{ // 200 um Si
  Surface{6};
  Layers{elzWire};
  Recombine;
}
Extrude {0,0,0.02} { // 20 um 
  Surface{28};
  Layers{elzCh};
  Recombine;
}
Extrude {0, 0, 0.02} { // say, metal contact
  Surface{50};
  Layers{elzWire};
  Recombine;
}
Physical Surface("hot") = {72};
Physical Surface("cold") = {6};
Physical Surface("Usource") = {28};
Physical Surface("Udrain") = {50};
Physical Volume("vWire1") = {1};
Physical Volume("vChannel") = {2};
Physical Volume("vWire2") = {3};

Point(299) = {2, 2, 0.0, 0.4};
Point(300) = {3, 2, 0.0, 0.4};
Point(301) = {2, 3, 0.0, 0.4};
Point(302) = {3, 3, 0.0, 0.4};
Line(73) = {299, 301};
Line(74) = {301, 302};
Line(75) = {302, 300};
Line(76) = {299, 300};
Transfinite Line {73,74,75,76} = 2 Using Progression 1;
Line Loop(77) = {75, -76, 73, 74};
Plane Surface(78) = {77};
Extrude {0, 0, .005} {
  Surface{78};
  Layers{1};
  Recombine;
}
Extrude {0, 0, .005} {
  Surface{100};
  Layers{1};
  Recombine;
}
Physical Surface("Ugate") =  {100};
Physical Volume("vGate") = {4,5};
