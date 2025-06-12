Mesh.SecondOrderLinear = 2;
Mesh.SecondOrderIncomplete = 1;
//Mesh.SubdivisionAlgorithm = 0;
//Mesh.RemeshAlgorithm = 0;
//Mesh.RemeshParametrization = 1;
//Mesh.RecombinationAlgorithm = 1;
Mesh.RecombineAll = 1;
//Mesh.Algorithm = 8;
//Mesh.Algorithm3D = 6;
Mesh.CharacteristicLengthFactor = 0.4;
Mesh.ElementOrder = 2;
//Mesh.Optimize = 1;
//Mesh.CharacteristicLengthFromCurvature = 1;
Mesh.Smoothing = 10;
// Mesh 3;
// Coherence Mesh;
Geometry.CopyMeshingMethod=1;

radius = 0.1*0.5;
lengthTube = 0.8;
volLayers = 160;
radiusCut = 2;

Point(2) = {0, 0, 0};
Point(1) = {radius, 0, 0};

Line(10) = {2, 1};
Transfinite Line {10} = radiusCut;
Extrude {0, lengthTube, 0} {
  Line{10};
  Layers{volLayers};
  Recombine;
}

Physical Line("excite") = {10};
//Physical Surface("surfInflow") = {6};
Physical Line("surfOutflow") = {11};
Physical Surface("acouVol") = {14};
