overall = 10;
pmlW = 1;
plateW = 0.1;
plateH = 3.5;
platePosX = 2;
platePosY = 2.25;

elemSize = 0.1;


//this is a bottom up modeling approach
//for structured mesh
//points->lines->surfaces->mesh

//fist we define all necesary points for a structured mesh just for one 
//layer at leftmost domain border

//define basepoint for plate
basePoint = -overall/2+pmlW+platePosY;

Point(1) = {-overall/2, -overall/2       , 0, 1.0};
Point(2) = {-overall/2, -overall/2+pmlW  , 0, 1.0};
Point(3) = {-overall/2, basePoint        , 0, 1.0};
Point(4) = {-overall/2, basePoint+plateH , 0, 1.0};
Point(5) = {-overall/2, overall/2-pmlW   , 0, 1.0};
Point(6) = {-overall/2, overall/2        , 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};

//using this way of extrusion will lead to misoriented elements
//CFS++ will correct this lateron
Extrude{ pmlW, 0, 0 }{ Line{1,2,3,4,5}; }
Extrude{ platePosX, 0, 0 }{ Line{6,10,14,18,22}; }
Extrude{ plateW, 0, 0 }{ Line{26,30,34,38,42}; }
Extrude{ overall-2*pmlW-platePosX-plateW, 0, 0 }{ Line{46,50,54,58,62}; }
Extrude{ pmlW, 0, 0 }{ Line{66,70,74,78,82}; }

//combine double entities
Coherence;

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// GENERATE REGIONS TO ADRESS LATER IN CFS++
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//give a name and assign the surface numbers
Physical Surface("prop") = {33,53,73,38,77,41,61,81,37};
Physical Surface("plate") = {57};
Physical Surface("pml") = {9,29,49,69,89,93,97,101,105,85,65,45,25,21,17,13};


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// MESHING
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//=====================================================================
// Structured Mesh generation
//=====================================================================
//order gmsh to to structured meshing
Transfinite Surface '*';
//recombine structured triangles to obtain quad mesh
Recombine Surface '*';

//assign line subdivisions, size on points is to be ignored
//X-Lines
Transfinite Line {7,8,12,16,20,24,87,88,92,96,100,104} = ((pmlW)/elemSize) Using Progression 1;
Transfinite Line {27,28,32,36,40,44} = ((platePosX)/elemSize) Using Progression 1;
Transfinite Line {47,48,52,56,60,64} = ((plateW)/elemSize) Using Progression 1;
Transfinite Line {67,68,72,76,80,84} = ((overall-2*pmlW-platePosX-plateW)/elemSize) Using Progression 1;
//Y-Lines
Transfinite Line {1,6,26,46,66,86,5,22,42,62,82,102} = ((pmlW)/elemSize) Using Progression 1;
Transfinite Line {2,10,30,50,70,90} = ((platePosY)/elemSize) Using Progression 1;
Transfinite Line {3,14,34,54,74,94} = ((plateH)/elemSize) Using Progression 1;
Transfinite Line {4,18,38,58,78,98} = ((overall-2*pmlW-platePosY-plateH)/elemSize) Using Progression 1;

////======================================================================
//// Unstructured quadrilateral mesh
////======================================================================
////NOTE: element sizes defined while generating points
//
////recombine structured triangles to obtain quad mesh
//Recombine Surface '*';
//
////select meshing algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad)
//Mesh.Algorithm = 1;
//
////Sizing
//Mesh.CharacteristicLengthFromPoints = 0;
//Mesh.CharacteristicLengthMin = 0.1;
//Mesh.CharacteristicLengthMax = 0.15;
//
//CharacteristicLength = 0.1;
//
//// Apply an elliptic smoother to the grid
//Mesh.Smoothing = 100;
//
////subdivide elements if you whant to (0=none, 1=all quadrangles, 2=all hexahedra)
//Mesh.SubdivisionAlgorithm=0;
//
//// Generate the dual mesh (with triangles you get polyhedral cells)
//Mesh.Dual = 0;




//======================================================================
// DO MESH 
//======================================================================

//OPTIONAL set element order
Mesh.ElementOrder = 2;
//OPTIONAL use quadratic elements of serendepety class
Mesh.SecondOrderIncomplete = 1;

//unkomment this to generate the mesh immediately
Mesh 2;

//YOU ARE DONE. SAVE THE MESH USING THE GUI AND LOAD IT TO YOUR CFS SIMULATION
