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
Coherence;
