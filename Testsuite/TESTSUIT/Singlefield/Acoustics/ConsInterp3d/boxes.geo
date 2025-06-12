propBox = 1;
numDivInner = 6;
numDivPML = 3;
pmlW = 0.3;

//create inner points
box2 = propBox/2;
Point (1) = {-box2,-box2,-box2};
Point (2) = {box2,-box2,-box2};
Point (3) = {box2,box2,-box2};
Point (4) = {-box2,box2,-box2};
Point (5) = {-box2,-box2,box2};
Point (6) = {box2,-box2,box2};
Point (7) = {box2,box2,box2};
Point (8) = {-box2,box2,box2};

Line (1) = {1,2};
Line (2) = {2,3};
Line (3) = {4,3};
Line (4) = {1,4};
Line (5) = {1,5};
Line (6) = {2,6};
Line (7) = {3,7};
Line (8) = {4,8};
Line (9) = {5,6};
Line (10) = {6,7};
Line (11) = {8,7};
Line (12) = {5,8};

Line Loop(13) = {9, 10, -11, -12};
Plane Surface(14) = {13};
Line Loop(15) = {6, 10, -7, -2};
Plane Surface(16) = {15};
Line Loop(17) = {11, -7, -3, 8};
Plane Surface(18) = {17};
Line Loop(19) = {5, 12, -8, -4};
Plane Surface(20) = {19};
Line Loop(21) = {9, -6, -1, 5};
Plane Surface(22) = {21};
Line Loop(23) = {1, 2, -3, -4};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 22, 16, 18, 24, 20};
Volume(26) = {25};
Extrude {0, pmlW, 0} {
  Surface{18};
}
Extrude {0, -pmlW, 0} {
  Surface{22};
}
Extrude {pmlW, 0, 0} {
  Surface{16, 39, 61};
}
Extrude {-pmlW, 0, 0} {
  Surface{69, 20, 47};
}
Extrude {0, 0, pmlW} {
  Surface{57, 149, 171, 14, 83, 113, 35, 193, 135};
}
//Extrude {0, 0, -pmlW} {
//  Surface{91, 24, 43, 105, 127, 65, 157, 179, 201};
//}
Coherence;



Transfinite Line {140, 138, 5, 53, 6, 118, 72, 51, 9, 50, 204, 206, 161, 250, 12, 248, 163, 4, 162, 184, 8, 31, 3, 30, 52, 1, 28, 338, 11, 272, 29, 29, 74, 96, 7, 2, 75, 10, 271, 73, 294} = numDivInner Using Progression 1;

Transfinite Line {152, 143, 141, 64, 174, 196, 42, 185, 192, 368, 350, 254, 359, 170, 258, 183, 360, 249, 229, 144, 240, 209, 207, 55, 148, 139, 236, 228, 227, 218, 210, 214, 205, 56, 119, 306, 295, 382, 383, 394, 130, 78, 319, 108, 328, 302, 317, 82, 280, 293, 77, 60, 126, 117, 104, 86, 95, 38, 34, 314, 97, 316} = numDivPML Using Progression 1;
Transfinite Surface '*';
Transfinite Volume{42} = {84, 9, 8, 74, 156, 150, 110, 114};
Transfinite Volume{41} = {9, 10, 7, 8, 150, 131, 120, 110};
Transfinite Volume{40} = {10, 48, 34, 7, 131, 140, 126, 120};
Transfinite Volume{39} = {34, 7, 6, 30, 126, 120, 90, 130};
Transfinite Volume{38} = {7, 8, 5, 6, 120, 110, 89, 90};
Transfinite Volume{37} = {8, 74, 60, 5, 110, 114, 108, 89};
Transfinite Volume{43} = {30, 6, 20, 58, 130, 90,94,166};
Transfinite Volume{35} = {6, 5, 19, 20, 90, 89, 98, 94};
Transfinite Volume{36} = {5, 60, 64, 19, 89, 108, 104, 98};


Transfinite Volume{30} = {44, 14, 3, 38, 48, 10, 7, 34};
Transfinite Volume{27} = {14, 18, 4, 3, 10, 9, 8, 7};
Transfinite Volume{34} = {18, 88, 78, 4, 9, 84, 74, 8};
Transfinite Volume{29} = {38, 3, 2, 29, 34, 7, 6, 30};
Transfinite Volume{26} = {3, 4, 1, 2, 7, 8, 5, 6};
Transfinite Volume{33} = {4, 78, 59, 1, 8, 74, 60, 5};
Transfinite Volume{31} = {29, 2, 24, 54, 30, 6, 20, 58};
Transfinite Volume{28} = {2, 1, 28, 24, 6, 5, 19, 20};
Transfinite Volume{32} = {1, 59, 68, 28, 5, 60, 64, 19};

Physical Volume ("inner") = {26};
Physical Volume ("pml") = {27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52};
Physical Surface ("excite") = {24};


Mesh.RecombineAll=1;
