boxi = 0.2;
ndiv = 3;

box2 = boxi/2;
offsetx = 0.123;
offsety = 0.123;
offsetz = 0.123;

Point (1) = {-box2+offsetx, -box2+offsety,-box2+offsetz};
Point (2) = { box2+offsetx, -box2+offsety,-box2+offsetz};
Point (3) = { box2+offsetx,  box2+offsety,-box2+offsetz};
Point (4) = {-box2+offsetx,  box2+offsety,-box2+offsetz};
Point (5) = {-box2+offsetx, -box2+offsety, box2+offsetz};
Point (6) = { box2+offsetx, -box2+offsety, box2+offsetz};
Point (7) = { box2+offsetx,  box2+offsety, box2+offsetz};
Point (8) = {-box2+offsetx,  box2+offsety, box2+offsetz};

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

Transfinite Line {9, 10, 8, 11, 12, 1, 6, 2, 3, 7, 4, 5} = ndiv Using Progression 1;
Transfinite Surface '*';
Transfinite Volume{26} = {5, 6, 7, 8, 1, 2, 3, 4};

Physical Volume("srcReg") = {26};

Mesh.RecombineAll =1;
Mesh 3;
