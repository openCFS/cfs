radius = 0.3;
box = 1;

px = Asin(45*Pi/180)*radius;
py = Acos(45*Pi/180)*radius;

Point(1) = {0,0,0,1};
Point(2) = {px,py,0,1};
Point(3) = {px,-py,0,1};
Point(4) = {-px,-py,0,1};
Point(5) = {-px,py,0,1};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

//create rectangle surrounding
Point (6) = {-box,-box,0,1};
Point (7) = {box,-box,0,1};
Point (8) = {box,box,0,1};
Point (9) = {-box,box,0,1};

Line(5) = {6, 9};
Line(6) = {9, 8};
Line(7) = {8, 7};
Line(8) = {7, 6};
Line(9) = {3, 7};
Line(10) = {4, 6};
Line(11) = {5, 9};
Line(12) = {2, 8};

Line Loop(13) = {9, -7, -12, 1};
Ruled Surface(14) = {13};
Line Loop(15) = {4, 12, -6, -11};
Ruled Surface(16) = {15};
Line Loop(17) = {10, 5, -11, -3};
Ruled Surface(18) = {17};
Line Loop(19) = {9, 8, -10, -2};
Ruled Surface(20) = {19};

Physical Line ("innerR") = {7};
Physical Line ("innerL") = {5};
Physical Line ("innerT") = {6};
Physical Line ("innerB") = {8};
Physical Line ("source") = {1,2,3,4};

//Transfinite Line {1,2,3,4} = 10 Using Progression 1;
//Transfinite Line {5,6,7,8} = 20 Using Progression 1;
Physical Surface ("inner") = {16,14,-18,-20};
//Mesh.CharacteristicLengthFactor = 0.01
//Mesh.CharacteristicLengthMax = 0.02
