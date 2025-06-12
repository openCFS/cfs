h = 1.0/10.0; // mesh size
L = 1.0; // domain length
L_pml = 2.0; // PML length
a = 2.0; // gross section length
N_pml = L_pml/h; // number of elements in PML
N_L = L/h; // number of elements in length direction
//+
Point(1) = {a/2, a/2, 0, h};
//+
Point(2) = {0, a/2, 0, h};
//+
Point(3) = {0, 0, 0, h};
//+
Point(4) = {a/2, 0, 0, h};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Line{1,2,3,4} = 1;
prop[] = Extrude {0, 0, L} {
  Surface{1}; Layers{N_L}; Recombine;
};
//+
pml[] = Extrude {0, 0, L_pml} {
  Surface{prop[0]}; Layers{N_pml}; Recombine;
};
Mesh.RecombineAll = 1;
Mesh.SecondOrderIncomplete = 1;
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.ElementOrder = 1; // use second order elements

Transfinite Surface{1,prop[0],prop[2],prop[3],prop[4],prop[5]};
//+
Physical Surface("sym-xz") = {prop[4],pml[4]};
Physical Surface("sym-yz") = {prop[3],pml[3]};
Physical Surface("outer-xz") = {prop[2],pml[2]};
Physical Surface("outer-yz") = {prop[5],pml[5]};
//+
Physical Volume("prop") = {prop[1]};
//+
Physical Volume("pml") = {pml[1]};
//+
Physical Surface("excite") = {1};
//+
Physical Surface("end") = {48};
//+
Physical Line("length") = {11, 33};
//+
Physical Point("fix-xy") = {1};
//+
Physical Point("fix-x") = {4};
//+
Physical Point("fix-y") = {2};
