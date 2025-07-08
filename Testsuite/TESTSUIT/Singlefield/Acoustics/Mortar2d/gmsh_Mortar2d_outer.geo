// Gmsh  (http://geuz.org/gmsh)  model  for  an acoustic propagation region
// PML region.
//
// To create the mesh run
// gmsh gmsh_Mortar2d_other.geo -2 -v 0 -format msh -order 2 -o l2d_outer.msh

esize = 25e-2;

Point(1) = {  0, 0, 0, esize};
Point(2) = {0.5, 0, 0, esize};
Point(3) = {1.0, 0, 0, esize};
Point(4) = {1.5, 0, 0, esize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Extrude {0, 0.5, 0} {
  Line{1, 2, 3};
}
Extrude {0, 0.5, 0} {
  Line{4, 8, 12};
}
Extrude {0, 0.5, 0} {
  Line{16, 20, 24};
}

// Physical Surface("inner") = {7};
Physical Surface("outer") = {11, 23, 19};
Physical Surface("pml") = {31, 35, 39, 27, 15};
Physical Line("outer_iface1") = {6};
Physical Line("outer_iface2") = {4};
