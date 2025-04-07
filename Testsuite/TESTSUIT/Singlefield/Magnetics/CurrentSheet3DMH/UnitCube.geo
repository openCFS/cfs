l = 1.0;
Nx = 1; // elements in x direction
Ny = 1; // elements in y direction
Nz = 10; // elements in z direction

Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Point(1) = {-0.5*l, -0.5*l, -0.5*l, 1};
Point(2) = { 0.5*l, -0.5*l, -0.5*l, 1};
Point(3) = { 0.5*l,  0.5*l, -0.5*l, 1};
Point(4) = {-0.5*l,  0.5*l, -0.5*l, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line{1,3} = Nx+1;
Transfinite Line{2,4} = Ny+1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1};

Extrude {0, 0, l} {Surface{1};Layers{Nz};Recombine;}

// corner points
Physical Point('P_NOB') = {3};
Physical Point('P_NWB') = {4};
Physical Point('P_SOB') = {2};
Physical Point('P_SWB') = {1};
Physical Point('P_NOT') = {10};
Physical Point('P_NWT') = {14};
Physical Point('P_SOT') = {6};
Physical Point('P_SWT') = {5};
// edges
Physical Line('L_NW') = {20};
Physical Line('L_SW') = {11};
Physical Line('L_SE') = {12};
Physical Line('L_NO') = {16};
Physical Line('L_NB') = {3};
Physical Line('L_OB') = {2};
Physical Line('L_SB') = {1};
Physical Line('L_WB') = {4};
Physical Line('L_NT') = {8};
Physical Line('L_OT') = {7};
Physical Line('L_ST') = {6};
Physical Line('L_WT') = {9};
// faces
Physical Surface('S_N') = {21};
Physical Surface('S_O') = {17};
Physical Surface('S_S') = {13};
Physical Surface('S_W') = {25};
Physical Surface('S_B') = {1};
Physical Surface('S_T') = {26};
// volume
Physical Volume('V') = {1};
// move such that z starts at zero
Translate {0, 0, l/2} {
  Volume{1}; 
}
