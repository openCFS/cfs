// mesh of a quater pipe
R = 0.001;
L = 0.1;
h = R/10;
hL = 0.03/15; 
// for single element
// h = 2*R;
// hL = L;
//+
Point(1) = {0, 0, 0, h};
Point(2) = {R, 0, 0, h};
Point(3) = {0, R, 0, h};
//+
Circle(1) = {2, 1, 3};
//+
Line(2) = {3, 1};
//+
Line(3) = {1, 2};
//+
Curve Loop(1) = {1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Recombine Surface {1};
out[] = Extrude {0, 0, L} {
  Surface{1}; Layers{ Round(L/hL) }; Recombine;
};

// regions
Physical Surface("S_in") = {1};
Physical Surface("S_out") = { out[0] };
Physical Surface("S_wall") = { out[2] };
Physical Surface("S_sym-x") = { out[3] };
Physical Surface("S_sym-y") = { out[4] };

Physical Volume("V_pipe") = { out[1] };

// acou domain
outA[] = Extrude {0, 0, L/2} {
  Surface{ out[0] }; Layers{ Round(L/2/hL) }; Recombine;
};
Physical Surface("S_end") = { outA[0] };
Physical Surface("S_Awall") = { outA[2] };
Physical Surface("S_Asym-x") = { outA[3] };
Physical Surface("S_Asym-y") = { outA[4] };
Physical Volume("V_Apipe") = { outA[1] };


