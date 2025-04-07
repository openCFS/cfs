// kate: indent-width 2; mixedindent on; indent-mode cstyle; line-numbers on;
// kate: syntax c; replace-tabs on;
// kate: remove-trailing-space on; bracket-highlight-color #ff00ff;

Mesh.Smoothing = 10;
//Mesh.SecondOrderLinear = 1;
//Mesh.SubdivisionAlgorithm = 0;
//Mesh.RemeshAlgorithm = 0;
//Mesh.RemeshParametrization = 1;
//Mesh.RecombineAll = 1;
Mesh.Algorithm3D = 6;
Mesh.ElementOrder = 1;
Mesh.Optimize = 1;
//Mesh 3;
Geometry.CopyMeshingMethod=1;


nu = 1e-4;
f = 3.1000000000e+01;
pi = 4*Atan(1);
omega = 2*pi*f;
// Grenzschicht
G = 6*Sqrt(nu / omega);
G_min = 2.0000000000e-03;

// If(G < G_min)
//  G = G_min;
// EndIf


Printf("Kin. Viskosität nu: %f", nu);
Printf("Frequenz f: %f", f);
Printf("pi: %f", pi);
Printf("omega: %f", omega);
Printf("Grenzschichtdicke G: %f", G);

// Radius
r = 1.3150000000e-02;

AllQuadNodesOnRadius = 1;

// Rohrinnenradius
// R_i = R_a - W;
R_i = r;

// Rohraußenradius
R_a = 1.4250000000e-02;

// Rohrwandstärke
W = R_a - R_i;

// Grenzschichtradius
R_g = R_i - G;

// Radius für inneres Quadrat
R_q = 3 * R_g / 4;

// Rohrlaenge
l = 3.0000000000e-01;

// Verfeinerung
feinheit = 2;

N_radial = (2*feinheit)+1;
N_grenz = (6+2*feinheit)+1;
N_laenge = Ceil((l/1e-3) / 10)+1;
N_rohr = 3;
// Floor(feinheit/2)+1;
// prog = 1.0 + 0.75^feinheit;
prog = 1.4;

Printf("Progression prog: %f", prog);

N_ExtrudeZ = N_laenge/2;
ExtrudeHeight = l/2;

Point(1) = {0, -R_a, 0};
Point(2) = {l, -R_a, 0};

Line(1) = {1, 2};
Extrude {0, W, 0} {
  Line{1};
}

Extrude {0, G, 0} {
  Line{2};
}

Extrude {0, 2*R_g, 0} {
  Line{6};
}

Extrude {0, G, 0} {
  Line{10};
}

Extrude {0, W, 0} {
  Line{14};
}

Transfinite Line {3,4,19,20} = N_rohr Using Progression 1.0;
Transfinite Line {1,2,6,10,14,18} = N_laenge Using Progression 1.0;

Transfinite Line {7, 8, 15, 16} = N_grenz Using Progression 1.0;

Transfinite Line {11, 12} = N_radial Using Progression 1.0;

Transfinite Surface{5,9,13,17,21} = {};
Recombine Surface{5,9,13,17,21};

Physical Surface("pipe") = {5,21};
Physical Surface("water") = {9,13,17};

Physical Line("face_water_inflow") = {7,11,15};
Physical Line("face_water_outflow") = {8,12,16};

Physical Line("face_pipe_inflow") = {3,19};
Physical Line("face_pipe_outflow") = {4,20};

Physical Line("pipe_inner") = {2,14};

Physical Line("pipe_outer") = {1,18};
