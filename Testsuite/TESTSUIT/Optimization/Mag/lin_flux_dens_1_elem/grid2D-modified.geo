// GMSH geometry definition
// vim: set filetype=c : //


// set view options
Geometry.LineNumbers = 0;
Geometry.PointNumbers = 0;
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 0; // show surface numbers
Geometry.VolumeNumbers = 0; // show volume numbers

// Dimension

r_1 = 0.005;
r_2 = 0.01;
l_1 = 0.008;
d = 0.02; // distance to outer domain

h1 = 50e-3;
h2 = h1+15e-3;
h3 = h2+1e-3;
h4 = h3+5e-3;
h5 = h4+10e-3;
h6 = h5+20e-3;
h7 = h6+10e-3;
h8 = h7+50e-3;
b1 = 55e-3;
b2 = b1+20e-3;
b3 = b2+50e-3;
b4 = b3+20e-3;
b5 = b4+55e-3;

Geometry.PointNumbers = 0;


// coordinates must be given in sorted vectors
x[] = {0.0, b1, b2, b3, b4, b5}; // x-ccordinates
y[] = {0, h1, h2, h3, h4, h5, h6, h7, h8}; // y-coodinates

h = 0.005; // mesh size

// create points of the grid
For i In {1:(#x[])}
  For j In {1:(#y[])}
    // The 'newp' command returns the next free Point-ID
    n = newp;
    Point(n) = {x[i-1],y[j-1],0,h};
    // If a string-expression is followd by "~{exp}" is appended by "_val" where "val" is the evaluation of the expression "exp" 
    P~{i}~{j} = n; // creates the variables P_1_1, P_1_2, ...
    Printf("Point %g = (%g,%g)",n,x[i-1],y[j-1]);
  EndFor
EndFor

// create y-|| lines
For i In {1:(#x[])}
  For j In {1:(#y[]-1)}
    n = newl;
    Line(n) = {P~{i}~{j},P~{i}~{j+1}};
    Ly~{i}~{j} = n;
    Ly~{i}[j] = n;
  EndFor
EndFor

// create x-|| lines
For j In {1:(#y[])}
  For i In {1:(#x[]-1)}
    n = newl;
    Line(n) = {P~{i}~{j},P~{i+1}~{j}};
    Lx~{i}~{j} = n;
    Lx~{j}[i] = n;
  EndFor
EndFor

// create surfaces
For i In {1:(#x[]-1)}
  For j In {1:(#y[]-1)}
    n = news;
    Line Loop(n) = { Lx~{i}~{j}, Ly~{i+1}~{j}, -Lx~{i}~{j+1}, -Ly~{i}~{j} };
    Plane Surface(n) = {n};
    S~{i}~{j} = n;
    Transfinite Surface{n};
  EndFor
EndFor


Mesh.RecombineAll = 1;
Mesh.ElementOrder = 1;

Transfinite Line {4,12,20,28,36,44} = 1;
Transfinite Line {50,55,60,65,70,75,80,85,90} = 1;

Physical Surface("joch") = {S_2_5, S_2_6, S_3_6, S_4_6, S_4_5,S_4_4};
Physical Surface("anker") = {S_2_2, S_3_2, S_4_2};
Physical Surface("spule1_in") = {S_3_7};
Physical Surface("spule2_out") = {S_3_5};
Physical Surface("air") = {S_1_1, S_1_2, S_1_3, S_1_4, S_1_5, S_1_6, S_1_7, S_1_8, S_2_1, S_2_3, S_2_7, S_2_8, S_3_1, S_3_3, S_3_4, S_3_7, S_3_8, S_4_1, S_4_3, S_4_7, S_4_8, S_5_1, S_5_2, S_5_3, S_5_4, S_5_5, S_5_6, S_5_7, S_5_8};
Physical Surface("OptimierungY") = {S_2_4};

Physical Line("FluxParallel_V") = { Lx_1[{1:(#x[]-1)}], Lx~{#y[]}[{1:(#x[]-1)}] }; // horizontal lines
Physical Line("FluxParallel_h") = { Ly_1[{1:(#y[]-1)}], Ly~{#x[]}[{1:(#y[]-1)}] }; // the vertical lines at the left (x=x1)
Physical Line("Force") = { Lx_2_3, Lx_3_3, Lx_4_3, Lx_2_2, Lx_3_2, Lx_4_2, Ly_5_2, Ly_2_2};
//Physical Line("OptimierungY") = {60};

//Physical Line("Ly_3") = { Ly_3_1, Ly_3_2, Ly_3_3 };
