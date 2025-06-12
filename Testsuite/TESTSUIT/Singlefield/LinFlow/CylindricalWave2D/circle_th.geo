// GMSH geometry definition
// vim: set filetype=c : //

// specify the dimensions & discretisation


// set view options
Geometry.LineNumbers = 1;
Geometry.PointNumbers = 1;
Geometry.Surfaces = 1; // show surfaces
Geometry.SurfaceNumbers = 1; // show surface numbers
Geometry.VolumeNumbers = 1; // show volume numbers

// radii
r[0] = 10e-2;
//r[1] = 3e-2;
//r[2] = 20e-2;
r[1] = 30e-2;

h = 2.0e-3; // mesh size

// create the points
Point(1) = { 0,  0,  0, h};
For i In { 0 : (#r[]-1) }
  px~{i} = newp;
  Point(px~{i}) = { r[i],  0,  0, h};
  py~{i} = newp;
  Point(py~{i}) = { 0, r[i],  0, h};
EndFor

// create the lines and surfaces
lx_0 = newl;
Line(lx_0) = {1,px_0}; // from origin in x-> 
lr[] = { lx_0 };
ly_0 = newl;
Line(ly_0) = {1,py_0};
lr[] += { ly_0 };
c_0 = newl;
Circle(c_0) = {px_0, 1, py_0};
lc[] = { c_0 };
// first surface
s_0 = news;
Line Loop(s_0) = {lx_0,c_0,-ly_0};
Plane Surface(s_0) = {s_0};
For i In { 1 : (#r[]-1) }
  c~{i} = newl;
  Circle(c~{i}) = {px~{i}, 1, py~{i}};
  lx~{i} = newl;
  Line(lx~{i}) = {px~{i-1},px~{i}};
  ly~{i} = newl;
  Line(ly~{i}) = {py~{i-1},py~{i}};
  ll = newl;
  Line Loop(ll) = {lx~{i}, c~{i}, -ly~{i}, -c~{i-1}};
  s~{i} = news;
  Plane Surface(s~{i}) = {ll};
  Physical Surface(Sprintf("S_%.0f",i)) = {s~{i}};
  Transfinite Line{ lx~{i}, ly~{i} } = (r[i]-r[i-1])/h;
  lc[] += { c~{i} };
EndFor

Transfinite Line{ lc[] } = 16;

For i In { 0 : (#r[]-1) }
  Physical Line(Sprintf("L_c%.0f",i)) = {c~{i}};
  Transfinite Surface{ s~{i} };
  Physical Line(Sprintf("L_x%.0f",i)) = {lx~{i}};
  Physical Line(Sprintf("L_y%.0f",i)) = {ly~{i}};  
EndFor

// specify mesh parameters
Mesh.RecombineAll = 1; // recombine triangles to quads
Mesh.RecombinationAlgorithm=1; // Blossom
Mesh.ElementOrder = 2; // use second order elements
Mesh.SecondOrderIncomplete = 1; // use 20-noded hexas
