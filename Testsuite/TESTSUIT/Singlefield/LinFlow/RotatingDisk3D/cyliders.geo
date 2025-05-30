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
r[0] = 1;
//r[1] = 3e-2;
//r[2] = 20e-2;
r[1] = 2;

h = 0.1; // mesh size

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
  Transfinite Surface{ s~{i} };
  Transfinite Line{ lx~{i}, ly~{i} } = (r[i]-r[i-1])/h;
  lc[] += { c~{i} };
  out[] = Extrude {0,0,0.2} { Surface{s~{i}}; Layers{1}; Recombine; };
  Physical Surface(Sprintf("S_%.0ff",i)) = { s~{i} };
  Physical Surface(Sprintf("S_%.0fb",i)) = { out[0] };
  Physical Volume(Sprintf("V_%.0f",i)) = { out[1] };
  Physical Surface(Sprintf("S_%.0fx",i)) = { out[2] };
  Physical Surface(Sprintf("S_%.0fy",i)) = { out[4] };
  Physical Surface(Sprintf("S_%.0f",i)) = { out[3] };
  Physical Surface(Sprintf("S_%.0f",i-1)) = { out[5] };
EndFor

Transfinite Line{ lc[] } = 16;

// specify mesh parameters
Mesh.RecombineAll = 1; // recombine triangles to quads
Mesh.RecombinationAlgorithm = 1; // Blossom
Mesh.ElementOrder = 2; // use second order elements
Mesh.SecondOrderIncomplete = 1; // use 20-noded hexas
