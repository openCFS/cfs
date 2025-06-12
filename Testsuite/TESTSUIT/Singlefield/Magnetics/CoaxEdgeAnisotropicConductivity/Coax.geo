t = 2; // model thickness
h = 0.1; // element size

Mesh.ElementOrder = 1;
Mesh.RecombineAll = 1;

Geometry.LineNumbers = 1;
Geometry.PointNumbers = 1;
Geometry.SurfaceNumbers = 1;

// radial divisions
R[0] = 0;
R[1] = 4;
R[2] = 4.5;
R[3] = 8;
Np = 10; // elements in circumferential direction
Nt = 1; // elements in thickness direction
N[0] = 15;
N[1] = 2;
N[2] = 5;
p[0] = 0.9;
p[1] = 1.0;
p[2] = 1.0;


// points
For i In {0 : #R[]-1 }
  Point(i) = {R[i],0,t/2,h};
EndFor

// build model
For i In {0 : #R[]-2 }
  nl = newl;
  Line(nl) = {i,i+1};
  Transfinite Line {nl} = N[i]+1 Using Progression p[i];
  outS[] = Extrude { {0,0,1} , {0,0,0} , Pi/2 } {
   Line{nl}; Layers{Np}; Recombine;
    };
  //Physical Region{V~{i+1}} = { out[1] };
  //Physical Line(Sprintf("L_%g-x",i+1)) = { nl };
  //Physical Line(Sprintf("L_%g-y",i+1)) = { outS[0] };
  //Physical Line(Sprintf("L_%g-r",i+1)) = { outS[2] };
  Physical Surface(Sprintf("S_%g-f",i+1)) = { outS[1] };
  // Extrude the volume
  nv = newv;
  outV[] = Extrude { 0,0,-t } {
   Surface{ outS[1] }; Layers{Nt}; Recombine;
    };
  Physical Surface(Sprintf("S_%g-b",i+1)) = { outV[0] };
  Physical Surface(Sprintf("S_%g-x",i+1)) = { outV[2] };
  Physical Surface(Sprintf("S_%g-r",i+1)) = { outV[3] };
  Physical Surface(Sprintf("S_%g-y",i+1)) = { outV[4] };
  Physical Volume(Sprintf("V_%g",i+1)) = { outV[1] };
EndFor
