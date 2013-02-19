// Merge the image (this will create a new post-processing view,
// View[0])
Merge "photo2.jpg";

// Modify the normalized pixel values to make them reasonnable
// characteristic lengths for the mesh
Plugin(Evaluate).Expression = "v0 * 10";
Plugin(Evaluate).Run;

// Apply the view as the current background mesh
Background Mesh View[0];

// Build a simple geometry on top of the background mesh
w = View[0].MaxX;
h = View[0].MaxY;
Point(1)={0,0,0,w};
Point(2)={w,0,0,w};
Point(3)={w,h,0,w};
Point(4)={0,h,0,w};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {3,4,1,2};
Plane Surface(6) = {5};

Physical Surface ( "Foto" ) = { 6 };
