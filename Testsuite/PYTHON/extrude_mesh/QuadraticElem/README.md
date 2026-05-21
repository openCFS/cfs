Extrude a 2D Mesh with Quadratic Element (BAR3, TRI6 and QUAD8) to a 3D Mesh with Quadratic Elements (QUAD8, WEDGE15 and HEX20) for z = "-1.0 -0.5 0.5 1.0"

Test the extruded 3D mesh with h5diff

h5diff usage:
### Use h5diff
-DTEST_H5DIFF:BOOL=ON
### which object in the hdf5 shoule be compared
### In this example we compare the /Mesh 
-DH5DIFF_OBJ:STRING=/Mesh
