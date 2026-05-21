Extrude a 2D Mesh with linear line/surface elements (BAR2, TRI3 and QUAD4) to a 3D Mesh with linear surface/volume elements (QUAD4, WEDGE6 and HEX8) for z = "-1.0 -0.5 0.5 1.0"

Test the extruded 3D mesh with h5diff

h5diff usage:
### Use h5diff
-DTEST_H5DIFF:BOOL=ON
### which object in the hdf5 shoule be compared
### In this example we compare the /Mesh 
-DH5DIFF_OBJ:STRING=/Mesh
