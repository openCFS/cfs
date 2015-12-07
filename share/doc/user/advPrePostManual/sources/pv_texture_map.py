# Obtain poly data input and output 
pdi = self.GetPolyDataInput()
pdo = self.GetPolyDataOutput()

# Determine number of points
numPts = pdi.GetPoints().GetNumberOfPoints()

# Texture coordinates may also be loaded from data arrays
#Uarray = pdi.GetPointData().GetArray('U')
#Varray = pdi.GetPointData().GetArray('V')

# Create new array for texture coordinates
tc = vtk.vtkDoubleArray()
tc.SetNumberOfComponents(2)
tc.SetNumberOfTuples(numPts)
tc.SetName('TCoords')

# Loop over all points and assign new texture coordinates
for ii in range(numPts):
#  Uval = Uarray.GetTuple1(ii)
#  Vval = Varray.GetTuple1(ii)
  coord = pdi.GetPoint(ii)
  x, y, z = coord[:3]
  Uval = x*y
  Vval = z*y
  tc.SetTuple2(ii,Uval,Vval)

# Copy input data set to output and attach texture coordinate array to it
pdo.ShallowCopy(pdi)
pdo.GetPointData().AddArray(tc)
pdo.GetPointData().SetActiveTCoords('TCoords')
