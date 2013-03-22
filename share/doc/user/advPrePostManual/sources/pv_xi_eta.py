def map_quadrilateral(input, output):
    # Copy the cells etc.
    output.ShallowCopy(input)
    newPoints = vtk.vtkPoints()
    numPoints = input.GetNumberOfPoints()
    values = input.GetPointData().GetScalars()
    ca = vtk.vtkFloatArray()
    ca.SetName('MapVector')
    ca.SetNumberOfComponents(3)
    ca.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(ca)

    # Create vectors for basis functions on mapped element
    nv1 = vtk.vtkFloatArray()
    nv1.SetName('N1')
    nv1.SetNumberOfComponents(1)
    nv1.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(nv1)

    nv2 = vtk.vtkFloatArray()
    nv2.SetName('N2')
    nv2.SetNumberOfComponents(1)
    nv2.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(nv2)

    nv3 = vtk.vtkFloatArray()
    nv3.SetName('N3')
    nv3.SetNumberOfComponents(1)
    nv3.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(nv3)

    nv4 = vtk.vtkFloatArray()
    nv4.SetName('N4')
    nv4.SetNumberOfComponents(1)
    nv4.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(nv4)

    # Create vectors for the mapped basis vectors
    xi = vtk.vtkFloatArray()
    xi.SetName('XiVector')
    xi.SetNumberOfComponents(3)
    xi.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(xi)

    eta = vtk.vtkFloatArray()
    eta.SetName('EtaVector')
    eta.SetNumberOfComponents(3)
    eta.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(eta)

    zeta = vtk.vtkFloatArray()
    zeta.SetName('ZetaVector')
    zeta.SetNumberOfComponents(3)
    zeta.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(zeta)

    # Create vector for Jacobian determinant of mapping
    jacDet = vtk.vtkFloatArray()
    jacDet.SetName('jacDet')
    jacDet.SetNumberOfComponents(1)
    jacDet.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(jacDet)

    # Define the four corner coordinates
    p1 = [0, 1, 0.75]
    p2 = [0.25, 1, 0.5749999881]
    p3 = [0.2797259986, 1, 0.5952895284]
    p4 = [0.2693634927, 1, 0.6263765097]

    for i in range(0, numPoints):
        # Get coordinates from reference element
        coord = input.GetPoint(i)
        x, y, z = coord[:3]
 
        # Define nodal basis functions
        N1 = 0.25 * (-x + 1) * (-y + 1);
        N2 = 0.25 * (x + 1) * (-y + 1);
        N3 = 0.25 * (x + 1) * (y + 1);
        N4 = 0.25 * (-x + 1) * (y + 1);

        # Define derivatives of basis functions
        N1_dxi = -0.25 * (-y + 1);
        N2_dxi = 0.25 * (-y + 1);
        N3_dxi = 0.25 * (y + 1);
        N4_dxi = -0.25 * (y + 1);

        N1_deta = -0.25 * (-x + 1);
        N2_deta = -0.25 * (x + 1);
        N3_deta = 0.25 * (x + 1);
        N4_deta = 0.25 * (-x + 1);
 
        # Compute new coordinates
        x_new = p1[0]*N1 + p2[0]*N2 + p3[0]*N3 + p4[0]*N4
        y_new = p1[1]*N1 + p2[1]*N2 + p3[1]*N3 + p4[1]*N4
        z_new = p1[2]*N1 + p2[2]*N2 + p3[2]*N3 + p4[2]*N4
 
        # Insert new coordinates into newPoints vector
        newPoints.InsertPoint(i, x_new, y_new, z_new)
 
        # Compute inverse mapping vector from global to reference point 
        ca.SetTuple3(i, -(x_new - x), -(y_new - y), -(z_new - z))

        # Fill vectors for basis functions
        nv1.SetValue(i, N1)
        nv2.SetValue(i, N2)
        nv3.SetValue(i, N3)
        nv4.SetValue(i, N4)

        # Compute column vectors of Jacabian matrix
        J_x_dxi = p1[0]*N1_dxi + p2[0]*N2_dxi + p3[0]*N3_dxi + p4[0]*N4_dxi
        J_y_dxi = p1[1]*N1_dxi + p2[1]*N2_dxi + p3[1]*N3_dxi + p4[1]*N4_dxi
        J_z_dxi = p1[2]*N1_dxi + p2[2]*N2_dxi + p3[2]*N3_dxi + p4[2]*N4_dxi

        J_x_deta = p1[0]*N1_deta + p2[0]*N2_deta + p3[0]*N3_deta + p4[0]*N4_deta
        J_y_deta = p1[1]*N1_deta + p2[1]*N2_deta + p3[1]*N3_deta + p4[1]*N4_deta
        J_z_deta = p1[2]*N1_deta + p2[2]*N2_deta + p3[2]*N3_deta + p4[2]*N4_deta

        normal_x = + J_y_dxi * J_z_deta - J_z_dxi * J_y_deta;
        normal_y = - J_x_dxi * J_z_deta + J_z_dxi * J_x_deta;
        normal_z = + J_x_dxi * J_y_deta - J_y_dxi * J_x_deta;

        # Compute Jacobian determinant
        jD = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)
        jacDet.SetValue(i, jD)

        # Fill mapped vectors
        xi.SetTuple3(i, J_x_dxi, J_y_dxi, J_z_dxi)
        eta.SetTuple3(i, J_x_deta, J_y_deta, J_z_deta)
        zeta.SetTuple3(i, normal_x/jD, normal_y/jD, normal_z/jD)

    output.SetPoints(newPoints)
 
input = self.GetInputDataObject(0, 0)
output = self.GetOutputDataObject(0)
 
if input.IsA("vtkMultiBlockDataSet"):
    output.CopyStructure(input)
    iter = input.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    while not iter.IsDoneWithTraversal():
        curInput = iter.GetCurrentDataObject()
        curOutput = curInput.NewInstance()
        curOutput.UnRegister(None)
        output.SetDataSet(iter, curOutput)
        map_quadrilateral(curInput, curOutput)
        iter.GoToNextItem();
else:
  print("Cannot do that!")

