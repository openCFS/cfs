# -*- coding: utf-8 -*-
 
def map_serendipity_triangle_to_sphere(input, output):
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
 
    # Mittelpunkt und Radius der Kugel definieren
    center = [1.1, 1.2, 0]
    radius = 1.2
 
    # Koordinaten der sechs Stützpunkte Punkte auf Kugel definieren
    p1 = [radius * 1 + center[0], radius * 0 + center[1], radius * 0 + center[2]]
    p2 = [radius * 0 + center[0], radius * 1 + center[1], radius * 0 + center[2]]
    p3 = [radius * 0 + center[0], radius * 0 + center[1], radius * 1 + center[2]]
    p4 = [radius * 0.70711 + center[0], radius * 0.70711 + center[1], radius * 0 + center[2]]
    p5 = [radius * 0 + center[0], radius * 0.70711 + center[1], radius * 0.70711 + center[2]]
    p6 = [radius * 0.70711 + center[0], radius * 0 + center[1], radius * 0.70711 + center[2]]
 
    for i in range(0, numPoints):
        # Koordinaten von Referenzelement holen
        coord = input.GetPoint(i)
        x, y, z = coord[:3]
 
        # Hilfsvariable und Serendipity Basisfunktionen definieren
        t = 1.0 - x - y
        print(t)
        N1 = t * (2*t - 1);
        N2 = x*(2*x - 1);
        N3 = y*(2*y - 1);
        N4 = 4 * x * t;
        N5 = 4 * x * y;
        N6 = 4 * y * t;
 
        # Neue Koordinaten ausrechnen
        x_new = p1[0]*N1 + p2[0]*N2 + p3[0]*N3 + p4[0]*N4 + p5[0]*N5 + p6[0]*N6
        y_new = p1[1]*N1 + p2[1]*N2 + p3[1]*N3 + p4[1]*N4 + p5[1]*N5 + p6[1]*N6
        z_new = p1[2]*N1 + p2[2]*N2 + p3[2]*N3 + p4[2]*N4 + p5[2]*N5 + p6[2]*N6
 
        # Neue Koordinaten in Ausgabe schreiben
        newPoints.InsertPoint(i, x_new, y_new, z_new)
 
        # Vektorergebnis definieren. Vektor zeigt vom neuen Punkt zum Referenzpunkt
        ca.SetTuple3(i, -(x_new - x), -(y_new - y), -(z_new - z))
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
        map_serendipity_triangle_to_sphere(curInput, curOutput)
        iter.GoToNextItem();
else:
  print("Cannot do that!")
