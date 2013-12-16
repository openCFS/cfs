import vtk
from vtk.util.colors import *
from numpy import *
from matviz_rot import *

## creates 3D data to vtkPolyData
def create_vtk_poly_data(angle, data):
  
  assert(len(angle) == len(data))
  samples = int(sqrt(len(angle)))
  assert(samples**2 == len(angle))
  
  pts = vtk.vtkPoints()
  
  for i in range(len(angle)):
    radius = data[i]
    
    # we use not the usual conversion for the spherical coordinate system but the first row of the rotation matrix!
    p = radius * to_vector(angle[i])
    
    # print "phi="  + str(phi) + " theta=" + str(theta) + " radius=" + str(radius) + " -> x=" + str(x) + " y=" + str(y) + " z=" + str(z)
    pts.InsertNextPoint(p[0], p[1], p[2])

  pdo = vtk.vtkPolyData()
  pdo.SetPoints(pts)
  #Allocate memory for triangles
  pdo.Allocate((samples-1) * (samples-1) * 2)
  
  #Generate two triangles per quad
  for i in range(0, samples-1):
    for j in range(0, samples-1):
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, i * samples + j)
      tri.GetPointIds().SetId(1, i * samples + j+1)
      tri.GetPointIds().SetId(2, (i+1) * samples + j)
      # print "i=" + str(i) + " j=" + str(j) + " ids=" + str(tri.GetPointIds())
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
      
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, (i+1) * samples + j+1)
      tri.GetPointIds().SetId(1, (i+1) * samples + j)
      tri.GetPointIds().SetId(2, i * samples + j+1)
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
  
  return pdo

## create a list of vtk actors displaying symmetry planesfrom vtk.util.colors import
# @param minim
def create_symmety_planes(minima, scale, add_planes):
  # code source: http://code.google.com/p/pythonxy/source/browse/src/python/vtk/DOC/Examples/Rendering/Cylinder.py
 
  minima = []
  minima.append(((0, 0), 1))
  minima.append(((0, numpy.pi/2), 1))
  minima.append(((numpy.pi/2, numpy.pi/2), 1))
 
  actors = [] 
   
  for i in range(len(minima)):
     
    # Create cylinder
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(scale / 200)
    cylinder.SetCenter( 0, 0, 0 )
    cylinder.SetHeight(2 * scale)
    cylinder.SetResolution(10)
    cylinder.Update()
     
    # The mapper is responsible for pushing the geometry into the graphics
    # library. It may also do color mapping, if scalars or other attributes are defined.
    mapper_c = vtk.vtkPolyDataMapper()
    mapper_c.SetInputConnection(cylinder.GetOutputPort())

    # The actor is a grouping mechanism: besides the geometry (mapper), it
    # also has a property, transformation matrix, and/or texture map.
    actor_c = vtk.vtkActor()
    actor_c.SetMapper(mapper_c)
        
    actor_c.GetProperty().SetColor(black)
    angle = minima[i][0]
    phi   = angle[0]
    theta = angle[1]

    
    print "angle: " + str(angle) + " -> " + str(to_vector(angle)) + " = " + str(minima[i][1])
    
    #actor_c.RotateX(90)
    #actor_c.RotateY(angle[0] * 180/numpy.pi)
    #actor_c.RotateZ(angle[1] * 180/numpy.pi)


    # cylinders!!!
    actor_c.RotateZ(90)
    actor_c.RotateY(angle[0] * 180/numpy.pi)
    actor_c.RotateX(angle[1] * 180/numpy.pi)


    actors.append(actor_c)

    if add_planes:
      disk = vtk.vtkDiskSource()
      disk.SetInnerRadius(0.1)
      disk.SetOuterRadius(0.9 * scale)
      disk.SetRadialResolution(30)
      disk.SetCircumferentialResolution(60)
      disk.Update()
  
      mapper_d = vtk.vtkPolyDataMapper()
      mapper_d.SetInputConnection(disk.GetOutputPort())
      
      actor_d = vtk.vtkActor()
      actor_d.SetMapper(mapper_d)

      actor_d.GetProperty().SetColor(banana)
      actor_d.GetProperty().SetOpacity(0.5)

      #actor_d.RotateX(angle[0] * 180/numpy.pi)
      #actor_d.RotateY(angle[1] * 180/numpy.pi)
      #actor_d.RotateZ(0)

      # plane!
      actor_d.RotateY(90)
      actor_d.RotateZ(90)
      actor_d.RotateY(angle[0] * 180/numpy.pi)
      actor_d.RotateX(angle[1] * 180/numpy.pi)
      

      
      actors.append(actor_d)
    
  return actors


# show the data on the screen
# @planes list of vtk actors containing symmetry planes 
def show_vtk(polydata, res, planes = []):
  #Create a mapper and actor
  mapper = vtk.vtkPolyDataMapper()
  if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(polydata)
  else:
    mapper.SetInputData(polydata)
  
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(0.5,0.5,0.5) # (R,G,B)
  # Setup a renderer, render window, and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  #renderWindow.SetWindowName("Test")
   
  renderWindow.AddRenderer(renderer);
  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)
   
  renderWindow.SetSize(res,res)
   
  #Add the actor to the scene
  for i in range(len(planes)):
    renderer.AddActor(planes[i])

  renderer.AddActor(actor)
  renderer.SetBackground(1,1,1) # Background color white
   

  #Render and interact
  renderWindow.Render()
   
  #*** SetWindowName after renderWindow.Render() is called***
  #renderWindow.SetWindowName("Test")
   
  renderWindowInteractor.Start()


# helper for create_frame
# @param cells  vtk.vtkCellArray() where cells are added via InsertNextCell
# @param points vtk.vtkPoints() where the points are added
# @param dir 'x', 'y', 'z'
def create_centered_bar(cells, points, center, length, width, height):
  
  #print "ccb: c=" + str(center) + " l=" + str(length) + " t=" + str(thick) + " dir=" + dir  
  
  # Add the points to a vtkPoints object
  dx = 0.5 * length
  dy = 0.5 * width
  dz = 0.5 * height  

  cx, cy, cz = center
  #points = vtk.vtkPoints()
  
  base = points.GetNumberOfPoints()
  
  points.InsertNextPoint([cx-dx,cy-dy,cz-dz]) # 0
  points.InsertNextPoint([cx+dx,cy-dy,cz-dz]) # 1
  points.InsertNextPoint([cx+dx,cy+dy,cz-dz]) # 2
  points.InsertNextPoint([cx-dx,cy+dy,cz-dz]) # 3
  points.InsertNextPoint([cx+dx,cy-dy,cz+dz]) # 4
  points.InsertNextPoint([cx-dx,cy-dy,cz+dz]) # 5
  points.InsertNextPoint([cx+dx,cy+dy,cz+dz]) # 6
  points.InsertNextPoint([cx-dx,cy+dy,cz+dz]) # 7
  
  # Create a cell array to store the quad in
  #quads = vtk.vtkCellArray()
  
  # Create a quad on the four points
  # a
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 0)
  quad.GetPointIds().SetId(1, base + 1)
  quad.GetPointIds().SetId(2, base + 2)
  quad.GetPointIds().SetId(3, base + 3)
  cells.InsertNextCell(quad)
  
  # b
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 1)
  quad.GetPointIds().SetId(1, base + 4)
  quad.GetPointIds().SetId(2, base + 6)
  quad.GetPointIds().SetId(3, base + 2)
  cells.InsertNextCell(quad)
  
  # c
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 4)
  quad.GetPointIds().SetId(1, base + 6)
  quad.GetPointIds().SetId(2, base + 7)
  quad.GetPointIds().SetId(3, base + 5)
  cells.InsertNextCell(quad)
  
  # d
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 2)
  quad.GetPointIds().SetId(1, base + 6)
  quad.GetPointIds().SetId(2, base + 7)
  quad.GetPointIds().SetId(3, base + 3)
  cells.InsertNextCell(quad)
  
  # e
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 0)
  quad.GetPointIds().SetId(1, base + 5)
  quad.GetPointIds().SetId(2, base + 7)
  quad.GetPointIds().SetId(3, base + 3)
  cells.InsertNextCell(quad)
  
  # f
  quad = vtk.vtkQuad()
  quad.GetPointIds().SetId(0, base + 0)
  quad.GetPointIds().SetId(1, base + 1)
  quad.GetPointIds().SetId(2, base + 4)
  quad.GetPointIds().SetId(3, base + 5)
  cells.InsertNextCell(quad)

## without rotation and shearing
def create_3d_frame(coords, s1, s2, s3, dir, scale):

  centers, min, max, elem = coords
  
  dx = elem[0]
  dy = elem[1] 
  dz = elem[2]

  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0

  for i in range(len(s1)):
    coord = centers[i]
    #print "s1=" + str(s1[i]) + " s2=" + str(s2[i]) + " s3=" + str(s3[i])
    #s1[i] *= 0.5
    #s2[i] *= 0.5
    #s3[i] *= 0.5
    if dir == 'horizontal' or dir == 'all':
      create_centered_bar(cells, points, coord, scale * dx, scale * s1[i] * dy, scale * s1[i] * dz)
    if dir == 'vertical' or dir == 'all':
      create_centered_bar(cells, points, coord, scale * s2[i] * dx, scale * dy, scale * s2[i] * dz)
    if dir == 'sagittal' or dir == 'all':
      create_centered_bar(cells, points, coord, scale * s3[i] * dx, scale * s3[i] * dy, scale * dz)
    
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)
    
  return polydata


## litte helper
# @param save filename or none
# @param list which might be empty
def show_write_vtk(poly, res, save, actors = []):
  if save:
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(poly)
    else:
      writer.SetInputData(poly)
    writer.SetFileName(save)
    writer.Write()    
  else:
    show_vtk(poly, res, actors)  
