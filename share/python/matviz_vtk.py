import vtk
from vtk.util.colors import *
from numpy import *
from matviz_rot import *
import scipy.interpolate as ip 

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
# @param dim list of length, width, height
# @param angle list of angle_x, angle_y, angle_z or None
def create_centered_bar(cells, points, center, dim, angle = None):
  
  #print "ccb: c=" + str(center) + " l=" + str(length) + " t=" + str(thick) + " dir=" + dir  
  length, width, height = dim
  angle_x, angle_y, angle_z = (0.0, 0.0, 0.0) if angle is None else angle
  
  
  cthetax = cos(angle_x)
  sthetax = sin(angle_x)
  cthetay = cos(angle_y)
  sthetay = sin(angle_y)
  cthetaz = cos(angle_z)
  sthetaz = sin(angle_z)
  
  # see DesignMaterial::SetRotationMatrix()
  R = zeros((3,3))
  R[0][0] =  cthetay * cthetaz
  R[0][1] =  -cthetay * sthetaz
  R[1][0] =  cthetax * sthetaz + sthetax * sthetay * cthetaz
  R[1][1] =  cthetax * cthetaz - sthetax * sthetay * sthetaz
  R[0][2] =  sthetay
  R[1][2] = -sthetax * cthetay
  R[2][0] =  sthetax * sthetaz - cthetax * sthetay * cthetaz
  R[2][1] =  sthetax * cthetaz + cthetax * sthetay * sthetaz
  R[2][2] =  cthetax * cthetay

  cx, cy, cz = center
  #points = vtk.vtkPoints()
  
  base = points.GetNumberOfPoints()
  
  
  for x in [(-1.0,-1.0,-1.0),(1.0,-1.0,-1.0),(1.0,1.0,-1.0),(-1.0,1.0,-1.0),(1.0,-1.0,1.0),(-1.0,-1.0,1.0),(1.0,1.0,1.0),(-1.0,1.0,1.0)]:
    p = (x[0] * 0.5 * length, x[1] * 0.5 * width, x[2] * 0.5 * height)
    r = R.dot(p)
    n = [float(cx + r[0]), float(cy + r[1]),float(cz + r[2])]
    points.InsertNextPoint(n) # 0 ... 7
  
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
def create_3d_frame(coords, s1, s2, s3, angles, dir, scale):

  centers, min, max, elem_dim = coords

  dx = elem_dim[0]
  dy = elem_dim[1] 
  dz = elem_dim[2]

  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0

  for i in range(len(s1)):
    coord = centers[i]
    angle = angles[i] if angles is not None else None
    if dir == 'horizontal' or dir == 'all':
      create_centered_bar(cells, points, coord, (scale * dx, scale * s1[i] * dy, scale * s1[i] * dz), angle)
    if dir == 'vertical' or dir == 'all':
      create_centered_bar(cells, points, coord, (scale * s2[i] * dx, scale * dy, scale * s2[i] * dz), angle)
    if dir == 'sagittal' or dir == 'all':
      create_centered_bar(cells, points, coord, (scale * s3[i] * dx, scale * s3[i] * dy, scale * dz), angle)
    
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)
    
  return polydata


## for the robot arm we have check for the two nondesign holes as they are within the
## convex hull of the design :(
def valid_position(pos, coords):
  mi, ma = coords[1:3]
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  if int(delta[0]) == 508 and int(delta[2]) == 126:
    if (pos[0] + 147.0)**2 + pos[2]**2 < 40.0**2: # center -147, 0, 0
      return False 
    if (pos[0] - 250.0)**2 + pos[2]**2 < 37.0**2: # center 250, 0, 0
      return False 
  

  return True   

## without rotation and shearing
def create_3d_frame_ip(coords, s1, s2, s3, angles, ip_nx, grad, dir, scale):
  centers, min, max = coords[0:3] # we cannot use the first region element element dimensions 
  
  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0

  ip_data, ip_near, out, ndim = get_interpolation(coords, grad, s1, s2, s3, ip_nx, angles)

  dx = (max[0] - min[0]) / ip_nx

  within = 0
  invalid = 0
  for i in range(len(out)):
    coord = out[i]
    s1, s2, s3 = ip_data[i][0:3]
    angle = None if angles is None else ip_data[i][3:6]
    
    if s1 > 0.0:
      within += 1
      
      if not valid_position(coord, coords):
        invalid += 1
        continue
      
      if dir == 'horizontal' or dir == 'all':
        create_centered_bar(cells, points, coord, (scale * dx, scale * s1 * dx, scale * s1 * dx), angle)
      if dir == 'vertical' or dir == 'all':
        create_centered_bar(cells, points, coord, (scale * s2 * dx, scale * dx, scale * s2 * dx), angle)
      if dir == 'sagittal' or dir == 'all':
        create_centered_bar(cells, points, coord, (scale * s3 * dx, scale * s3 * dx, scale * dx), angle)
    
  if grad <> 'nearest':  
    print str(within) + ' elements out of ' + str(len(out)) + ' in convex hull'  
  if invalid > 0:  
    print str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)'
    
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)
    
  return polydata

# this is copy & paste from matviz_2d but extendet to 3D
# @param nx_ip number of interpolations within x
def get_interpolation(coords, grad, s1, s2, s3, nx, angle = None):
  # we make our own elem
  centers, mi, ma = coords[0:3] # skip elem
 
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  # where we want nodes
  ny = int(delta[1] / (delta[0] / nx))
  nz = int(delta[2] / (delta[0] / nx))
 
  if ny == 0 or nz == 0 or nx == 0:
    print 'chose a higher hom_samples such that also the smallest side gets discretized'
    exit()

  out = numpy.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3))
  idx = 0
  for z in range(nz+1):
    for y in range(ny+1):
      for x in range(nx+1):
        out[idx] = ((mi[0] + float(x)/nx * delta[0], mi[1] + float(y)/ny * delta[1], mi[2] + float(z)/nz * delta[2]))
        idx += 1

  v = numpy.zeros((len(s1),  3 if angle == None else 6))
  v[:,0] = s1[:,0]
  v[:,1] = s2[:,0]
  v[:,2] = s3[:,0]
  if angle <> None:
    v[:,3:6] = angle[:,:]
  
  ip_data = ip.griddata(centers, v, out, grad, -1.0)
  # any interpolation but nearest neighbor can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad <> 'nearest' else None
  
  return ip_data, ip_near, out, (nx, ny, nz) 


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
