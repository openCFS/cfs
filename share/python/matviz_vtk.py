import vtk
from vtk.util.colors import *
from numpy import *
from matviz_rot import *
import scipy.interpolate as ip 

# # creates 3D data to vtkPolyData
def create_vtk_poly_data(angle, data):
  
  assert(len(angle) == len(data))
  samples = int(sqrt(len(angle)))
  assert(samples ** 2 == len(angle))
  
  pts = vtk.vtkPoints()
  
  for i in range(len(angle)):
    radius = data[i]
    
    # we use not the usual conversion for the spherical coordinate system but the first row of the rotation matrix!
    p = radius * to_vector(angle[i])
    
    # print "phi="  + str(phi) + " theta=" + str(theta) + " radius=" + str(radius) + " -> x=" + str(x) + " y=" + str(y) + " z=" + str(z)
    pts.InsertNextPoint(p[0], p[1], p[2])

  pdo = vtk.vtkPolyData()
  pdo.SetPoints(pts)
  # Allocate memory for triangles
  pdo.Allocate((samples - 1) * (samples - 1) * 2)
  
  # Generate two triangles per quad
  for i in range(0, samples - 1):
    for j in range(0, samples - 1):
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, i * samples + j)
      tri.GetPointIds().SetId(1, i * samples + j + 1)
      tri.GetPointIds().SetId(2, (i + 1) * samples + j)
      # print "i=" + str(i) + " j=" + str(j) + " ids=" + str(tri.GetPointIds())
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
      
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, (i + 1) * samples + j + 1)
      tri.GetPointIds().SetId(1, (i + 1) * samples + j)
      tri.GetPointIds().SetId(2, i * samples + j + 1)
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
  
  return pdo

# # create a list of vtk actors displaying symmetry planesfrom vtk.util.colors import
# @param minim
def create_symmety_planes(minima, scale, add_planes):
  # code source: http://code.google.com/p/pythonxy/source/browse/src/python/vtk/DOC/Examples/Rendering/Cylinder.py
 
  minima = []
  minima.append(((0, 0), 1))
  minima.append(((0, numpy.pi / 2), 1))
  minima.append(((numpy.pi / 2, numpy.pi / 2), 1))
 
  actors = [] 
   
  for i in range(len(minima)):
     
    # Create cylinder
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(scale / 200)
    cylinder.SetCenter(0, 0, 0)
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
    phi = angle[0]
    theta = angle[1]

    
    print "angle: " + str(angle) + " -> " + str(to_vector(angle)) + " = " + str(minima[i][1])
    
    # actor_c.RotateX(90)
    # actor_c.RotateY(angle[0] * 180/numpy.pi)
    # actor_c.RotateZ(angle[1] * 180/numpy.pi)


    # cylinders!!!
    actor_c.RotateZ(90)
    actor_c.RotateY(angle[0] * 180 / numpy.pi)
    actor_c.RotateX(angle[1] * 180 / numpy.pi)


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

      # actor_d.RotateX(angle[0] * 180/numpy.pi)
      # actor_d.RotateY(angle[1] * 180/numpy.pi)
      # actor_d.RotateZ(0)

      # plane!
      actor_d.RotateY(90)
      actor_d.RotateZ(90)
      actor_d.RotateY(angle[0] * 180 / numpy.pi)
      actor_d.RotateX(angle[1] * 180 / numpy.pi)
      

      
      actors.append(actor_d)
    
  return actors


# show the data on the screen
# @planes list of vtk actors containing symmetry planes 
def show_vtk(polydata, res, planes=[],show_edges=False):
  # Create a mapper and actor
  mapper = vtk.vtkPolyDataMapper()
  if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(polydata)
  else:
    mapper.SetInputData(polydata)
  
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(0.5, 0.5, 0.5)  # (R,G,B)
  
  if show_edges: # show surface with edges
    actor.GetProperty().EdgeVisibilityOn()
    
  # Setup a renderer, render window, and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  # renderWindow.SetWindowName("Test")
   
  renderWindow.AddRenderer(renderer);
  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)
   
  renderWindow.SetSize(res, res)
   
  # Add the actor to the scene
  for i in range(len(planes)):
    renderer.AddActor(planes[i])

  renderer.AddActor(actor)
  
  renderer.SetBackground(1, 1, 1)  # Background color white
   

  # Render and interact
  renderWindow.Render()
   
  #*** SetWindowName after renderWindow.Render() is called***
  # renderWindow.SetWindowName("Test")
   
  renderWindowInteractor.Start()

def create_point_vector_centered_bar(center, dim, angle=None):
  
  # print "ccb: c=" + str(center) + " l=" + str(length) + " t=" + str(thick) + " dir=" + dir  
  length, width, height = dim
  angle_x, angle_y, angle_z = (0.0, 0.0, 0.0) if angle is None else angle
  
  
  cthetax = cos(angle_x)
  sthetax = sin(angle_x)
  cthetay = cos(angle_y)
  sthetay = sin(angle_y)
  cthetaz = cos(angle_z)
  sthetaz = sin(angle_z)
  
  # see DesignMaterial::SetRotationMatrix()
  R = zeros((3, 3))
  R[0][0] = cthetay * cthetaz
  R[0][1] = -cthetay * sthetaz
  R[1][0] = cthetax * sthetaz + sthetax * sthetay * cthetaz
  R[1][1] = cthetax * cthetaz - sthetax * sthetay * sthetaz
  R[0][2] = sthetay
  R[1][2] = -sthetax * cthetay
  R[2][0] = sthetax * sthetaz - cthetax * sthetay * cthetaz
  R[2][1] = sthetax * cthetaz + cthetax * sthetay * sthetaz
  R[2][2] = cthetax * cthetay

  cx, cy, cz = center
  # points = vtk.vtkPoints()
  points = []
  for x in [(-1.0, -1.0, -1.0), (1.0, -1.0, -1.0), (1.0, 1.0, -1.0), (-1.0, 1.0, -1.0), (1.0, -1.0, 1.0), (-1.0, -1.0, 1.0), (1.0, 1.0, 1.0), (-1.0, 1.0, 1.0)]:
    p = (x[0] * 0.5 * length, x[1] * 0.5 * width, x[2] * 0.5 * height)
    r = R.dot(p)
    n = [float(cx + r[0]), float(cy + r[1]), float(cz + r[2])]
    points.append(n)  # 0 ... 7
  return points

def create_centered_bars(cells, points, coords, dim, angle=None, not_drawn = None):
  # helper for create_frame
  # @param cells  vtk.vtkCellArray() where cells are added via InsertNextCell
  # @param points vtk.vtkPoints() where the points are added
  # @param dim list of length, width, height
  # 
  # @param cells, points: cells and points array with all current VTK elements
  # @param center: center of current cell
  # @param dim: (length, width, height) of the current cell
  
  #optional parameters: @param  not_drawn (faces which are not drawn)
  #                     @param angle list of angle_x, angle_y, angle_z or None
  for i in range(len(coords)):
    create_centered_bar(cells, points, coords[i],dim, angle,not_drawn)

def create_centered_bar(cells, points, center, dim, angle=None,not_drawn = None):
  # helper for create_cross and create_frame
  # @param cells  vtk.vtkCellArray() where cells are added via InsertNextCell
  # @param points vtk.vtkPoints() where the points are added
  # @param dim list of length, width, height
  # 
  # @param cells, points: cells and points array with all current VTK elements
  # @param center: center of current cell
  # @param dim: (length, width, height) of the current cell
  
  #optional parameters: @param  not_drawn (faces which are not drawn)
  #                     @param angle list of angle_x, angle_y, angle_z or None
  
  base = points.GetNumberOfPoints()
  # calculate corner points of quad and add them to global points list
  point_vec = create_point_vector_centered_bar(center, dim, angle)
  for i in range(len(point_vec)):
    points.InsertNextPoint(point_vec[i])  # 0 ... 7
  
  # Create a cell array to store the quad in
  # quads = vtk.vtkCellArray()
  
  # cell is created out of different VTK quads
  # Create a quad on the four points
  if not_drawn is None or not_drawn[0] != 'front' or not_drawn[1] != 'front':
    # front face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 0)
    quad.GetPointIds().SetId(1, base + 1)
    quad.GetPointIds().SetId(2, base + 2)
    quad.GetPointIds().SetId(3, base + 3)
    cells.InsertNextCell(quad)
  
  if not_drawn is None or not_drawn[0] != 'right' or not_drawn[1] != 'right':
    # right face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 1)
    quad.GetPointIds().SetId(1, base + 4)
    quad.GetPointIds().SetId(2, base + 6)
    quad.GetPointIds().SetId(3, base + 2)
    cells.InsertNextCell(quad)
  
  if not_drawn is None or not_drawn[0] != 'back' or not_drawn[1] != 'back':
    # back face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 4)
    quad.GetPointIds().SetId(1, base + 6)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 5)
    cells.InsertNextCell(quad)
  
  if not_drawn is None or not_drawn[0] != 'top' or not_drawn[1] != 'top':
    # top face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 2)
    quad.GetPointIds().SetId(1, base + 6)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 3)
    cells.InsertNextCell(quad)
  
  if not_drawn is None or not_drawn[0] != 'left' or not_drawn[1] != 'left':
    # left face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 0)
    quad.GetPointIds().SetId(1, base + 5)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 3)
    cells.InsertNextCell(quad)
  
  if not_drawn is None or not_drawn[0] != 'bottom' or not_drawn[1] != 'bottom':
    # bottom face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 0)
    quad.GetPointIds().SetId(1, base + 1)
    quad.GetPointIds().SetId(2, base + 4)
    quad.GetPointIds().SetId(3, base + 5)
    cells.InsertNextCell(quad)

# # without rotation and shearing
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

def sign(p1,p2,p3):
  return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);

def point_in_triangle(pt, v1,v2,v3):
  b1 = sign(pt, v1, v2) < 0.;
  b2 = sign(pt, v2, v3) < 0.;
  b3 = sign(pt, v3, v1) < 0.;
  return ((b1 == b2) and (b2 == b3));

def test_point_outside_circle(midpoint,radius, point):
  if (point[0] - midpoint[0]) ** 2 + (point[1] - midpoint[1]) ** 2 < radius ** 2:
    return False
  else:
    return True    

    
def valid_bar_position_apod6(coords,center, dim, angle=None):
  point_vec = create_point_vector_centered_bar(center, dim, angle)
  valid = True
  for i in range(len(point_vec)):
    if not valid_position_apod6(point_vec[i], coords):
      valid = False
  return valid

# # for the robot arm we have check for the two nondesign holes as they are within the
# # convex hull of the design :(
def valid_position_robot(pos, coords,opt=0):
  #mi, ma = coords[1:3]
  #delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
 # if int(delta[0]) == 508 and int(delta[2]) == 126:
  #if (pos[0] + 147.4) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center -147, 0, 0
  if (pos[0] + 147.4) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center -147, 0, 0
    return False 
  if (pos[0] - 250.0) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center 250, 0, 0
    return False
  return True

def valid_ring_position_apod6(pos, coords,opt = 0. ):
  # option opt: change cut out area for validation mesh
  # coordinates of the holes manually, returns False if point is inside a hole
  # mesh is rotated by Ry
  #ay = -0.084636333418591
  ay = 0.
  Ry = numpy.matrix(((math.cos(ay), 0., math.sin(ay)), (0., 1., 0.), (-math.sin(ay), 0., math.cos(ay))))
  tmp = Ry*numpy.matrix(((33.052), (-0.353), (-2.474))).T
  m1 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.046), (-0.353), (-2.518))).T
  m2 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.131), (-0.353), (-2.449))).T
  m3 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.124), (-0.353), (-2.498))).T
  m4 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((32.978), (-0.353), (-2.436))).T
  m5 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((32.971), (-0.353), (-2.485))).T
  m6 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.023), (-0.353), (-2.559))).T
  m7 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.004), (-0.353), (-2.443))).T
  m8 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.096), (-0.353), (-2.450))).T
  m9 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.004), (-0.353), (-2.468))).T
  m10 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.066), (-0.353), (-2.495))).T
  m11 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.096), (-0.353), (-2.475))).T
  m12 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.042), (-0.353), (-2.548))).T
  m13 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  
  #d: thickness of the ring
  d1 = 0.0045
  d2 = 0.0053
  d3 = 0.002#0.00842
  d4 = 0.001#0.01289
  
  r1 = 0.01501  + d1 #0.0158
  r2 = 0.01201  + d2 #0.0128
  r3 = 0.00501  + d3#0.0057
  r4 = 0.002401 + d4#0.0028
  if (pos[0] - m1[0]) ** 2 + (pos[2] - m1[2]) ** 2 < r1 ** 2:
    return True 
  elif (pos[0] - m2[0]) ** 2 + (pos[2] - m2[2]) ** 2 < r2 ** 2:
    return True
  elif (pos[0] - m3[0]) ** 2 + (pos[2] - m3[2]) ** 2 < r3 ** 2:
    return True 
  elif (pos[0] - m4[0]) ** 2 + (pos[2] - m4[2]) ** 2 < r3 ** 2:
    return True
  elif (pos[0] - m5[0]) ** 2 + (pos[2] - m5[2]) ** 2 < r3 ** 2:
    return True
  elif (pos[0] - m6[0]) ** 2 + (pos[2] - m6[2]) ** 2 < r3 ** 2:
    return True
  elif (pos[0] - m7[0]) ** 2 + (pos[2] - m7[2]) ** 2 < r3 ** 2:
    return True
  elif (pos[0] - m8[0]) ** 2 + (pos[2] - m8[2]) ** 2 < r1 ** 2:
    return True
  elif (pos[0] - m9[0]) ** 2 + (pos[2] - m9[2]) ** 2 < r1 ** 2:
    return True
  elif (pos[0] - m10[0]) ** 2 + (pos[2] - m10[2]) ** 2 < r4 ** 2:
    return True
  elif (pos[0] - m11[0]) ** 2 + (pos[2] - m11[2]) ** 2 < r4 ** 2:
    return True
  elif (pos[0] - m12[0]) ** 2 + (pos[2] - m12[2]) ** 2 < r4 ** 2:
    return True
  elif (pos[0] - m13[0]) ** 2 + (pos[2] - m13[2]) ** 2 < r3 ** 2:
    return True


# # for the apod6 part we have check for the holes in nondesign region as they are within the
# # convex hull of the design :(
def valid_position_apod6(pos, coords,opt = 0. ):
  # option opt: change cut out area for validation mesh
  # coordinates of the holes manually, returns False if point is inside a hole
  # mesh is rotated by Ry
  #ay = -0.084636333418591
  ay = 0.
  Ry = numpy.matrix(((math.cos(ay), 0., math.sin(ay)), (0., 1., 0.), (-math.sin(ay), 0., math.cos(ay))))
  tmp = Ry*numpy.matrix(((33.052), (-0.353), (-2.474))).T
  m1 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.046), (-0.353), (-2.518))).T
  m2 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.131), (-0.353), (-2.449))).T
  m3 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.124), (-0.353), (-2.498))).T
  m4 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((32.978), (-0.353), (-2.436))).T
  m5 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((32.971), (-0.353), (-2.485))).T
  m6 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.023), (-0.353), (-2.559))).T
  m7 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.004), (-0.353), (-2.443))).T
  m8 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.096), (-0.353), (-2.450))).T
  m9 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.004), (-0.353), (-2.468))).T
  m10 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.066), (-0.353), (-2.495))).T
  m11 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.096), (-0.353), (-2.475))).T
  m12 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  tmp = Ry*numpy.matrix(((33.042), (-0.353), (-2.548))).T
  m13 = [tmp[0][0],tmp[1][0],tmp[2][0]]
  
  r1 = 0.01501
  r2 = 0.01201
  r3 = 0.00501
  r4 = 0.002401
  if (pos[0] - m1[0]) ** 2 + (pos[2] - m1[2]) ** 2 < r1 ** 2:
    return False 
  elif (pos[0] - m2[0]) ** 2 + (pos[2] - m2[2]) ** 2 < r2 ** 2:
    return False
  elif (pos[0] - m3[0]) ** 2 + (pos[2] - m3[2]) ** 2 < r3 ** 2:
    return False 
  elif (pos[0] - m4[0]) ** 2 + (pos[2] - m4[2]) ** 2 < r3 ** 2:
    return False
  elif (pos[0] - m5[0]) ** 2 + (pos[2] - m5[2]) ** 2 < r3 ** 2:
    return False
  elif (pos[0] - m6[0]) ** 2 + (pos[2] - m6[2]) ** 2 < r3 ** 2:
    return False
  elif (pos[0] - m7[0]) ** 2 + (pos[2] - m7[2]) ** 2 < r3 ** 2:
    return False
  elif (pos[0] - m8[0]) ** 2 + (pos[2] - m8[2]) ** 2 < r1 ** 2:
    return False
  elif (pos[0] - m9[0]) ** 2 + (pos[2] - m9[2]) ** 2 < r1 ** 2:
    return False
  elif (pos[0] - m10[0]) ** 2 + (pos[2] - m10[2]) ** 2 < r4 ** 2:
    return False
  elif (pos[0] - m11[0]) ** 2 + (pos[2] - m11[2]) ** 2 < r4 ** 2:
    return False
  elif (pos[0] - m12[0]) ** 2 + (pos[2] - m12[2]) ** 2 < r4 ** 2:
    return False
  elif (pos[0] - m13[0]) ** 2 + (pos[2] - m13[2]) ** 2 < r3 ** 2:
    return False
  
  # big right triangle
  tmp = [Ry*numpy.matrix(((33.070), (-0.353), (-2.563))).T, Ry*numpy.matrix(((33.078),(-0.353),(-2.503))).T, Ry*numpy.matrix(((33.127+opt), (-0.353), (-2.508))).T]
  corners_b = [[tmp[0][0][0], tmp[0][2][0]], [tmp[1][0][0], tmp[1][2][0]], [tmp[2][0][0], tmp[2][2][0]]]                         
  # small right triangle
  tmp = [Ry*numpy.matrix(((33.077), (-0.353), (-2.513))).T, Ry*numpy.matrix(((33.078),(-0.353),(-2.503))).T, Ry*numpy.matrix(((33.087), (-0.353), (-2.504))).T]
  corners_s = [[tmp[0][0][0], tmp[0][2][0]], [tmp[1][0][0], tmp[1][2][0]], [tmp[2][0][0], tmp[2][2][0]]]
  #corners_s = [[33077., -2513.], [33078.,-2503.], [33087., -2504.]]
  # small circle
  tmp = Ry*numpy.matrix(((33.087), (-0.353),(-2.515))).T
  mid_s = [tmp[0][0],tmp[2][0]]
  r_s = 0.010
  if point_in_triangle([pos[0],pos[2]], corners_b[0],corners_b[1],corners_b[2]):
    if (test_point_outside_circle(mid_s,r_s, [pos[0],pos[2]])) and point_in_triangle([pos[0],pos[2]], corners_s[0],corners_s[1],corners_s[2]):
      return True
    else:
      return False
   # big left triangle
  tmp = [Ry*numpy.matrix(((32.964),(-0.353), (-2.494))).T, Ry*numpy.matrix(((33.018),(-0.353),(-2.500))).T, Ry*numpy.matrix(((33.010),(-0.353), (-2.559))).T]
  corners_b = [[tmp[0][0][0], tmp[0][2][0]], [tmp[1][0][0], tmp[1][2][0]], [tmp[2][0][0], tmp[2][2][0]]]
  # small right triangle
  tmp = [Ry*numpy.matrix(((33.008),(-0.353), (-2.498))).T, Ry*numpy.matrix(((33.018),(-0.353),(-2.500))).T, Ry*numpy.matrix(((33.017), (-0.353), (-2.509))).T]
  corners_s = [[tmp[0][0][0], tmp[0][2][0]], [tmp[1][0][0], tmp[1][2][0]], [tmp[2][0][0], tmp[2][2][0]]]
  # small circle
  tmp = Ry*numpy.matrix(((33.007),(-0.353),(-2.508))).T
  mid_s = [tmp[0][0],tmp[2][0]]
  r_s = 0.010
  if point_in_triangle([pos[0],pos[2]], corners_b[0],corners_b[1],corners_b[2]):
    if (test_point_outside_circle(mid_s,r_s, [pos[0],pos[2]])) and point_in_triangle([pos[0],pos[2]], corners_s[0],corners_s[1],corners_s[2]):
      return True
    else:
      return False
  return True
  
# # without rotation and shearing
def create_3d_frame_ip(coords, s1, s2, s3, angles, ip, grad, scale, valid_position, thres=0.0, csize = None):
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # NOT tested with angles
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # valid_position: returns false if point inside the convex hull of the part, should be excluded, otherwise true.
  #                 Needs to be implemented for every mechanical part, currently available for robot and apod6.
  #                 If part is not implemented valid_position is None and no cells inside the convex hull are removed from the structure
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]
  
  # point coordinates from h5 file
  centers, min, max = coords[0:3] 
  
  # create vtk cells and points
  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0
  
  # set size dx/dy/dz of one cell
  if csize is None:
    dx = (max[0] - min[0]) / ip[0]
    dy = (max[1] - min[1]) / ip[1]
    dz = (max[2] - min[2]) / ip[2]
  else:
    dx = csize[0]
    dy = csize[1]
    dz = csize[2]
    
  # calculate 3d volume of the structure
  vol = calc_cross_elem_vol_3D(s1,s2,s3)
  
  # calculate interpolated values of the design variables s1,s2,s3 for a uniform 3d grid 
  ip_data, ip_near, out, ndim, scale_ = get_interpolation(coords, grad, s1, s2, s3, dx, dy, dz, angles)

  #scales the lattice cells to fit in the design domain exactly
  #scale = scale_max 
  
  # counters for visualized or non-visualized cells inside the convex hull of the part
  within = 0
  invalid = 0
  real_volume = 0.
  for i in range(len(out)):
    coord = out[i]
    # get interpolated design variables on uniform grid
    s1, s2, s3 = ip_data[i][0:3]
    angle = None if angles is None else ip_data[i][3:6]
    # if s1 < 0 point is out of the convex hull of the part
    if s1 > 0.0:
      if not valid_position is None and not valid_position(coord, coords):
        invalid += 1
        continue
      within += 1
      if s1 >= thres or s2 >= thres or s3 >= thres:
        # draw each bar of 3D cross for s1 > s2,s3
        if s1 >= s2 and s1 >= s3:
          if s1 >= thres:#valid_bar_position_apod6(points,coord, (scale * scale_[0] * dx, scale * s1 * dx, scale * s1 * dx), angle):
            # draw thickest bars first 
            coords = []
            for i  in range(4):
              coords.append(coord)
            # add offset for s1 for all s1-bars
            coords[0] += [0.,-s1/4.,s1/4.]
            coords[1] += [0.,s1/4.,s1/4.]
            coords[2] += [0.,-s1/4.,-s1/4.]
            coords[3] += [0.,s1/4.,-s1/4.]        
            create_centered_bars(cells, points, coords, (scale * scale_[0] * dx, scale * 0.5 * s1 * dx, scale * 0.5 * s1 * dx), angle,['right','left'])
            real_volume += scale * scale_[0] * dx * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,scale* dx * s1 * 0.5 + scale * 0.25 * (scale_[1] * dy - dx * s1),0.]
          dy_offset = scale * 0.5 * (scale_[1]*dy-s1*dx)
          #add two parts of s2-bar, two parts are necessary that it doesn't intersect the s1-bar
          if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            coords = []
            for i  in range(4):
              coords.append(coord)
            # add offset for s2 for all s2-bars
            coords[0] += [-s1/4.,0., s1/4.]
            coords[1] += [s1/4., 0., s1/4.]
            coords[2] += [-s1/4.,0., -s1/4.]
            coords[3] += [s1/4., 0.,  -s1/4.]  
            create_centered_bars(cells, points, coords, (scale * 0.5 * s2 * dy, dy_offset, scale * 0.5 * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy* dy_offset* scale * s2 * dy
          
          coord_offset = [0.,0.,scale * dx * s1 * 0.5 + scale * 0.25 * (scale_[2] * dz - dx * s1)]
          dz_offset = scale * 0.5 * (scale_[2]*dz-s1*dx)
          #add two parts of s3-bar, two parts are necessary that it doesn't intersect the s1-bar
          if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            coords = []
            for i  in range(4):
              coords.append(coord)
            # add offset for s3 for all s3-bars
            coords[0] += [-s1/4.,s1/4.,0.]
            coords[1] += [s1/4., s1/4.,0.]
            coords[2] += [-s1/4., -s1/4., 0.]
            coords[3] += [s1/4., -s1/4.,0.]   
            create_centered_bars(cells, points, coords, (scale * 0.5 * s3 * dz, scale * 0.5 * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset

        # draw each bar of 3D cross for s2 > s1,s3
        elif s2 >= s1 and s2 >= s3:
          if s2 >= thres:#valid_bar_position_apod6(points,coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy * scale * scale_[1]* dy * scale * s2 * dy
          coord_offset = [scale* dy * s2 * 0.5 + scale * 0.25 * (scale_[0] * dx - dy * s2),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx-s2*dy)
          if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,0.,scale * dy * s2 * 0.5 + scale * 0.25 * (scale_[2] * dz - dy * s2)]
          dz_offset = scale * 0.5 * (scale_[2] * dz-s2*dy)
          if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
          if s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
        # draw each bar of 3D cross for s3 > s1,s2
        elif s3 >= s1 and s3 >= s2:  
          if s3 >= thres:#valid_bar_position_apod6(points,coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle):
            create_centered_bar(cells, points, coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * scale * scale_[2] * dz
          coord_offset = [scale* dz * s3 * 0.5 + scale * 0.25 * (scale_[0] * dx - dz * s3),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx-s3*dz)
          if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,scale * dz * s3 * 0.5 + scale * 0.25 * (scale_[1] * dy - dz * s3),0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy-s3*dz)
          if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy * dy_offset * scale * s2 * dy
          if s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy * dy_offset * scale * s2 * dy

  real_volume /= within * dx * dy * dz
  print 'volume of 3D Two-scale result = ' + str(vol)
  print 'real volume of 3D lattice = ' + str(real_volume)  
  if grad <> 'nearest':  
    print str(within) + ' elements out of ' + str(len(out)) + ' in convex hull'  
  if invalid > 0:  
    print str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)'
    
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)
    
  return polydata

# # without rotation and shearing
def create_3d_cross_ip(coords, s1, s2, s3, angles, ip_nx, grad, scale,valid_position,thres=0.0,csize = None):
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # valid_position: returns false if point inside the convex hull of the part, should be excluded, otherwise true.
  #                 Needs to be implemented for every mechanical part, currently available for robot and apod6.
  #                 If part is not implemented valid_position is None and no cells inside the convex hull are removed from the structure
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]
  
  # point coordinates from h5 file
  centers, min, max = coords[0:3] 
  
  # create vtk cells and points
  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0
  
  # set size dx/dy/dz of one cell
  if csize is None:
    dx = (max[0] - min[0]) / ip_nx
    dy = dx
    dz = dx
  else:
    dx = csize[0]
    dy = csize[1]
    dz = csize[2]
    
  # calculate 3d volume of the structure
  vol = calc_cross_elem_vol_3D(s1,s2,s3)
  
  # calculate interpolated values of the design variables s1,s2,s3 for a uniform 3d grid 
  ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1, s2, s3, dx,dy,dz, angles)

  #scales the lattice cells to fit in the design domain exactly
  #scale = scale_max 
  
  # counters for visualized or non-visualized cells inside the convex hull of the part
  within = 0
  invalid = 0
  real_volume = 0.
  for i in range(len(out)):
    coord = out[i]
    # get interpolated design variables on uniform grid
    s1, s2, s3 = ip_data[i][0:3]
    angle = None if angles is None else ip_data[i][3:6]
    # if s1 < 0 point is out of the convex hull of the part
    if s1 > 0.0:
      if not valid_position is None and not valid_position(coord, coords):
        invalid += 1
        continue
      within += 1
      if s1 >= thres or s2 >= thres or s3 >= thres:
        # draw each bar of 3D cross for s1 > s2,s3
        if s1 >= s2 and s1 >= s3:
          if s1 >= thres:#valid_bar_position_apod6(points,coord, (scale * scale_[0] * dx, scale * s1 * dx, scale * s1 * dx), angle):
            # draw thickest bar first      
            create_centered_bar(cells, points, coord, (scale * scale_[0] * dx, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += scale * scale_[0] * dx * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., scale * s1 * dy * 0.5 + scale * 0.25 * (scale_[1] * dy - s1 * dy), 0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy - s1 * dy)
          #add two parts of s2-bar, two parts are necessary that it doesn't intersect the s1-bar
          if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * dy_offset * scale * s2 * dz
          if s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * dy_offset * scale * s2 * dz
          
          coord_offset = [0., 0., scale * s1 * dz * 0.5 + scale * 0.25 * (scale_[2] * dz - s1 * dz)]
          dz_offset = scale * 0.5 * (scale_[2] * dz - s1 * dz)
          #add two parts of s3-bar, two parts are necessary that it doesn't intersect the s1-bar
          if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
          if s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
        # draw each bar of 3D cross for s2 > s1,s3
        elif s2 >= s1 and s2 >= s3:
          if s2 >= thres:#valid_bar_position_apod6(points,coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord, (scale * s2 * dx, scale * scale_[1] * dy, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * scale * scale_[1] * dy * scale * s2 * dz
          coord_offset = [scale * s2 * dx * 0.5 + scale * 0.25 * (scale_[0] * dx - s2 * dx), 0., 0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx - s2 * dx)
          if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., 0., scale * s2 * dz * 0.5 + scale * 0.25 * (scale_[2] * dz - s2 * dz)]
          dz_offset = scale * 0.5 * (scale_[2] * dz - s2 * dz)
          if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
          if s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
        # draw each bar of 3D cross for s3 > s1,s2
        elif s3 >= s1 and s3 >= s2:
          if s3 >= thres:#valid_bar_position_apod6(points,coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle):
            create_centered_bar(cells, points, coord, (scale * s3 * dx, scale * s3 * dy, scale * scale_[2] * dz), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * scale * scale_[2] * dz
          coord_offset = [scale * s3 * dx * 0.5 + scale * 0.25 * (scale_[0] * dx - s3 * dx),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx - s3 * dx)
          if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., scale * dy * s3 * 0.5 + scale * 0.25 * (scale_[1] * dy - s3 * dy), 0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy - s3 * dy)
          if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dy * dy_offset * scale * s2 * dy
          if s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * dy_offset * scale * s2 * dz

  real_volume /= within * dx * dy * dz
  print 'volume of 3D Two-scale result = ' + str(vol)
  print 'real volume of 3D lattice = ' + str(real_volume)  
  if grad <> 'nearest':  
    print str(within) + ' elements out of ' + str(len(out)) + ' in convex hull'  
  if invalid > 0:  
    print str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)'
    
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)
    
  return polydata

# this is copy & paste from matviz_2d but extended to 3D
# @param nx_ip number of interpolations within x
def get_interpolation(coords, grad, s1, s2, s3, dx, dy, dz, angle=None):
  # we make our own regular element grid
  centers, mi, ma = coords[0:3]  # skip elem
 
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  # where we want nodes
  nx = int(delta[0] / dx)
  ny = int(delta[1] / dy)
  nz = int(delta[2] / dz)
  
  scale_x = delta[0]/(nx*dx)
  scale_y = delta[1]/(ny*dy)
  scale_z = delta[2]/(nz*dz)
  
  print "delta: " + str(delta)
  print "dx,dy,dz: " + str(dx) + ", "+ str(dy) + ", " + str(dz) 
  if ny == 0 or nz == 0 or nx == 0:
    print 'chose a higher hom_samples such that also the smallest side gets discretized'
    exit()

  out = numpy.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3))
  idx = 0
  for z in range(nz + 1):
    for y in range(ny + 1):
      for x in range(nx + 1):
        out[idx] = ((mi[0] + 0.5*dx + float(x) / nx * delta[0], mi[1] + 0.5*dy +  float(y) / ny * delta[1], mi[2] + 0.5*dx + float(z) / nz * delta[2]))
        idx += 1
  if s2 is None and s3 is None:
      v = numpy.zeros((len(s1), 1))
      v[:, 0] = s1[:, 0]
  else:
    v = numpy.zeros((len(s1), 3 if angle is None else 6))
    v[:, 0] = s1[:, 0]
    v[:, 1] = s2[:, 0]
    v[:, 2] = s3[:, 0]
    if angle is not None:
	     v[:, 3:6] = angle[:, :]
  
  ip_data = ip.griddata(centers, v, out, grad, -1.0)
  # any interpolation, ie. linear interpolation can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation which can also interpolate values outside the convex hull
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad <> 'nearest' else None
  
  return ip_data, ip_near, out, (nx, ny, nz), (scale_x, scale_y, scale_z)  


# # litte helper
# @param save filename or none
# @param list which might be empty
def show_write_vtk(poly, res, save, actors=[]):
  if save:
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(poly)
    else:
      writer.SetInputData(poly)
    #writer.SetDataModeToAscii()
    writer.SetFileName(save)
    writer.Write()
    print "saved polydata to file", save
  else:
    show_vtk(poly, res, actors)  
    
def calc_cross_elem_vol_3D(s1,s2,s3):
  # calculates element volume of cross structure in 3D
  vol = 0.
  for i in range(len(s1[:,0])):
    stiff1 = s1[i,0]
    stiff2 = s2[i,0]
    stiff3 = s3[i,0]
    if stiff1 >= stiff2 and stiff1 >= stiff3:
      vol += stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff1*stiff3*stiff3 - stiff1*stiff2*stiff2
    elif stiff2 >= stiff1 and stiff2 >= stiff3:
      vol += stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff2*stiff3*stiff3 - stiff2*stiff1*stiff1
    elif stiff3 >= stiff2 and stiff3 >= stiff2:
      vol += stiff1*stiff1 + stiff2*stiff2 + stiff3*stiff3 - stiff3*stiff2*stiff2 - stiff3*stiff1*stiff1
    else:
      vol += 0.
  return vol/len(s1[:,0])

# writes polydata to file in STL format
# @param save filename or none
def write_stl(polydata,save=None):
  stlWriter =  vtk.vtkSTLWriter()
  fName = save if save else 'surface.stl'
  stlWriter.SetFileName(fName)
  stlWriter.SetInputData(polydata)
  stlWriter.Write()
  
  print "saved polydata to file", fName  
