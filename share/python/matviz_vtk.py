import basecell
import draw_profile_functions
from matviz_rot import *
from numpy import *
import pymp
import scipy.interpolate as ip 
import time
import vtk
import gc
# from mpi_class import * 
try:
  from mpi4py import MPI
  use_mpi = True if MPI.COMM_WORLD.Get_size() > 1 else False
except:
  use_mpi = False

#from vtk.util.colors import *
try:
  import meshpy.triangle as triangle
  from meshpy.tet import MeshInfo, build
except:
  print("Failed to load meshpy - need it for tetrahedralized basecell mesh")

#computation for worker nodes MPI


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

    
    print("angle: " + str(angle) + " -> " + str(to_vector(angle)) + " = " + str(minima[i][1]))
    
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
  print('volume of 3D Two-scale result = ' + str(vol))
  print('real volume of 3D lattice = ' + str(real_volume))  
  if grad != 'nearest':  
    print(str(within) + ' elements out of ' + str(len(out)) + ' in convex hull')  
  if invalid > 0:  
    print(str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)')
    
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
  print('volume of 3D Two-scale result = ' + str(vol))
  print('real volume of 3D lattice = ' + str(real_volume))  
  if grad != 'nearest':  
    print(str(within) + ' elements out of ' + str(len(out)) + ' in convex hull')  
  if invalid > 0:  
    print(str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)')
    
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
  
  if ny == 0 or nz == 0 or nx == 0:
    print('chose a higher hom_samples such that also the smallest side gets discretized')
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
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad != 'nearest' else None
  
  return ip_data, ip_near, out, (nx, ny, nz), (scale_x, scale_y, scale_z)  

def get_interpolation_row_major(coords, grad, s1, s2, s3, dx, dy, dz, angle=None):
  # we make our own regular element grid
  centers, mi, ma = coords[0:3]  # skip elem
 
  delta = (abs(ma[0] + mi[0]), abs(ma[1] + mi[1]), abs(ma[2] + mi[2]))
  # where we want nodes
  nx = int(delta[0] / dx)
  ny = int(delta[1] / dy)
  nz = int(delta[2] / dz)
  
  scale_x = delta[0]/(nx*dx)
  scale_y = delta[1]/(ny*dy)
  scale_z = delta[2]/(nz*dz)
  
  if ny == 0 or nz == 0 or nx == 0:
    print('chose a higher hom_samples such that also the smallest side gets discretized')
    exit()

  out = numpy.zeros(((nx + 1), (ny + 1), (nz + 1), 3))
  idx = 0
  
  for x in range(nx + 1):
    for y in range(ny + 1):
      for z in range(nz + 1):
        out[x,y,z] = [x*dx,y*dy,z*dz]

  assert(s1 is not None and s2 is not None and s3 is not None)      
  v = numpy.zeros((len(s1),3))
  v[:, 0] = s1[:, 0]
  v[:, 1] = s2[:, 0]
  v[:, 2] = s3[:, 0]
  
#   test = [(0.8,0.95,0.025),(0.8,0.05,0.025)]
#   test_data = ip.griddata(centers, v, test, grad, -1.0)
#   print(test_data)
#   sys.exit()
  
  ip_data = ip.griddata(centers, v, out, grad, -1.0)
  # any interpolation, ie. linear interpolation can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation which can also interpolate values outside the convex hull
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad != 'nearest' else None
  
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
    print("saved polydata to file", save)
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
  
  print("saved polydata to file " + fName)  

# similar to create_3d_cross_ip; # without rotation and shearing
# returns
def create_3d_interpretation_ortho(args,coords,s1,s2,s3,scale,samples,grad,thresh):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]
  
  # point coordinates from h5 file
  centers, minimum, maximum = coords[0:3]
  
  print("min:",minimum," max:",maximum)
  
  # appendind cells
  appends = vtk.vtkAppendPolyData()
  
  if scale <= 0:
    scale = 1.0
  
  # assume we always start at (0,0,0)
  # order: min_x,min_y,min_z,max_x,max_y,max_z
  bounds = numpy.ones(6) * (-1)
  bounds[0] = bounds[1] = bounds[2] = 0
  bounds[3] = maximum[0] + minimum[0]
  bounds[4] = maximum[1] + minimum[1]
  bounds[5] = maximum[2] + minimum[2]
  
  print("dimensions:",bounds)
  
  # set size dx/dy/dz of one cell
  dx = bounds[3] / samples[0]
  dy = bounds[4] / samples[1]
  dz = bounds[5] / samples[2]

  min_thresh = 3.0/args.bc_res
  max_thresh = 0.82
    
  delta = (abs(maximum[0] + minimum[0]), abs(maximum[1] + minimum[1]), abs(maximum[2] + minimum[2]))
  # where we want nodes
  nx = int(delta[0] / dx)
  ny = int(delta[1] / dy)
  nz = int(delta[2] / dz)    
  
  print("before ip data")
  start = time.time()
  data_grid, data_grid_near, sample_coords, ndim, scale_ = get_interpolation_row_major(coords, grad, s1, s2, s3, dx, dy, dz)
  end = time.time()
  print("got ip data. elapsed time:",end-start," s")
  
  basecells = pymp.shared.list()
  boundary_flags = pymp.shared.list()
  start = time.time()
  with pymp.Parallel(args.bc_num_threads) as p:
    for id in p.range(nx*ny*nz):
      i, j, k = get_3d_grid_coords(id,nx,ny,nz)
      
      this = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k))
      east = get_interp_3darray_elem(data_grid,data_grid_near,(i+1,j,k))
      top = get_interp_3darray_elem(data_grid,data_grid_near,(i,j+1,k))
      front = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k+1))
      
      assert(this is not None and east is not None and top is not None and front is not None)

      # if one of the values is < min_thresh, set it to min_thresh        
      # if one of the values is > max_thresh, set it to max_thresh
      x1 = min(max(this[0],min_thresh),max_thresh)
      x2 = min(max(east[0],min_thresh),max_thresh)
      y1 = min(max(this[1],min_thresh),max_thresh)
      y2 = min(max(top[1],min_thresh),max_thresh)
      z1 = min(max(this[2],min_thresh),max_thresh)
      z2 = min(max(front[2],min_thresh),max_thresh)
      
      bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta)
      bc_input.eta = 0.7
      bc_input.stiffness_as_diameter = True
      cell_obj = basecell.Basecell(bc_input)
      if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
       #if cell_obj.volume <= args.bc_volume_thresh:
        continue
      
      # translate cell to correct position
      coord  = numpy.asarray(sample_coords[i,j,k])
      
      # flags for meshing circles on the boundary
      flag = [None] * 6
      for c,bound in enumerate(bounds):
        if numpy.isclose(coord[int(c/2)],bound):
          flag[c] = True
      
      with p.lock:
        basecells.append((cell_obj.points,cell_obj.cells,coord))
        boundary_flags.append(flag)
        print("appended ",i,j,k,coord,x1,x2,y1,y2,z1,z2)
  end = time.time()
  print("Openmp Exection TIme",end-start," s")

  start = time.time()
  for i,bc in enumerate(basecells):
    if (i%100==0):
      end_vtk_iteration=time.time()
      print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )

    vtk_points = vtk.vtkPoints()
    for p in bc[0]:
      vtk_points.InsertNextPoint(p)
      
    vtk_cells = vtk.vtkCellArray()
    for ce in bc[1]:
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0,ce[0])
      tri.GetPointIds().SetId(1,ce[1])
      tri.GetPointIds().SetId(2,ce[2])
      vtk_cells.InsertNextCell(tri)
    
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(vtk_points)
    polydata.SetPolys(vtk_cells)
          
    # tranformation of vtk poly data
    transform = vtk.vtkTransform() 
    transform.Translate(bc[2])
    transform.Scale(dx,dy,dz)
       
    # filter to perform transformation
    transformFilter=vtk.vtkTransformPolyDataFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInputData(polydata)
    transformFilter.Update()
    
    # append cells in vtk
    appends.AddInputConnection(transformFilter.GetOutputPort())
    appends.Update() # not sure if we have to do this in each loop iteration
    
    dict = {0:"x_left", 1:"y_left", 2:"z_left", 3:"x_right", 4:"y_right", 5:"z_right"}
    # add meshed boundary circles if necessary
#     for idx,flag in enumerate(boundary_flags[i]):
#       if flag:
#         # returns vtk polydata object
#         pd = mesh_boundary_circle(bc[0], dict[idx], bounds[0:3], bounds[3:6])
#         
# #     appends.AddInputData(pd)
# #     appends.Update() # not sure if we have to do this in each loop iteration
#     
#       # tranformation of vtk poly data
#       transform = vtk.vtkTransform() 
#       transform.Translate(bc[2])
#       transform.Scale(dx,dy,dz)
#          
#       # filter to perform transformation
#       transformFilter=vtk.vtkTransformPolyDataFilter()
#       transformFilter.SetTransform(transform)
#       transformFilter.SetInputData(pd)
#       transformFilter.Update()    
#       appends.AddInputConnection(transformFilter.GetOutputPort())
#       appends.Update() # not sure if we have to do this in each loop iteration
  end=time.time()
  print("Writing to vtk",end-start," s")
  return appends.GetOutput()
    
# @param idx: tuple of three ints storing array indices(i,j,k)
# @param array: return element of array at position idx if exist
# if out of range, return None 
# @param fallback for array, if array elem at idx has value -1    


def create_3d_interpretation_ortho_mpi(args,coords,s1,s2,s3,scale,samples,grad,thresh):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]
  
  # point coordinates from h5 file
  comm = MPI.COMM_WORLD
  print (comm.Get_rank())
  if comm.Get_rank() == 0:

    centers, minimum, maximum = coords[0:3]
    
    print("min:",minimum," max:",maximum)
    
    # appendind cells
    appends = vtk.vtkAppendPolyData()
    
    if scale <= 0:
      scale = 1.0
    
    # assume we always start at (0,0,0)
    # order: min_x,min_y,min_z,max_x,max_y,max_z
    bounds = numpy.ones(6) * (-1)
    bounds[0] = bounds[1] = bounds[2] = 0
    bounds[3] = maximum[0] + minimum[0]
    bounds[4] = maximum[1] + minimum[1]
    bounds[5] = maximum[2] + minimum[2]
    
    print("dimensions:",bounds)
    
    # set size dx/dy/dz of one cell
    dx = bounds[3] / samples[0]
    dy = bounds[4] / samples[1]
    dz = bounds[5] / samples[2]

    min_thresh = 3.0/args.bc_res
    max_thresh = 0.82
      
    delta = (abs(maximum[0] + minimum[0]), abs(maximum[1] + minimum[1]), abs(maximum[2] + minimum[2]))
    # where we want nodes
    nx = int(delta[0] / dx)
    ny = int(delta[1] / dy)
    nz = int(delta[2] / dz)    
    
   
    print("before ip data")
    start = time.time()
    data_grid, data_grid_near, sample_coords, ndim, scale_ = get_interpolation_row_major(coords, grad, s1, s2, s3, dx, dy, dz)
    end = time.time()
    print("got ip data. elapsed time:",end-start," s")

    basecells = list()
    boundary_flags = list()
    args_for_worker =args,data_grid,data_grid_near,sample_coords,ndim,scale,min_thresh,max_thresh,id,nx,ny,nz,bounds 
    chief=Master(args_for_worker)
    # chief.initialize(args_for_worker)
    start=time.time()
    basecells,boundary_flags=chief.run()  
    end=time.time()

    print("MPI BLOCK EXECUTION TIME",end-start," s")

    start=time.time()
    
    for i,bc in enumerate(basecells): 
      
      if (i%200==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
      vtk_points = vtk.vtkPoints()

      for p in bc[0]:
        vtk_points.InsertNextPoint(p)
        
      vtk_cells = vtk.vtkCellArray()
      for ce in bc[1]:
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0,ce[0])
        tri.GetPointIds().SetId(1,ce[1])
        tri.GetPointIds().SetId(2,ce[2])
        vtk_cells.InsertNextCell(tri)
      
      polydata = vtk.vtkPolyData()
      polydata.SetPoints(vtk_points)
      polydata.SetPolys(vtk_cells)
            
      # tranformation of vtk poly data
      transform = vtk.vtkTransform() 
      transform.Translate(bc[2])
      transform.Scale(dx,dy,dz)
         
      # filter to perform transformation
      transformFilter=vtk.vtkTransformPolyDataFilter()
      transformFilter.SetTransform(transform)
      transformFilter.SetInputData(polydata)
      transformFilter.Update()
      
      # append cells in vtk
      appends.AddInputConnection(transformFilter.GetOutputPort())
      appends.Update() # not sure if we have to do this in each loop iteration
      
      dict = {0:"x_left", 1:"y_left", 2:"z_left", 3:"x_right", 4:"y_right", 5:"z_right"}
      # add meshed boundary circles if necessary
  #     for idx,flag in enumerate(boundary_flags[i]):
  #       if flag:
  #         # returns vtk polydata object
  #         pd = mesh_boundary_circle(bc[0], dict[idx], bounds[0:3], bounds[3:6])
  #         
  # #     appends.AddInputData(pd)
  # #     appends.Update() # not sure if we have to do this in each loop iteration
  #     
  #       # tranformation of vtk poly data
  #       transform = vtk.vtkTransform() 
  #       transform.Translate(bc[2])
  #       transform.Scale(dx,dy,dz)
  #          
  #       # filter to perform transformation
  #       transformFilter=vtk.vtkTransformPolyDataFilter()
  #       transformFilter.SetTransform(transform)
  #       transformFilter.SetInputData(pd)
  #       transformFilter.Update()    
  #       appends.AddInputConnection(transformFilter.GetOutputPort())
  #       appends.Update() # not sure if we have to do this in each loop iteration
    end=time.time()
    print("Writing to vtk",end-start," s")
    return appends.GetOutput()

  
  else :
    w=Worker()
    w.run()
  
       # tmp_basecells, tmp_boundary_flags = worker_computaion(args,data_grid,data_grid_near,sample_coords,ndim,scale,min_thresh,max_thresh,idx,nx,ny,nz,bounds)  
       # basecells.extend(tmp_basecells)
       # boundary_flags.extend(tmp_boundary_flags)
   
  
 
    
# @param idx: tuple of three ints storing array indices(i,j,k)
# @param array: return element of array at position idx if exist
# if out of range, return None 
# @param fallback for array, if array elem at idx has value -1      
def get_interp_3darray_elem(array,fallback,idx):
  # must be a ndarray so that indexing with tuples works
  assert(type(array) == numpy.ndarray)
  assert(type(fallback) == numpy.ndarray)
#   assert(array.ndim == 3)
  nx, ny, nz, _ = array.shape
  
  if idx[0] >= nx or idx[0] < 0  or idx[1] >= ny or idx[1] < 0  or idx[2] >= nz  or idx[2] < 0 :
    return None
  else:
    if array[idx][0] == -1 or array[idx][1] == -1 or array[idx][2] == -1:
      return fallback[idx]
    else:
      return array[idx]

# @param id: idx in 1D array (row major)
# @return i,j,k: indices of equivalent 3D array
def get_3d_grid_coords(index,nx,ny,nz):
  i = index % nx
  j = (index - i)/nx % ny
  k = ((index - i)/nx-j)/ny
  
  return int(i),int(j),int(k)

# returns a list with all points lying on given boundary
# @param location, e.g. x_0 for all x=mini[0] and x_1 for all x = maxi[0]
# @param mini/max: list with 3 min max coordinates of domain
def mesh_boundary_circle(points,location,mini,maxi):
  result = []
  major_dir = -1
  for p in points:
    if location == "x_left":
      if isclose(p[0],mini[0]):
        result.append((p[1],p[2]))
        major_dir = 0
    if location == "x_right":
      if isclose(p[0],maxi[0]):
        result.append((p[1],p[2]))
        major_dir = 0
        
    if location == "y_left":
      if isclose(p[1],mini[1]):
        result.append(p)
        major_dir = 1
    if location == "y_right":
      if isclose(p[1],maxi[1]):
        result.append(p)
        major_dir = 1
        
    if location == "z_left":
      if isclose(p[2],mini[2]):
        result.append(p)
        major_dir = 2
    if location == "z_right":
      if isclose(p[2],maxi[2]):            
        result.append(p)
        major_dir = 2
        
  # sort points in circle order
  result.sort(key=lambda c:math.atan2(c[0]-0.5, c[1]-0.5))
  
  info = triangle.MeshInfo()
  test = [ [elem[0],elem[1]] for elem in result]
  info.set_points(test)
  info.set_facets(draw_profile_functions.round_trip_connect(0,len(result)-1))
  mesh = triangle.build(info,generate_faces=True) 
  
  mesh_points = numpy.array(mesh.points)
  mesh_tris = numpy.array(mesh.elements)
  
  minor_dir_1, minor_dir_2 = draw_profile_functions.give_normal_plane_axes(major_dir)
  
  # check whether we're on the left or right side
  comp = mini[major_dir]
  if location.endswith("right"):
    comp = maxi[major_dir]
    
  
  vtk_points = vtk.vtkPoints()
  cells = vtk.vtkCellArray()
  polydata = vtk.vtkPolyData()
  # map from 2d point to 3d point
  for p in mesh_points:
    new_p = numpy.zeros(3)
    new_p[major_dir] = comp
    new_p[minor_dir_1] = p[0]
    new_p[minor_dir_2] = p[1]
    vtk_points.InsertNextPoint(new_p)
  for tri in mesh_tris:
    draw_profile_functions.add_triangle(tri[0], tri[1], tri[2], cells) 
    
  polydata.SetPoints(vtk_points)
  polydata.SetPolys(cells)
             
  return polydata












#MPI IMPLEMENTATION CLASS SHOULD BE MOVED TO SEPRATE FILE BUT UNABLE TO USE FUNCTIONS DEFINED OUTSIDE THE CLASSS :/
def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# these tags will be attached to all messages which are communicated
# between the master and workers
tags = enum('READY', 'DONE', 'EXIT', 'START')

# in this example we use the standard communicator
comm = MPI.COMM_WORLD

# worker class
class Worker():
    # this function is called when creating a worker instance
    def __init__(self):
        # the worker's unique ID
        self.id = comm.Get_rank()
    

    def run(self):
        # a status object containing source, tag and size of a message
        status = MPI.Status()
        # the worker runs forever
        # (until it receives a poison pill, i. e. an 'EXIT' tag)
        while True:
            # the worker tells the master that it is ready to do something.
            # the message data is the worker's ID, the destination is the
            # rank of the master (0).
            comm.send(self.id, dest=0, tag=tags.READY)
            # the worker waits for a message from the master with any tag.
            # the information about source, tag and size of the received
            # message are available from the status object.
            data_from_master = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
            tag = status.Get_tag()
            # if the tag is 'START' run computation
            if tag == tags.START:
                # do some long computation
                result = self.worker_computation(data_from_master)
                # send the result back to the master
                comm.send(result, dest=0, tag=tags.DONE)
            elif tag == tags.EXIT:
                # if the worker receives a message with
                # tag 'EXIT' the loop stops
                del result
                gc.collect()             
                break
        # just before the worker stops running it tells the master that 
        # it performed a clean exit and is no longer available
        comm.send(None, dest=0, tag=tags.EXIT)
        sys.exit(0)

    def worker_computation(self,data_from_master):

      arguments,job_index = data_from_master
      args,data_grid,data_grid_near,sample_coords,ndim,scale,min_thresh,max_thresh,id,nx,ny,nz,bounds = arguments  
      basecells = list()
      boundary_flags = list()
      jobs_per_data_transfer =int(ceil(nx*ny*nz/(comm.Get_size()-1)))
      #print (jobs_per_data_transfer)
      for index in range(jobs_per_data_transfer):
        index_in_worker = job_index+((comm.Get_size()-1)*index)
        if index_in_worker>= nx*ny*nz :
          continue 

        i, j, k = get_3d_grid_coords(index_in_worker,nx,ny,nz)
        this = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k))
        east = get_interp_3darray_elem(data_grid,data_grid_near,(i+1,j,k))
        top = get_interp_3darray_elem(data_grid,data_grid_near,(i,j+1,k))
        front = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k+1))
        
        assert(this is not None and east is not None and top is not None and front is not None)

        # if one of the values is < min_thresh, set it to min_thresh        
        # if one of the values is > max_thresh, set it to max_thresh
        x1 = min(max(this[0],min_thresh),max_thresh)
        x2 = min(max(east[0],min_thresh),max_thresh)
        y1 = min(max(this[1],min_thresh),max_thresh)
        y2 = min(max(top[1],min_thresh),max_thresh)
        z1 = min(max(this[2],min_thresh),max_thresh)
        z2 = min(max(front[2],min_thresh),max_thresh)
        
        bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta)
        bc_input.eta = 0.7
        bc_input.stiffness_as_diameter = True
        cell_obj = basecell.Basecell(bc_input)
        if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
      #       if cell_obj.volume <= args.bc_volume_thresh:
           continue
        
        # translate cell to correct position
        coord  = numpy.asarray(sample_coords[i,j,k])
        
        # flags for meshing circles on the boundary
        flag = [None] * 6
        for c,bound in enumerate(bounds):
          if numpy.isclose(coord[int(c/2)],bound):
            flag[c] = True
        
        #with p.lock:
        basecells.append((cell_obj.points,cell_obj.cells,coord))
        boundary_flags.append(flag)
        print("appended ",i,j,k,coord,x1,x2,y1,y2,z1,z2)

      del coord
      del bounds
      del bc_input
      del cell_obj
      gc.collect()
      return basecells, boundary_flags

# master class
class Master():
    def __init__ (self, args_for_worker):
      self.args_for_worker=args_for_worker
      args,data_grid,data_grid_near,sample_coords,ndim,scale,min_thresh,max_thresh,id,nx,ny,nz,bounds = self.args_for_worker

      # the master's ID
      self.id = comm.Get_rank()
      # ID should always be 0 due to the initialization in the main function
      assert(self.id==0)
      
      self.number_of_jobs = nx*ny*nz
      # a list where jobs will be put
      # self.tasks = [None] * self.number_of_jobs
      # a list where results will be put
      self.basecells =list()
      self.boundary_flags= list()
     
        
    def run(self):
        # put some jobs in the task list
        # for ii in range(self.number_of_jobs):
            # self.tasks[ii] = ii
        print('The no of jobs is {}.'.format(self.number_of_jobs))

        # a status object containing source, tag and size of a message
        status = MPI.Status()
        # the number of workers is the size of the communicator minus the master
        number_of_workers = comm.Get_size() - 1
        print('There are {} workers in this communicator.'.format(number_of_workers))
        # at start we have no closed workers
        closed_workers = 0

        # index for the current job to process
        job_index = 0
        # as long as there are workers available do something
        while closed_workers < number_of_workers:
            # wait for a message from any source with any tag
            data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            # get the source's ID / rank
            source = status.Get_source()
            # get the message tag
            tag = status.Get_tag()
            
            # if the received tag is 'READY' the worker waits for data to process
            if tag == tags.READY:
                if job_index < (number_of_workers):
                    # get a task from the list
                    # job_data = self.tasks[job_index]
                    # assemble the message to send to the worker
                    args_for_worker = self.args_for_worker,job_index
                    # send a message to the worker and tell it to process the
                    # contained data ('args')
                    comm.send(args_for_worker, dest=source, tag=tags.START)
                    # increase the index such that the following job is picked next
                    job_index += 1

                # if we have no job left we send a poison pill to the worker
                else:
                    comm.send(None, dest=source, tag=tags.EXIT)

            # if the received  tag is 'DONE' the worker finished computation
            # and sent a result
            elif tag == tags.DONE:
                # we store the result
                tmp_basecells, tmp_boundary_flags = data
                self.basecells.extend(tmp_basecells)
                self.boundary_flags.extend(tmp_boundary_flags)

            # if the received  tag is 'EXIT' the worker performed a clean exit
            elif tag == tags.EXIT:
                # we print a short message...
                print('Worker {} exited.'.format(source))
                # ...and increase the number of unavailable workers
                closed_workers += 1
              
        return self.basecells , self.boundary_flags