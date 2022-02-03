from matviz_rot import *
import numpy
import scipy.interpolate as ip
import time
import vtk
from os.path import splitext
import vtk.util.numpy_support
#from vtk.util.colors import *

# # creates 3D data to vtkPolyData
def create_vtk_poly_data(angle, data):

  assert(len(angle) == len(data))
  samples = int(sqrt(len(angle)))

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
def create_symmetry_planes(minima, scale, add_planes):
  # code source: http://code.google.com/p/pythonxy/source/browse/src/python/vtk/DOC/Examples/Rendering/Cylinder.py

  minima = []
  minima.append(((0, 0), 1))
  minima.append(((0, numpy.pi / 2), 1))
  minima.append(((numpy.pi / 2, numpy.pi / 2), 1))

  actors = []

  for mini in minima:

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
    nc = vtk.vtkNamedColors().GetColor3d('black')
    actor_c.GetProperty().SetColor(nc)
    angle = mini[0]
    phi = angle[0]
    theta = angle[1]


    print("angle: " + str(angle) + " -> " + str(to_vector(angle)) + " = " + str(mini[1]))

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

      nc = vtk.vtkNamedColors().GetColor3d('banana')
      actor_d.GetProperty().SetColor(nc)
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


# small helper to genetate a 3d wire frame for a box as polydata.
# see vtk.vtkAppendPolyData()
def generate_outline_box(size = [1,1,1], offset = [0,0,0]):

  # see https://www.vtk.org/Wiki/VTK/Examples/Python/GeometricObjects/Display/ColoredLines

  pts = vtk.vtkPoints()
  for x in [offset[0], offset[0] + size[0]]:
    for y in [offset[1], offset[1] + size[1]]:
     for z in [offset[2], offset[2] + size[2]]:
        pts.InsertNextPoint(x, y, z)
        ##print(x,y,z)
  #0: 0 0 0
  #1: 0 0 1
  #2: 0 1 0
  #3: 0 1 1
  #4: 1 0 0
  #5: 1 0 1
  #6: 1 1 0
  #7: 1  1 1

  # point codings for the lines
  edges = [[0,4], [4,6], [6,2], [2,0], [0,1], [4,5], [6, 7], [2,3], [1,5], [5, 7], [7,3], [3,1]]
  lines = vtk.vtkCellArray()
  for edge in edges:
    line = vtk.vtkLine()
    line.GetPointIds().SetId(0,edge[0])
    line.GetPointIds().SetId(1,edge[1])
    lines.InsertNextCell(line)

  poly = vtk.vtkPolyData()
  poly.SetPoints(pts)
  poly.SetLines(lines)

  colors = vtk.vtkUnsignedCharArray()
  colors.SetNumberOfComponents(3)
  colors.SetName("color")
  nc = vtk.vtkNamedColors().GetColor3d('black')
  for c in edges:
    colors.InsertNextTuple3(nc[0], nc[1], nc[2])

  poly.GetCellData().SetScalars(colors)

  return poly

# show the data on the screen
# @planes list of vtk actors containing symmetry planes
def show_vtk(polydata, res, save = None, planes=[], show_edges=False, show_axes=False, camera_settings=None):
  # Create a mapper and actor
  mapper = vtk.vtkPolyDataMapper()
  if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(polydata)
  else:
    if show_axes:
        # Scale (normalize) the data. Scaling the axes actor instead would destroy
        # the rendered coordinate axes for very small scaling factors.
        bounds = polydata.GetBounds()
        scale = 1.0/numpy.min([bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]])
        transform = vtk.vtkTransform()
        transform.Scale(scale, scale, scale)
        filter = vtk.vtkTransformFilter()
        filter.SetInputData(polydata)
        filter.SetTransform(transform)

        mapper.SetInputConnection(filter.GetOutputPort())
    else:
      mapper.SetInputData(polydata)

  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(0.5, 0.5, 0.5)  # (R,G,B)

  if show_edges: # show surface with edges
    actor.GetProperty().EdgeVisibilityOn()

  if show_axes:
      # Create axes actor
      axes = vtk.vtkAxesActor()
      bounds = polydata.GetBounds()
      # Scale axes
      scale *= 1.5
      length = numpy.array([bounds[1]*scale, bounds[3]*scale, bounds[5]*scale])
      tipLength = axes.GetNormalizedTipLength() / length * length[numpy.argmin(length)]
      axes.SetTotalLength(length)
      axes.SetNormalizedTipLength(tipLength)
      axes.SetNormalizedShaftLength(1-tipLength)
      # Set axis label color
      axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,0,0)
      axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,1,0)
      axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,1)

  # Setup a renderer, render window and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  # renderWindow.SetWindowName("Test")

  renderWindow.AddRenderer(renderer);
  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)

  renderWindow.SetSize(res, res)

  # Add the actor to the scene
  for plane in planes:
    renderer.AddActor(plane)

  renderer.AddActor(actor)
  if show_axes:
    renderer.AddActor(axes)

  renderer.SetBackground(1, 1, 1)  # Background color white
#  renderWindow.SetAlphaBitPlanes(1)

  # Render and interact
  renderWindow.Render()

  if camera_settings:
    cam = renderer.GetActiveCamera();
    cam.SetPosition(camera_settings[0], camera_settings[1], camera_settings[2])
    cam.SetFocalPoint(camera_settings[3], camera_settings[4], camera_settings[5])
    cam.SetRoll(camera_settings[6])

  if save:
    save_screenshot(renderWindow, save)
  else:
    #*** SetWindowName after renderWindow.Render() is called***
    # renderWindow.SetWindowName("Test")
    renderWindowInteractor.Start()

  # Display camera settings
  cam = renderer.GetActiveCamera()
  cam_settings = '{:f} {:f} {:f}'.format(*cam.GetPosition())
  cam_settings += ' {:f} {:f} {:f}'.format(*cam.GetFocalPoint())
  cam_settings += ' {:f}'.format(cam.GetRoll())
  print('camera when closing window: ' + cam_settings)


def save_screenshot(renderWindow, save):
  renderWindow.SetOffScreenRendering(True)
  renderWindow.Render()

  windowToImageFilter = vtk.vtkWindowToImageFilter()
  windowToImageFilter.SetInput(renderWindow)
  windowToImageFilter.Update()

  fileType = splitext(save)[1]
  if fileType == '.png':
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(save)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
  elif fileType == '.pdf':
    writer = vtk.vtkGL2PSExporter()
    writer.SetSilent(True)
    writer.SetFileFormatToPDF()
    writer.SetTitle(save)
    writer.SetRenderWindow(renderWindow)
    writer.SetFilePrefix(splitext(save)[0])
    writer.Write()
  elif fileType == '.eps':
    writer = vtk.vtkGL2PSExporter()
    writer.SetSilent(True)
    writer.SetFileFormatToEPS()
    writer.SetTitle(save)
    writer.SetRenderWindow(renderWindow)
    writer.SetFilePrefix(splitext(save)[0])
    writer.Write()
  else:
    print('File type {:s} not supported. Choose one out of (png).'.format(fileType))


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
  R = numpy.zeros((3, 3))
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

  # use normal list with tuples of coords/cell vertices as we cannot parallelize vtk objects
  points_list = []
  cells_list = []

  #optional parameters: @param  not_drawn (faces which are not drawn)
  #                     @param angle list of angle_x, angle_y, angle_z or None
  for cc in coords:
    p, c = create_centered_bar(cells, points, cc, dim, angle, not_drawn)
    points_list.extend(p)
    cells_list.extend(c)

  return points_list, cells_list

def create_centered_bar(cells, points, center, dim, angle=None, not_drawn = None):
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

  # use normal list with tuples of coords/cell vertices as we cannot parallelize vtk objects
  points_list = []
  cells_list = []

  base = points.GetNumberOfPoints()
  # calculate corner points of quad and add them to global points list
  point_vec = create_point_vector_centered_bar(center, dim, angle)
  for p in point_vec:
    points.InsertNextPoint(p)  # 0 ... 7
    points_list.append(p)

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
    cells_list.append((base + 0,base + 1,base + 2,base + 3))

  if not_drawn is None or not_drawn[0] != 'right' or not_drawn[1] != 'right':
    # right face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 1)
    quad.GetPointIds().SetId(1, base + 4)
    quad.GetPointIds().SetId(2, base + 6)
    quad.GetPointIds().SetId(3, base + 2)
    cells.InsertNextCell(quad)
    cells_list.append((base + 1,base + 4,base + 6,base + 2))

  if not_drawn is None or not_drawn[0] != 'back' or not_drawn[1] != 'back':
    # back face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 4)
    quad.GetPointIds().SetId(1, base + 6)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 5)
    cells.InsertNextCell(quad)
    cells_list.append((base + 4,base + 6,base + 7,base + 5))

  if not_drawn is None or not_drawn[0] != 'top' or not_drawn[1] != 'top':
    # top face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 2)
    quad.GetPointIds().SetId(1, base + 6)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 3)
    cells.InsertNextCell(quad)
    cells_list.append((base + 2,base + 6,base + 7,base + 3))

  if not_drawn is None or not_drawn[0] != 'left' or not_drawn[1] != 'left':
    # left face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 0)
    quad.GetPointIds().SetId(1, base + 5)
    quad.GetPointIds().SetId(2, base + 7)
    quad.GetPointIds().SetId(3, base + 3)
    cells.InsertNextCell(quad)
    cells_list.append((base + 0,base + 5,base + 7,base + 3))

  if not_drawn is None or not_drawn[0] != 'bottom' or not_drawn[1] != 'bottom':
    # bottom face
    quad = vtk.vtkQuad()
    quad.GetPointIds().SetId(0, base + 0)
    quad.GetPointIds().SetId(1, base + 1)
    quad.GetPointIds().SetId(2, base + 4)
    quad.GetPointIds().SetId(3, base + 5)
    cells.InsertNextCell(quad)
    cells_list.append((base + 0,base + 1,base + 4,base + 5))

  return points_list, cells_list

def create_block(coords, design, scale, thres = 0.0, elems_in_regions = None):

  s1 = design['s1']

  centers, min, max, elem_dim = coords

  dx = elem_dim[0]
  dy = elem_dim[1]
  dz = elem_dim[2]

  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()

  if scale <= 0:
    scale = 1.0

  for s, coord in zip(s1,centers):
    if s > thres:
      create_centered_bar(cells, points, coord, (scale * dx, scale * dy, scale * dz))

  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)

  colors = vtk.vtkUnsignedCharArray()
  colors.SetNumberOfComponents(3)
  colors.SetName("color")
  if not elems_in_regions:
    # only one color
    nc = vtk.vtkNamedColors().GetColor3d('black')
    for _ in s1:
      colors.InsertNextTuple3(nc[0], nc[1], nc[2])
  else:
    # one color for each region
    # generate list of RGB colors
    n_regions = size(elems_in_regions)
    n_color = int(ceil(pow(n_regions,1./3.)))
    color1D = numpy.linspace(0, 255, n_color)
    nc = numpy.array(numpy.meshgrid(color1D, color1D, color1D)).T.reshape(-1,3)
    for index, s in enumerate(s1):
      for ri, elems in enumerate(elems_in_regions):
        if index+1 in elems: # elems are 1 based
          regindex = ri
          break
      if s > thres:
        colors.InsertNextTuple3(nc[regindex][0], nc[regindex][1], nc[regindex][2])

    polydata.GetCellData().SetScalars(colors)

  return polydata

# # without rotation and shearing
def create_3d_frame(coords, design, dir, scale):

  s1 = design['s1']
  s2 = design['s2']
  s3 = design['s3']
  angles = design['angle']

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

def sign(p1, p2, p3):
  return (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1]);

def point_in_triangle(pt, v1, v2, v3):
  b1 = sign(pt, v1, v2) < 0.;
  b2 = sign(pt, v2, v3) < 0.;
  b3 = sign(pt, v3, v1) < 0.;
  return ((b1 == b2) and (b2 == b3));

def test_point_outside_circle(midpoint, radius, point):
  if (point[0] - midpoint[0]) ** 2 + (point[1] - midpoint[1]) ** 2 < radius ** 2:
    return False
  else:
    return True


def valid_bar_position_apod6(coords, center, dim, angle=None):
  point_vec = create_point_vector_centered_bar(center, dim, angle)
  valid = True
  for p in point_vec:
    if not valid_position_apod6(p, coords):
      valid = False
  return valid

# # for the robot arm we have check for the two nondesign holes as they are within the
# # convex hull of the design :(
def valid_position_robot(pos, coords, opt=0):
  #mi, ma = coords[1:3]
  #delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
 # if int(delta[0]) == 508 and int(delta[2]) == 126:
  #if (pos[0] + 147.4) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center -147, 0, 0
  if (pos[0] + 147.4) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center -147, 0, 0
    return False
  if (pos[0] - 250.0) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center 250, 0, 0
    return False
  return True

def valid_position_lufo(pos, coords):
  dx = 2.
  dy = 2.
  dz = 2.
  sup5 = [24.,2.,-12.]
  sup4 = [24.,2.,-47.]
  sup3 = [24.,2.,-82.]
  sup2 = [24.,2.,-117.]
  sup1 = [24.,2.,-152.]
  force = [0.0, 37.49281,-30.1963]
  r = 8. - dx/2.

  # cut out support ring
  if (pos[1] < 6.00001):
    if (pos[0] - sup1[0]) ** 2 + (pos[2]-sup1[2]) ** 2 < r ** 2:
      return False
    if (pos[0] - sup2[0]) ** 2 + (pos[2]-sup2[2])  ** 2 < r ** 2:
      return False
    if (pos[0] - sup3[0]) ** 2 + (pos[2]-sup3[2])  ** 2 < r ** 2:
      return False
    if (pos[0] - sup4[0]) ** 2 + (pos[2]-sup4[2])  ** 2 < r ** 2:
      return False
    if (pos[0] - sup5[0]) ** 2 + (pos[2]-sup5[2])  ** 2 < r ** 2:
      return False
  # cut out force ring and surrounding area
  if (pos[0] < 33.91 and pos[0] > -33.89 and pos[1] > 37.783 and pos[1] < 66.01 and pos[2] > -28.99 and pos[2] < 0.01):
    return False
  m2 = [0,37.49281,-30.1963]
  r2 = 27.01 - dx/2.
  # cut out big cylinder near force
  if (pos[1] - m2[1]) ** 2 + (pos[2] - m2[2]) ** 2 < r2 ** 2:
    return False
  m3 = [ 24.000299, 8.000000, -12.000068]
  r3 = 12.51 - 2.5*dx
  # cut out box above support5
  if (pos[0] < 35.01- dx/2. and pos[0] > 22.99 + dx/2. and pos[1] > 5.99 + dx/2. and pos[1] < 10.01 - dx/2. and pos[2] > -24.51 + dx/2.and pos[2] < -0.49 - dx/2.) or ((pos[0] - m3[0]) ** 2 + (pos[2]-m3[2]) ** 2 < r3 ** 2):
    return False
  m4 = [24.0, 8.0, -152.0]
  r4 = 12.51 - 2.5*dx
  # cut out box above support1
  if (pos[0] < 35.01 - dx/2. and pos[0] > 22.99 + dx/2. and pos[1] > 5.99 + dx/2. and pos[1] < 10.01 - dx/2. and pos[2] > -164.51 + dx/2. and pos[2] < -139.49 - dx/2.) or ((pos[0] - m4[0]) ** 2 + (pos[2]-m4[2]) ** 2 < r4 ** 2):
    return False
  m5 = [24.,  8., -117.]
  r5 = 10.01 - 1.5* dx
  # cut out box above support2
  if (pos[0] < 34.01 and pos[0] > 23.99 and pos[1] > 5.99 and pos[1] < 10.01 and pos[2] > -127.01 and pos[2] < -106.99) or ((pos[0] - m5[0]) ** 2 + (pos[2]-m5[2]) ** 2 < r5 ** 2):
    return False
  m6 = [ 24.0, 8.0, -82.]
  r6 = 10.01 - 1.5* dx
  # cut out box above support3
  if (pos[0] < 34.01 and pos[0] > 23.99 and pos[1] > 5.99 and pos[1] < 10.01 and pos[2] > -92.01 and pos[2] < -71.99) or ((pos[0] - m6[0]) ** 2 + (pos[2]-m6[2]) ** 2 < r6 ** 2):
    return False
  m7 = [ 24.0, 8.0, -47.]
  r7 = 10.01 - 1.5* dx
  # cut out box above support3
  if (pos[0] < 34.01 and pos[0] > 23.99 and pos[1] > 5.99 and pos[1] < 10.01 and pos[2] > -57.01 and pos[2] < -36.99) or ((pos[0] - m7[0]) ** 2 + (pos[2]-m7[2]) ** 2 < r7 ** 2):
    return False

  return True

def valid_ring_position_apod6(pos, coords, opt = 0. ):
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
def valid_position_apod6(pos, coords, opt = 0. ):
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
def create_3d_frame_ip(coords, design, ip, grad, scale, valid_position, thres=0.0, csize = None):
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

  s1 = design['s1']
  s2 = design['s2']
  s3 = design['s3']
  angles = design['angle']

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
          if True:# if s1 >= thres:#valid_bar_position_apod6(points,coord, (scale * scale_[0] * dx, scale * s1 * dx, scale * s1 * dx), angle):
            # draw thickest bars first
            coords = []
            for _ in range(4):
              coords.append(coord)
            # add offset for s1 for all s1-bars
            coords[0] = coords[0] + [0.,-s1/4.,s1/4.]
            coords[1] = coords[1] + [0.,s1/4.,s1/4.]
            coords[2] = coords[2] + [0.,-s1/4.,-s1/4.]
            coords[3] = coords[3] + [0.,s1/4.,-s1/4.]
            create_centered_bars(cells, points, coords, (scale * scale_[0] * dx, scale * 0.5 * s1 * dx, scale * 0.5 * s1 * dx), angle,['right','left'])
            real_volume += scale * scale_[0] * dx * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,scale* dx * s1 * 0.5 + scale * 0.25 * (scale_[1] * dy - dx * s1),0.]
          dy_offset = scale * 0.5 * (scale_[1]*dy-s1*dx)
          #add two parts of s2-bar, two parts are necessary that it doesn't intersect the s1-bar
          if True:# if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            coords = []
            for i  in range(4):
              coords.append(coord)
            # add offset for s2 for all s2-bars
            coords[0] = coords[0] + [-s1/4.,0., s1/4.]
            coords[1] = coords[1] + [s1/4., 0., s1/4.]
            coords[2] = coords[2] + [-s1/4.,0., -s1/4.]
            coords[3] = coords[3] + [s1/4., 0.,  -s1/4.]
            create_centered_bars(cells, points, coords, (scale * 0.5 * s2 * dy, dy_offset, scale * 0.5 * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy* dy_offset* scale * s2 * dy

          coord_offset = [0.,0.,scale * dx * s1 * 0.5 + scale * 0.25 * (scale_[2] * dz - dx * s1)]
          dz_offset = scale * 0.5 * (scale_[2]*dz-s1*dx)
          #add two parts of s3-bar, two parts are necessary that it doesn't intersect the s1-bar
          if True: # if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            coords = []
            for i  in range(4):
              coords.append(coord)
            # add offset for s3 for all s3-bars
            coords[0] = coords[0] + [-s1/4.,s1/4.,0.]
            coords[1] = coords[1] + [s1/4., s1/4.,0.]
            coords[2] = coords[2] + [-s1/4., -s1/4., 0.]
            coords[3] = coords[3] + [s1/4., -s1/4.,0.]
            create_centered_bars(cells, points, coords, (scale * 0.5 * s3 * dz, scale * 0.5 * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset

        # draw each bar of 3D cross for s2 > s1,s3
        elif s2 >= s1 and s2 >= s3:
          if True: #if s2 >= thres:#valid_bar_position_apod6(points,coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy * scale * scale_[1]* dy * scale * s2 * dy
          coord_offset = [scale* dy * s2 * 0.5 + scale * 0.25 * (scale_[0] * dx - dy * s2),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx-s2*dy)
          if True: #if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          if True: # if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,0.,scale * dy * s2 * 0.5 + scale * 0.25 * (scale_[2] * dz - dy * s2)]
          dz_offset = scale * 0.5 * (scale_[2] * dz-s2*dy)
          if True: #if s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
          if True: # if s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
        # draw each bar of 3D cross for s3 > s1,s2
        elif s3 >= s1 and s3 >= s2:
          if True: #if s3 >= thres:#valid_bar_position_apod6(points,coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle):
            create_centered_bar(cells, points, coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * scale * scale_[2] * dz
          coord_offset = [scale* dz * s3 * 0.5 + scale * 0.25 * (scale_[0] * dx - dz * s3),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx-s3*dz)
          if True: #if s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          if True:# if s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dx * scale * s1 * dx
          coord_offset = [0.,scale * dz * s3 * 0.5 + scale * 0.25 * (scale_[1] * dy - dz * s3),0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy-s3*dz)
          if True: #if s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle,['top','bottom'])
            real_volume += scale * s2 * dy * dy_offset * scale * s2 * dy
          if True: #if s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
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
def create_3d_cross_ip(coords, design, ip, grad, scale, valid_position, thres=0.0, csize = None):
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip: number of uniform cells in (x,y,z)-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # valid_position: returns false if point inside the convex hull of the part, should be excluded, otherwise true.
  #                 Needs to be implemented for every mechanical part, currently available for robot and apod6.
  #                 If part is not implemented valid_position is None and no cells inside the convex hull are removed from the structure
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]

  s1 = design['s1']
  s2 = design['s2']
  s3 = design['s3']
  angles = design['angle']

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
          if True:#s1 >= thres:#valid_bar_position_apod6(points,coord, (scale * scale_[0] * dx, scale * s1 * dx, scale * s1 * dx), angle):
            # draw thickest bar first
            create_centered_bar(cells, points, coord, (scale * scale_[0] * dx, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += scale * scale_[0] * dx * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., scale * s1 * dy * 0.5 + scale * 0.25 * (scale_[1] * dy - s1 * dy), 0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy - s1 * dy)
          #add two parts of s2-bar, two parts are necessary that it doesn't intersect the s1-bar
          if True:#s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * dy_offset * scale * s2 * dz
          if True:#s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * dy_offset * scale * s2 * dz

          coord_offset = [0., 0., scale * s1 * dz * 0.5 + scale * 0.25 * (scale_[2] * dz - s1 * dz)]
          dz_offset = scale * 0.5 * (scale_[2] * dz - s1 * dz)
          #add two parts of s3-bar, two parts are necessary that it doesn't intersect the s1-bar
          if True:#s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
          if True:#s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
        # draw each bar of 3D cross for s2 > s1,s3
        elif s2 >= s1 and s2 >= s3:
          if True:#s2 >= thres:#valid_bar_position_apod6(points,coord, (scale * s2 * dy, scale * scale_[1]* dy, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord, (scale * s2 * dx, scale * scale_[1] * dy, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dx * scale * scale_[1] * dy * scale * s2 * dz
          coord_offset = [scale * s2 * dx * 0.5 + scale * 0.25 * (scale_[0] * dx - s2 * dx), 0., 0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx - s2 * dx)
          if True:#s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          if True:#s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset, scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., 0., scale * s2 * dz * 0.5 + scale * 0.25 * (scale_[2] * dz - s2 * dz)]
          dz_offset = scale * 0.5 * (scale_[2] * dz - s2 * dz)
          if True:#s3 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dz * scale * s3 * dz * dz_offset
          if True:#s3 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s3 * dz, scale * s3 * dz,dz_offset), angle):
            create_centered_bar(cells, points, coord - coord_offset, (scale * s3 * dx, scale * s3 * dy, dz_offset), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * dz_offset
        # draw each bar of 3D cross for s3 > s1,s2
        elif s3 >= s1 and s3 >= s2:
          if True:#s3 >= thres:#valid_bar_position_apod6(points,coord, (scale * s3 * dz, scale * s3 * dz, scale * scale_[2] * dz), angle):
            create_centered_bar(cells, points, coord, (scale * s3 * dx, scale * s3 * dy, scale * scale_[2] * dz), angle,['front','back'])
            real_volume += scale * s3 * dx * scale * s3 * dy * scale * scale_[2] * dz
          coord_offset = [scale * s3 * dx * 0.5 + scale * 0.25 * (scale_[0] * dx - s3 * dx),0.,0.]
          dx_offset = scale * 0.5 * (scale_[0] * dx - s3 * dx)
          if True:#s1 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord + coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          if True:#s1 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (dx_offset,scale * s1 * dx, scale * s1 * dx), angle):
            create_centered_bar(cells, points, coord - coord_offset, (dx_offset, scale * s1 * dy, scale * s1 * dz), angle,['right','left'])
            real_volume += dx_offset * scale * s1 * dy * scale * s1 * dz
          coord_offset = [0., scale * dy * s3 * 0.5 + scale * 0.25 * (scale_[1] * dy - s3 * dy), 0.]
          dy_offset = scale * 0.5 * (scale_[1] * dy - s3 * dy)
          if True:#s2 >= thres:#valid_bar_position_apod6(points,coord + coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
            create_centered_bar(cells, points, coord + coord_offset, (scale * s2 * dx, dy_offset, scale * s2 * dz), angle,['top','bottom'])
            real_volume += scale * s2 * dy * dy_offset * scale * s2 * dy
          if True:#s2 >= thres:#valid_bar_position_apod6(points,coord - coord_offset, (scale * s2 * dy, dy_offset, scale * s2 * dy), angle):
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

# this is copy & paste from matviz_2d but extended to 3D
# @param nx_ip number of interpolations within x
def get_interpolation_natural_neighbor(coords, s1, s2, s3, dx, dy, dz):
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
        out[idx] = ((mi[0] + float(x) / nx * delta[0], mi[1] +  float(y) / ny * delta[1], mi[2] + float(z) / nz * delta[2]))
        idx += 1
  if s2 is None and s3 is None:
      v = numpy.zeros((len(s1), 1))
      v[:, 0] = s1[:, 0]
  else:
    v = numpy.zeros((len(s1), 3 ))
    v[:, 0] = s1[:, 0]
    v[:, 1] = s2[:, 0]
    v[:, 2] = s3[:, 0]

  import naturalneighbor
  l1 = numpy.array([s[0] for s in s1])
  l2 = numpy.array([s[0] for s in s2])
  l3 = numpy.array([s[0] for s in s3])

  l = numpy.array([numpy.array((s1[i],s2[i],s3[i])) for i in range(len(s1))])

  # y-z face: x has offset 0
  tmp_1 = naturalneighbor.griddata(numpy.array(centers), l1, [ [mi[0],mi[0]+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_2 = naturalneighbor.griddata(numpy.array(centers), l2, [ [mi[0],mi[0]+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_3 = naturalneighbor.griddata(numpy.array(centers), l3, [ [mi[0],mi[0]+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  yz = numpy.stack((tmp_1,tmp_2,tmp_3), axis=-1)

  # x-z face: y has offset 0
  tmp_1 = naturalneighbor.griddata(numpy.array(centers), l1, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1],mi[1]+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_2 = naturalneighbor.griddata(numpy.array(centers), l2, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1],mi[1]+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_3 = naturalneighbor.griddata(numpy.array(centers), l3, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1],mi[1]+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2]+ 0.5*dz,mi[2]+ 0.5*dz+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  xz = numpy.stack((tmp_1,tmp_2,tmp_3), axis=-1)

  # x-y face: z has offset 0
  tmp_1 = naturalneighbor.griddata(numpy.array(centers), l1, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2],mi[2]+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_2 = naturalneighbor.griddata(numpy.array(centers), l2, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2],mi[2]+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  tmp_3 = naturalneighbor.griddata(numpy.array(centers), l3, [ [mi[0]+ 0.5*dx,mi[0]+ 0.5*dx+delta[0]+1.0/nx,1.0/nx*delta[0]], [mi[1] + 0.5*dy,mi[1]+ 0.5*dy+delta[1]+1.0/ny,1.0/ny*delta[1]], [mi[2],mi[2]+delta[2]+1.0/nz,1.0/nz*delta[2]]])
  xy = numpy.stack((tmp_1,tmp_2,tmp_3), axis=-1)

  return [yz, xz, xy], out

# this is copy & paste from matviz_2d but extended to 3D
# @param nx_ip number of interpolations within x
# assume we get element data at barycenters (coords)
def get_3d_interpolation_at_faces(coords, bounds, grad, s1, s2, s3, nx, ny, nz, dx, dy, dz):
  # we make our own regular element grid
  centers, _, _ = coords[0:3]  # skip elem

  mi = bounds[0:3]
  ma = bounds[3:6]


  # the coordinate system in Paraview is a right-hand sided coodrdinate system with z pointing to the viewer
  #
  #  y ^
  #    |
  # z (.)--> x
  #
  # This are the node numbers if we have only one element. The .mesh file will be transformed to 1-based
  # x is the fastet variable, z is the slowest variable
  #
  #        -----------
  #      /|          /|      # left face   (x=0):  id = 0
  #     / |    3    / |      # right face  (x=1):  id = 1
  #     -----------   |      # bottom face (y=0):  id = 2
  #    |  |        |  |      # top face    (y=1):  id = 3
  #    |0 |    4   | 1|      # back face   (z=0):  id = 4
  #    |   --------|--       # front face  (z=1):  id = 5
  #    | /         | /
  #    |/     2    |/
  #     -----------  front face(z=1) has id 5

  # loop over elements, for each element we store 3 faces (left,bottom,back) and
  # each face is described by its barycenter
  out = []
  for x in range(nx+1):
    for y in range(ny+1):
      for z in range(nz+1):
        # center of left face
        out.append((mi[0] + x*dx, mi[1] + y * dy + 0.5*dy, mi[2] + z*dz + 0.5*dz))
        # center of bottom face
        out.append((mi[0]+x*dx + 0.5*dx, mi[1] + y*dy, mi[2] + z * dz + 0.5*dz))
        # center of back face
        out.append((mi[0] + x*dx + 0.5*dx, mi[1] + y * dy + 0.5*dy, mi[2] + z*dz))

  if s2 is None and s3 is None:
      v = numpy.zeros((len(s1), 1))
      v[:, 0] = s1[:, 0]
  else:
    v = numpy.zeros((len(s1), 3))
    v[:, 0] = s1#[:]#, 0]
    v[:, 1] = s2#[:]#, 0]
    v[:, 2] = s3#[:]#, 0]

  ip_data = ip.griddata(centers, v, out, grad, -1.0)

  # any interpolation, ie. linear interpolation can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation which can also interpolate values outside the convex hull
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad != 'nearest' else None

  # for each element, store center coordinates of its three faces
  ip_data = numpy.reshape(ip_data, ((nx+1),(ny+1),(nz+1),3,3))
  ip_near = numpy.reshape(ip_near, ((nx+1),(ny+1),(nz+1),3,3))
  out = numpy.reshape(out,((nx+1),(ny+1),(nz+1),3,3))


  return ip_data, ip_near, out

# interpolate given data at cell barycenters to cell vertices
# elems: elements of desgin domain described by its vertices
# reg_info: info on design region such as node coords, elements and connectivity
def interp_cell_to_point_data(center_coords, reg_info, bounds, grad, s1, s2, s3, nx, ny, nz, dx, dy, dz):
  # we make our own regular element grid
  # coordinates of barycenters at which the values for s1, s2, s3 are given
  centers = center_coords
  mi = bounds[0:3]
  ma = bounds[3:6]

  # coords of vertices defining elements
  elems = reg_info['elements']
  nodes = list(reg_info['nodes'])
  # stores vertex ids of each element
  connectivity = reg_info['connectivity']
  # region_map: maps local node id (list idx) in given design region (e.g. "mech") to global node id in all regions
  local_to_global_id = reg_info['region_map']
  if any(isinstance(el, numpy.ndarray) for el in local_to_global_id):
    local_to_global_id = [el[0] for el in local_to_global_id]
  global_to_local_id = dict()
  for i,id in enumerate(local_to_global_id):
    global_to_local_id[id-1] = i
#   print("nodes:",nodes)

  assert(len(elems) == len(connectivity))
  assert(len(elems) == len(s1))
  assert(len(elems) == len(s2))
  assert(len(elems) == len(s3))

  # for each global node: calc weighted sum of design variables of adjacent elements
  s1_nodes = numpy.zeros(len(nodes))
  s2_nodes = numpy.zeros(len(nodes))
  s3_nodes = numpy.zeros(len(nodes))
  neighbors = numpy.zeros(len(nodes))

  for ie,e in enumerate(connectivity): # each element
    for v in e: # each element vertex
      s1_nodes[global_to_local_id[v-1]] += s1[ie]#/numpy.linalg.norm(numpy.array(centers[ie]) - numpy.array(nodes[v-1])) # cfs is 1-based
      s2_nodes[global_to_local_id[v-1]] += s2[ie]#/numpy.linalg.norm(numpy.array(centers[ie]) - numpy.array(nodes[v-1]))
      s3_nodes[global_to_local_id[v-1]] += s3[ie]#/numpy.linalg.norm(numpy.array(centers[ie]) - numpy.array(nodes[v-1]))
      neighbors[global_to_local_id[v-1]] += 1

  # normalize
  for i,_ in enumerate(neighbors):
    assert(neighbors[i] > 0)
    s1_nodes[i] /= neighbors[i]
    s2_nodes[i] /= neighbors[i]
    s3_nodes[i] /= neighbors[i]

  return nodes, s1_nodes, s2_nodes, s3_nodes


  # loop over elements, for each element we store 3 faces (left,bottom,back) and
  # each face is described by its barycenter
  out = []
  for x in range(nx+1):
    for y in range(ny+1):
      for z in range(nz+1):
        # center of left face
        out.append((mi[0] + x*dx, mi[1] + y * dy + 0.5*dy, mi[2] + z*dz + 0.5*dz))
        # center of bottom face
        out.append((mi[0]+x*dx + 0.5*dx, mi[1] + y*dy, mi[2] + z * dz + 0.5*dz))
        # center of back face
        out.append((mi[0] + x*dx + 0.5*dx, mi[1] + y * dy + 0.5*dy, mi[2] + z*dz))

  ip_data= naturalneighbor.griddata(centers, v, out)

# # for a rectangular grid calculate the left, bottom and back face center coords of one element
# # barycenters:
# def face_center_coords_cube(barycenters,bounds, nx,ny,nz,dx,dy,dz):

# # litte helper
# @param save filename or none
# @param list which might be empty
def show_write_vtk(poly, res, save, actors=[], show_axes=False, camera_settings=None):
  if save and save.split('.')[-1] == 'vtp':
    writer = vtk.vtkXMLPolyDataWriter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(poly)
    else:
      writer.SetInputData(poly)
    #writer.SetDataModeToAscii()
    writer.SetFileName(save)
    writer.Write()
    print("saved polydata to '" + save + "'")
  else:
    show_vtk(poly, res, save, actors, show_axes=show_axes, camera_settings=camera_settings)

def calc_cross_elem_vol_3D(s1, s2, s3):
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
def write_stl(polydata, save=None):
  stlWriter =  vtk.vtkSTLWriter()
  fName = save if save else 'surface.stl'
  if not fName.endswith(".stl"):
    fName += ".stl"
  stlWriter.SetFileName(fName)
  stlWriter.SetInputData(polydata)
  stlWriter.SetFileTypeToASCII()
  stlWriter.Write()

  print("saved polydata to file " + fName)

# take points and cells list/arrays and create polydata
# stores shortest edge length in 'short'
# optional pointIds, in case point ids don't match with list index of 'points'
def fill_vtk_polydata(points,cells,pointIds=None):
  vtk_points = vtk.vtkPoints()
  vtk_cells = vtk.vtkCellArray()
  polydata = vtk.vtkPolyData()

  if not pointIds:
    for p in points:
      vtk_points.InsertNextPoint(p)
  else:
    vtk_points.SetNumberOfPoints(len(points))
    for p,pid in zip(points,pointIds):
      vtk_points.InsertPoint(pid, p)

  for ce in cells:
    if len(ce) == 3: # triangle
      add_triangle(ce[0],ce[1],ce[2],vtk_cells)
    elif len(ce) == 4: # quad
      add_triangle(ce[0],ce[1],ce[2],vtk_cells)
      add_triangle(ce[2],ce[3],ce[0],vtk_cells)
#       quad = vtk.vtkQuad()
#       quad.GetPointIds().SetId(0, ce[0])
#       quad.GetPointIds().SetId(1, ce[1])
#       quad.GetPointIds().SetId(2, ce[2])
#       quad.GetPointIds().SetId(3, ce[3])
#       vtk_cells.InsertNextCell(quad)
    else:
      print("fill_vtk_polydata: Ohoh ERROR")
      sys.exit()

  polydata.SetPoints(vtk_points)
  polydata.SetPolys(vtk_cells)

  return polydata

def add_triangle(id1,id2,id3,cells):
  if id1 == id2 or id1 == id3:
    return

  i1 = id1
  i2 = id2
  i3 = id3

  tri = vtk.vtkTriangle()
  tri.GetPointIds().SetId(0, i1)
  tri.GetPointIds().SetId(1, i2)
  tri.GetPointIds().SetId(2, i3)
  cells.InsertNextCell(tri)

def vtk_polydata_to_numpy(polydata):
  from vtk.util.numpy_support import vtk_to_numpy
  # get point data
  points = vtk_to_numpy(polydata.GetPoints().GetData())
#   print("number of points:",len(points))

  polys = polydata.GetPolys()
  ne = polys.GetNumberOfCells()
#   print("number of polys:",polydata.GetNumberOfPolys())

  cells = numpy.zeros((ne,3),dtype=int)
  polys.InitTraversal()

  for i in range(ne):
    ids = vtk.vtkIdList()
    polys.GetNextCell(ids)
    ni = ids.GetNumberOfIds()
    for j in range(ni):
      cells[i,j] = ids.GetId(j)

  return points, cells
