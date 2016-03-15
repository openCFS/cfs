#!/usr/bin/env python
import platform
from PIL import Image
import sys, os, copy, numpy, math
from hdf5_tools import *
import scipy.interpolate as ip
from matviz_vtk import *
from scipy.spatial import Delaunay
from PyQt4.Qt import left
from matplotlib.sankey import RIGHT


# writes a dense two region mesh
# def write_dense_mesh(pixels, size, file, threshold):

# element types as in gid (simInputMESH.cc -> AnsysType2ElemType)
TRIANGLE3 = 4
QUAD4 = 6
TET4 = 8
HEXA8 = 10
WEDGE6 = 14
LINE = 100


def nodes_by_type(type):
  if type == QUAD4:
    return 4
  if type == HEXA8:
    return 8
  if type == WEDGE6:
    return 6
  if type == TRIANGLE3:
    return  3
  if type == LINE:
    return 2
  if type == TET4:
    return 4
  assert(False)
 
def mesh_type_from_hdf5(type_id):
  if type_id == 6:
    return QUAD4
  if type_id == 16:
    return WEDGE6
  if type_id == 11:
    return HEXA8
  if type_id == 4:
    return TRIANGLE3
  if type_id == 8:
    return TET4
  assert(False)

def elem_dim(type):
  if type == HEXA8 or type == WEDGE6 or type == TET4:
    return 3
  elif type == LINE:
    return 1
  else: 
    return 2 
  


# gid element
class Element: 
  def __init__(self):
    self.nodes = []  # list of zero based node indices. counter-clock wise
    self.region = None  # region name
    self.density = 1  # from lower_bound to 1, not necessarily used
    self.stiff1 = 0
    self.stiff2 = 0
    self.stiff3 = 0
    self.rotAngle = 0
    self.rotX = 0
    self.rotY = 0
    self.rotZ = 0
    self.type = -1
    
    
  def dump(self):
    print self.nodes
    print self.region
    print self.density     


# gid Mesh
class Mesh:
  def __init__(self, nx = -1, ny = -1, nz = -1):
   self.nodes = []  # list 2d tupels (float, float) or 3d tuples
   self.elements = []  # list of Element
   # list of boundary conditon nodes
   self.bc = []  # list of tupel (name, <list of zero based nodes>)
   # list of named nodes (save element in gd)
   self.ne = []  # list of tupel (name, <list of zero based elements>)
   self.nx = nx
   self.ny = ny
   self.nz = nz
  
  def element(self, i, j):
    assert(self.nx > 0 and self.ny > 0)
    return self.elements[i * self.nx + j]
  
def show_dense_mesh_image(mesh, shape, binary, size):
  check_img = Image.new("RGB", shape, "white")
  check_pix = check_img.load()

  nx, ny = shape
  
  for x in range(nx):
    for y in range(ny):
      # print input_pix[x,y]
      e = mesh.elements[x * ny + y]
      val = 1 - e.density  # black is 0 in the image but 1 as density
      # print str(val) + " - " + str(barrier)
      show = (200, 10, 10) if binary else (int(val * 255), int(val * 255), int(val * 255)) 
      check_pix[x, ny - y - 1] = show if e.region == 'mech' else (10, 10, 200) if e.region == 'void' else (200, 10, 10)
      
  check_img = check_img.resize((size, ny * size / nx))    
  check_img.show()     


def create_dense_mesh_img(input_img, mesh, threshold, scale, rhomin, shearAngle, type = 1, color_mode="random"):
  input_pix = input_img.load()
  nx, ny = input_img.size
  create_dense_mesh(input_pix, nx, ny, mesh, threshold, scale, rhomin, 1, shearAngle, type, color_mode)

def create_dense_mesh_density(numpy_array, mesh, threshold, scale, rhomin, multi_d=1):
  # handle different types of numpy_array, further description in create_dense_mesh
  # only one design variable
  if multi_d == 1:
    nx, ny = numpy_array.shape
  # multiple design variables
  else:
    #m,n are dummy variables
    nx, ny, m, n = numpy_array.shape
  create_dense_mesh(numpy_array, nx, ny, mesh, threshold, scale, rhomin, multi_design=multi_d, shearAngle=0.0)
  
def create_dense_mesh(input_array, nx, ny, mesh, threshold, scale, rhomin, multi_design=1, shearAngle=0, type=1, color_mode="random"):
  # convert angle to rad and check for feasibility
  angle = shearAngle / 180 * math.pi
  if (abs(angle) > math.pi / 2 - 1e-6):
    print 'angle has to be between -pi/2 + 1e-6 and pi/2 - 1e-6'
    return 0 
  dx = scale / nx / math.cos(angle)
  # from daniel ?! dy = scale/ny
  dy = dx
  
  # input_array can be one of three cases: grayscale imgae (array of ints), color image (array of tuples (r,g,b(,a)), numpy.ndarray
  is_data = isinstance(input_array, numpy.ndarray)
  is_gray = True if color_mode == "L" and not is_data else False
  is_color = True if color_mode == "RGB" and not is_data else False

  assert(is_data or is_gray or is_color)
  
  # create mesh.nodes
  for y in range(ny + 1):
    for x in range(nx + 1):
      if angle == 0.0:
        mesh.nodes.append((x * dx, y * dy))
      else:
        x_Coord = round(x * dx - y * dy * math.tan(angle), 8)
        if abs(x_Coord) < 1e-8:
          x_Coord = 0.0
        mesh.nodes.append((x_Coord, y * dy))
  # print mesh.nodes 
  mech_count = 0
  colorful_count = 0
  for x in range(nx):
    for y in range(ny):
      e = Element()
      e.type = QUAD4
      # assign preliminary data value
      if is_gray:
        # convert to black is one and white = 0
        e.density = 1 - (input_array[x, y] / 255.0)
      if is_color:
        val = sum(input_array[x, y][0:3]) / 3.0
        e.density = 1.0 - (val / 255.0)
      if is_data:
        if multi_design == 1:
          e.density = input_array[x, y]
        else:
          e.stiff1 = input_array[x, y, 0, 0]
          e.stiff2 = input_array[x, y, 0, 1]
          if multi_design == 3:
            e.rotAngle = input_array[x, y, 0, 2]
      # compare data against threshold value
      if multi_design == 1:
        if e.density < rhomin:
          e.density = rhomin
      else:
        if e.stiff1 < rhomin:
          e.stiff1 = rhomin
        elif e.stiff2 < rhomin:
          e.stiff2 = rhomin
      # assign region    
      if multi_design == 1:
        # are we gray or not?
        if is_gray or (is_color and (input_array[x, y][0] == input_array[x, y][1] and input_array[x, y][1] == input_array[x, y][2])):  
          if e.density >= threshold:
            e.region = 'mech'
            mech_count += 1
          else:
            e.region = 'void'
        elif is_color:
          if input_array[x, y][0] > 0 and input_array[x, y][1] == 0 and input_array[x, y][1] == input_array[x, y][2] == 0:
            e.region = 'red'
            colorful_count += 1
          elif input_array[x, y][0] == 0 and input_array[x, y][1] > 0 and input_array[x, y][1] == input_array[x, y][2] == 0:
            e.region = 'green'
            colorful_count += 1
          elif input_array[x, y][0] == 0 and input_array[x, y][1] == 0 and input_array[x, y][1] == input_array[x, y][2] > 0:
            e.region = 'blue'
            colorful_count += 1
          else:
            e.region = 'colorful'
            colorful_count += 1
      else:
        if e.stiff1 >= threshold or float(e.stiff2) >= threshold:  
          e.region = 'mech'
          mech_count += 1
        else:
          e.region = 'void'
      # assign nodes
      ll = (nx + 1) * y + x  # lowerleft
      e.nodes = ((ll, ll + 1, ll + 1 + nx + 1, ll + nx + 1))
      mesh.elements.append(e)
      # e.dump()
#  if pressure:
#    y = ny - 1
#    for x in range(nx):
#      if (x >= int(0.8 * nx) and y == ny - 1):
#        b = Element()
#        b.type = LINE
#        ll = x
#        b.nodes = ((ll, ll + 1))
#        if e.region == 'mech' or e.region == 'colorful' or e.region == 'red': 
#          b.region = 'pressure2'
#        else:
#          b.region = 'void'
#        mesh.elements.append(b)
#    print 'Warning: pressure area has to be set manually in method create_dense_mesh.'
  mesh.bc.append(("south", range(0, nx + 1)))
  mesh.bc.append(("north", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
  mesh.bc.append(("west", range(0, (nx + 1) * ny + 1, nx + 1)))
  mesh.bc.append(("east", range(nx, (nx + 1) * (ny + 1), nx + 1)))
  if type == 2:
    # pressure boundary for a test case
    print 'Warning: pressure area has to be set manually in method create_dense_mesh.'
    mesh.bc.append(("pressure2", range(int(0.8 * nx), nx + 1)))
  elif type == 3:
    # other boundary conditions for msfem test
    # lower/upper loads
    mid = int((nx+1.)/2.)
    off_x = int(0.05 * nx)
    off_y = int(0.1 * ny)
    if (nx+1)/2. % 2 == 0:
      mesh.bc.append(("load1", range(mid-off_x,mid + off_x)))
      mesh.bc.append(("load2", range((nx+1)*ny + mid-off_x, (nx+1)*ny + mid + off_x)))
    else:
      mesh.bc.append(("load1", range(mid-off_x,mid + 1 + off_x)))
      mesh.bc.append(("load2", range((nx+1)*ny + mid-off_x,(nx+1)*ny + mid + off_x + 1)))
    # support lower left
    mesh.bc.append(("support", range(0,off_x+1)))
    mesh.bc.append(("support", range(0,(nx+1)*off_y+1,nx+1)))
    # support lower right
    mesh.bc.append(("support", range(nx - off_x,nx+1)))
    mesh.bc.append(("support", range(nx,(nx+1)*(off_y+1),nx+1)))
    # support upper left
    mesh.bc.append(("support", range((nx+1)*ny, (nx+1)*ny+off_x+1)))
    mesh.bc.append(("support", range((nx+1)*ny-(nx+1)*off_y,(nx+1)*ny+1,nx+1)))

    # support upper right
    mesh.bc.append(("support", range((nx+1)*(ny+1) -1 - off_x,(nx+1)*(ny+1))))
    mesh.bc.append(("support", range((nx+1)*(ny+1)-(nx+1)*off_y-1,(nx+1)*(ny+1),nx+1)))


  
  # print mesh.bc

  print "dense resolution: " + str(nx) + " x " + str(ny) + " elements (" + str(scale) + "m x " + str(float(ny) / nx * scale) + "m)",
  print " -> " + str(mech_count) + " mech elements out of " + str(nx * ny) + " (" + str(float(mech_count) / (nx * ny) * 100.0) + " %)",
  print " with threshold " + str(threshold) 
  if colorful_count > 0:
    print 'plus ' + str(colorful_count) + ' non gray elements'

# @param mesh dense mesh (input)
# @return sparse mesh
def convert_to_sparse_mesh(dense):
  sparse = Mesh(dense.nx, dense.ny, dense.nz)

  # necessary 0-based nodes as unique set
  nns = set()
  
  # copy element, the indices of the nodes will be replaced later
  for i in range(len(dense.elements)):
    e = dense.elements[i]
    if e.region <> 'void':
      sparse.elements.append(copy.deepcopy(e))
      for n in range(len(e.nodes)):
        nns.add(e.nodes[n])
        
  # nns contains the required nodes uniquely and ordered
  # convert to nnl to access the indices
  nnl = list(nns)
  sparse.nodes = []
  for i in range(len(nnl)):
    sparse.nodes.append(dense.nodes[nnl[i]])
    
  # next we need to correct the element node indices
  # the map is indexed by the dense numbering and contains the sparse indices or -1
  map = len(dense.nodes) * [-1]        
  
  for i in range(len(nnl)):
    map[nnl[i]] = i
      
  # now correct the element nodes
  for e in range(len(sparse.elements)):
     el = sparse.elements[e]
     newnodes = []  # el.nodes is a tuple, values cannot be replaces
     for n in range(len(el.nodes)):
       newnodes.append(map[el.nodes[n]])
       assert(el.nodes[n] <> -1)
     el.nodes = newnodes

  # finally handle the boundary conditions
  sparse.bc = []
  for b in range(len(dense.bc)):
    bc = dense.bc[b]
    dnn = bc[1]  # dense nodes
    nodes = []
    for n in range(len(dnn)):
      if map[dnn[n]] <> -1:
        nodes.append(map[dnn[n]])
    sparse.bc.append((bc[0], nodes))
  return sparse

def count_elements(elements, type):
  count = 0
  for e in elements:
    if e.type == type:
      count += 1
  return count

def write_gid_elements(out, elements, dim):
  for i in range(len(elements)):  # write one based!
    e = elements[i]
    if elem_dim(e.type) == dim:
      nodes = len(e.nodes)
      out.write(str(i + 1) + ' ' + str(e.type) + ' ' + str(nodes_by_type(e.type)) + ' ' + e.region + "\n")
    
      # prepare for second order elements
      for n in range(nodes_by_type(e.type)):
        out.write(str(e.nodes[n] + 1) + ("\n" if n == nodes_by_type(e.type) - 1 else " "))  # write one based node numbers


def write_gid_mesh(mesh, filename,scale = 1):
  # Warning: mesh dimensions should be in [m]
  print "Notification: Make sure that mesh dimensions are in meter [m]!! Otherwise use scale Parameter."
  quad4 = count_elements(mesh.elements, QUAD4)
  hexa8 = count_elements(mesh.elements, HEXA8)
  wedge6 = count_elements(mesh.elements, WEDGE6)
  line = count_elements(mesh.elements, LINE)
  tet4 = count_elements(mesh.elements, TET4)
  tri3 = count_elements(mesh.elements, TRIANGLE3)
  num_1d = line
  num_2d = quad4 + tri3
  num_3d = hexa8 + wedge6 + tet4
  assert(num_1d + num_2d + num_3d == len(mesh.elements))
  print 'number of elements ' + str(num_1d + num_2d + num_3d)
  dim = 3 if num_3d > 0 else 2
  
  out = open(filename, "w")
  
  out.write('[Info]\n')
  out.write('Version 1\n')
  out.write('Dimension ' + str(dim) + '\n')
  out.write('NumNodes ' + str(len(mesh.nodes)) + '\n')
  out.write('Num3DElements ' + str(num_3d) + '\n')
  out.write('Num2DElements ' + str(num_2d) + '\n')
  out.write('Num1DElements ' + str(num_1d) + '\n')
  bcn = 0
  for i in range(len(mesh.bc)):
    bcn += len(mesh.bc[i][1])
  out.write('NumNodeBC ' + str(bcn) + '\n')
  nen = 0
  for i in range(len(mesh.ne)):
    nen += len(mesh.ne[i][1])
  out.write('NumSaveNodes 0\n')
  out.write('NumSaveElements ' + str(nen) + '\n')
  out.write('Num 2d-line      : ' + str(num_1d) + '\n')
  out.write('Num 2d-line,quad : 0\n')
  out.write('Num 3d-line      : 0\n')
  out.write('Num 3d-line,quad : 0\n')
  out.write('Num triangle     : ' + str(tri3) + '\n')
  out.write('Num triangle,quad: 0\n')
  out.write('Num quadr        : ' + str(quad4) + '\n')
  out.write('Num quadr,quad   : 0\n')
  out.write('Num tetra        : ' + str(tet4) + '\n')
  out.write('Num tetra,quad   : 0\n')
  out.write('Num brick        : ' + str(hexa8) + '\n')
  out.write('Num brick,quad   : 0\n')
  out.write('Num pyramid      : 0\n')
  out.write('Num pyramid,quad : 0\n')
  out.write('Num wedge        : ' + str(wedge6) + '\n')
  out.write('Num wedge,quad   : 0\n')
  
  out.write('\n[Nodes]\n')
  out.write('#NodeNr x-coord y-coord z-coord\n')
  for i in range(len(mesh.nodes)):  # write one based!
    out.write(str(i + 1) + "  " + str(mesh.nodes[i][0]/scale) + "  " + str(mesh.nodes[i][1]/scale))
    if dim == 3:
      out.write("  " + str(mesh.nodes[i][2]/scale) + "\n")
    else:
      out.write("  0.0\n")

  out.write('\n[1D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_gid_elements(out, mesh.elements, 1)
  
  out.write('\n[2D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_gid_elements(out, mesh.elements, 2)
 
  out.write('\n[3D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_gid_elements(out, mesh.elements, 3)
  
  out.write('\n[Node BC]\n')
  out.write('#NodeNr Level\n')
  for b in range(len(mesh.bc)):
    bc = mesh.bc[b]
    for n in range(len(bc[1])):
      out.write(str(bc[1][n] + 1) + " " + bc[0] + "\n")
      
  out.write('\n[Save Nodes]\n')
  out.write('#NodeNr Level\n')
  out.write('\n[Save Elements]\n')
  out.write('#ElemNr Level\n')
  for e in range(len(mesh.ne)):
    ne = mesh.ne[e]
    for n in range(len(ne[1])):
      out.write(str(ne[1][n] + 1) + " " + ne[0] + "\n")
  
  out.write("\n \n")
  out.close()
  
  
## creates a 2D mesh of predefined geometry
def create_2d_mesh(type, x_res, y_res, width, opt_height = None, inclusion = None, inclusion_size = None, patch = None):
  
  assert(type == 'bulk2d' or type == 'cantilever2d' or type == 'cantilever2d_reinforced' or type == 'msfem_two_load')
  assert(inclusion == None or inclusion == "rect" or inclusion == "ball")
  assert(inclusion_size == None or inclusion_size <= 2.0)
  
  nx = x_res
  ny = x_res if y_res is None else y_res

  mesh = Mesh(nx, ny)
  
  # buld2d case
  ny = y_res if y_res <> None else x_res
  width = 1.0
  height = float(ny)/nx if opt_height is None else opt_height 
   
  offx = 0.
  offy = 0.
  if type.startswith('cantilever2d'):
    width = 3.0
    height = 2.0
    ny = int(nx * (2. / 3.))
  if type == 'ghost':
    # for virginie
    width = 2.0
    height = 1.0
    nx = nx + 1
    ny = int((nx - 1) * 0.5) + 1
  if type == 'msfem_test':
    width = 1.0
    height = 1.0
    ny = nx

  if type == 'triangle_msfem':
    offx = -1.
    offy = -1.
    width = 2. 
    height = 2.
  if type == 'msfem_two_load':
    width= 2.
    height = 1.
    
     
  dx = width / nx
  dy = height / ny

  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((offx + x * dx, offy + y * dy))
  
  scale = 2 if type == 'triangle_msfem' else 1

 
  # count second region
  second = 0 
 
  # inner to outer boundary interface
  left_iface = []
  right_iface = []
  upper_iface = []
  lower_iface = []
 
  # print mesh.nodes 
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.density = 1.0
      e.type = TRIANGLE3 if type == 'triangle_msfem' else QUAD4
      if type == 'cantilever2d_reinforced' and float(x) >= (28. / 30. * nx):
        e.region = 'reinforce'
      # strange: assure that x is meant to be 2.0 and y is meant to be 1.0 ?!  
      elif type == 'mbb_reinforced' and (x + 1 <= .015 * nx + 1e-5 or x >= 0.985 * nx - 1e-5 or y + 1 <= 0.03 * ny + 1e-5 or y >= 0.97 * ny - 1e-5):
        e.region = 'reinforce'
      elif type == 'ghost':
        if (x == nx - 1 or y == ny - 1):
          e.region = 'ghost'
        else:
          e.region = 'mech'
      elif type == 'msfem_test':
        if False:
          e.region = 'reinforce'
        else:
          e.region = 'mech'
        second += 1
      elif inclusion == 'rect' and x >= nx/2 * (1 - inclusion_size) and x < nx/2 * (1 + inclusion_size) \
                               and y >= ny/2 * (1 - inclusion_size) and y < ny/2 * (1 + inclusion_size):
        e.region = 'inner'
        second += 1        
      elif inclusion == 'ball' and numpy.sqrt((x-nx/2)**2 + (y-ny/2)**2) <= nx*0.5*inclusion_size:
        e.region = 'inner'
        second += 1
      elif patch:
        assert(patch == "3x3")
        e.region = 'reg_' + str((y % ny/3)+1) + '_' + str((x % nx/3)+1)                 
      else:
        e.region = 'mech'
      # assign nodes
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))

      # mark the interface of the inclusion with the outer boundary
      if inclusion <> None:
        if x==0 and e.region == 'inner':
          left_iface.append(ll)
          left_iface.append(ll+nx+1)
        if x == nx-1 and e.region == 'inner':  
          right_iface.append(ll+1)
          right_iface.append(ll+1+nx+1)
        if y==0 and e.region == 'inner':
          lower_iface.append(ll)
          lower_iface.append(ll+1)
        if y==ny-1 and e.region == 'inner':
          upper_iface.append(ll+nx+1)
          upper_iface.append(ll+1+nx+1)
            
      if type == 'triangle_msfem':
        e.nodes = ((ll, ll + 1, ll + 1 + nx + 1))
        e2 = Element()
        e2.density = e.density
        e2.region = e.region
        e2.type = e.type
        e2.nodes = ((ll + 1 + nx + 1, ll + nx + 1, ll))
        mesh.elements.append(e)
        mesh.elements.append(e2)
      else:  
        e.nodes = ((ll, ll + 1, ll + 1 + nx + 1, ll + nx + 1))    
        mesh.elements.append(e)
  if type == 'ghost':
    mesh.bc.append(("south", range(0, nx)))
    mesh.bc.append(("north", range((nx + 1) * (ny - 1), (nx + 1) * (ny - 1) + nx)))
    mesh.bc.append(("west", range(0, (nx + 1) * (ny - 1) + 1, nx + 1)))
    mesh.bc.append(("east", range(nx - 1, (nx + 1) * (ny - 1) + nx, nx + 1)))

    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx - 1]))
    mesh.bc.append(("left_upper", [(nx - 1) * (ny)]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny - 1) + nx - 1]))
    
    mesh.bc.append(("south_ghost", range(0, nx + 1)))
    mesh.bc.append(("north_ghost", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
    mesh.bc.append(("west_ghost", range(0, (nx + 1) * ny + 1, nx + 1)))
    mesh.bc.append(("east_ghost", range(nx, (nx + 1) * (ny + 1), nx + 1)))
  
    mesh.bc.append(("left_lower_ghost", [0]))
    mesh.bc.append(("right_lower_ghost", [nx]))
    mesh.bc.append(("left_upper_ghost", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper_ghost", [(nx + 1) * (ny + 1) - 1]))
    
    mesh.bc.append(("south_mid_ghost", [ int(round(0.5 * (nx)))]))
    mesh.bc.append(("north_mid_ghost", [ int(round(0.5 * ((nx + 1) * (ny + 1) - 1 - (nx + 1) * ny) + (nx + 1) * ny))]))
    mesh.bc.append(("west_mid_ghost", [(nx + 1) * ny - (nx + 1) * int(round(0.5 * ny))]))
    mesh.bc.append(("east_mid_ghost", [(nx + 1) * (ny + 1) - 1 - (nx + 1) * int(round(0.5 * ny))]))
    
    mesh.bc.append(("south_mid", [ int(round(0.5 * (nx - 1)))]))
    mesh.bc.append(("north_mid", [ int(round(0.5 * ((nx + 1) * (ny - 1) + nx - 1 - (nx + 1) * (ny - 1)) + (nx + 1) * (ny - 1)))]))
    mesh.bc.append(("west_mid", [int(round(0.5 * ((nx + 1) * (ny - 1))))]))
    mesh.bc.append(("east_mid", [int(round(0.5 * ((nx + 1) * (ny - 1) + nx - 1 - (nx - 1)) + (nx - 1)))]))
    return mesh, nx, ny, dx, dy
  elif type == 'pressure2':
    mesh.bc.append(("south", range(0, nx + 1)))
    mesh.bc.append(("north", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
    mesh.bc.append(("west", range(0, (nx + 1) * ny + 1, nx + 1)))
    mesh.bc.append(("east", range(nx, (nx + 1) * (ny + 1), nx + 1)))
    mesh.bc.append(("pressure2", range(int(0.8 * nx), nx + 1)))
  
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh
  elif type == 'triangle_msfem':
    mesh.bc.append(("bottom", range(0, nx + 1)))
    mesh.bc.append(("top", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
    mesh.bc.append(("left", range(0, (nx + 1) * ny + 1, nx + 1)))
    mesh.bc.append(("right", range(nx, (nx + 1) * (ny + 1), nx + 1)))
  
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh
  elif type == 'msfem_two_load':
    # lower/upper loads
    mid = int((nx+1.)/2.)
    off_x = int(0.05 * nx)
    off_y = int(0.1 * ny)
    if (nx+1)/2. % 2 == 0:
      mesh.bc.append(("load1", range(mid-off_x,mid + off_x)))
      mesh.bc.append(("load2", range((nx+1)*ny + mid-off_x, (nx+1)*ny + mid + off_x)))
    else:
      mesh.bc.append(("load1", range(mid-off_x,mid + 1 + off_x)))
      mesh.bc.append(("load2", range((nx+1)*ny + mid-off_x,(nx+1)*ny + mid + off_x + 1)))
    # support lower left
    mesh.bc.append(("support", range(0,off_x+1)))
    mesh.bc.append(("support", range(0,(nx+1)*off_y+1,nx+1)))
    # support lower right
    mesh.bc.append(("support", range(nx - off_x,nx+1)))
    mesh.bc.append(("support", range(nx,(nx+1)*(off_y+1),nx+1)))
    # support upper left
    mesh.bc.append(("support", range((nx+1)*ny, (nx+1)*ny+off_x+1)))
    mesh.bc.append(("support", range((nx+1)*ny-(nx+1)*off_y,(nx+1)*ny+1,nx+1)))
  else: 
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx+1)*ny]))
    mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))

  print 'width=' + str(width) + ' height=' + str(height) + ' dx=' + str(dx) + ' dy=' + str(dy)
  
  liu = numpy.unique(left_iface)
  
  if len(left_iface) > 0:
    mesh.bc.append(("left_iface", liu))
    print str(len(liu)) + ' nodes (' + str(len(liu) / float(ny) * 100.) + '%) at the left interface for radius ' + str(inclusion_size)
  if len(right_iface) > 0:
    mesh.bc.append(("right_iface", numpy.unique(right_iface)))
  if len(lower_iface) > 0:
    mesh.bc.append(("lower_iface", numpy.unique(lower_iface)))
  if len(upper_iface) > 0:
    mesh.bc.append(("upper_iface", numpy.unique(upper_iface)))
  
  if second > 0:
    print str(second) + ' elements of secondary region (' + str(100.0 * second / (nx * ny)) + '%)'
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh

  
  else:
    mesh.bc.append(("south", range(0, nx + 1)))
    mesh.bc.append(("north", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
    mesh.bc.append(("west", range(0, (nx + 1) * ny + 1, nx + 1)))
    mesh.bc.append(("east", range(nx, (nx + 1) * (ny + 1), nx + 1)))
  
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh

def create_regular3d_mesh(type, resolution):
  mesh = Mesh()
  
  width = 1.0
  height = 1.0
  
  nx = resolution
  ny = int(nx)
  nz = int(nx)
  
  dx = width / nx
  dy = dx
  dz = dx 
  
  e = 1e-4
  for z in range(nz + 1):
    for y in range(ny + 1):
      for x in range(nx + 1):
        mesh.nodes.append((x * dx, y * dy, z * dz))
 
  # print mesh.nodes
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.type = HEXA8
        e.density = 1.0
         # if type == 'mbb_reinforced' and (float(x) <= (.03 * nx + e) or float(x) >= (.97 * nx - e) or float(y) <= (.03 * ny + e) or float(y) >= (.97 * ny - e)):
        e.region = 'design'
        # assign nodes
        ll = (nx + 1) * (ny + 1) * z + (nx + 1) * y + x  # lowerleft
        e.nodes = ((ll + (nx + 1) * (ny + 1), ll + (nx + 1) * (ny + 1) + nx + 1, ll + (nx + 1) * (ny + 1) + nx + 1 + 1, ll + (nx + 1) * (ny + 1) + 1, ll, ll + nx + 1, ll + nx + 1 + 1, ll + 1))        
        mesh.elements.append(e)  
  return mesh

## creates a mesh of predefined geometry
def create_3d_mesh(type, x_res, y_res, z_res, inclusion, inclusion_size):
  assert(type == "bulk3d" or type == "cantilever3d")
   
  if (type == "bulk3d"): 
      nx = x_res 
      ny = y_res if y_res <> None else x_res 
      nz = z_res if z_res <> None else x_res 
       
      width = 1.0 
      dx = width / nx 
       
      height = float(ny)/nx 
      dy = height / ny 
       
      depth = float(nz)/nx 
      dz = depth / nz 
  elif (type == "cantilever3d"): 
    nx = x_res 
    ny = int(nx * (2./3.)) 
    nz = int(nx * (2./3.)) 
     
    width = 3.0 
    height = 2.0 
    depth = 2.0 
     
    dx = width / nx 
    dy = height / ny 
    dz = depth / nz 

  mesh = Mesh(nx, ny, nz)

  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1
    
  print 'width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz)
  
  
  # the coordinate system in Paraview is a right-hand sided coodrdinate system with z pointing to the viewer 
  #
  #  y ^ 
  #    |
  # z (.)--> x
  # 
  # This are the node numbers if we have only one element. The .mesh file will be transformed to 1-based              
  # x is the fastet variable, z is the slowest variable
  #
  #       2 --------- 3         
  #      /|          /|
  #     / |         / |
  #    6 --------- 7  |
  #    |  |        |  |
  #    |  |        |  |
  #    |  0 -------|- 1
  #    | /         | /
  #    |/          |/
  #    4 --------- 5
  

  for z in range(nnz):  # slowest variable
    for y in range(nny):
      for x in range(nnx):  # fastest variable
        mesh.nodes.append((x * dx, y * dy, z * dz))
  
  # count second region 
  second = 0
 
  # count second region
  second = 0 
  
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.density = 1.0
        e.type = HEXA8
        if inclusion == 'rect' and x >= nnx/2 * (1 - inclusion_size) and x < nnx/2 * (1 + inclusion_size)\
                  and y >= nny/2 * (1 - inclusion_size) and y < nny/2 * (1 + inclusion_size) \
                  and z >= nnz/2 * (1 - inclusion_size) and z < nnz/2 * (1 + inclusion_size) : 
                       e.region = 'inner' 
                       second += 1   
        elif inclusion == 'ball' and numpy.sqrt((x-nnx/2)**2 + (y-nny/2)**2 + (z-nnz/2)**2) <= nnx*inclusion_size: 
            e.region = 'inner' 
            second += 1     
        else: 
            e.region = 'mech' 
   
            # assign nodes 
            # ll = (nx+1)*y*(nx+1) * z + (nx+1) * y + x  # lowerleftfront 
            ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element 
            # start with upper-front-left counterclockwise in the x-z plane. Repeat in then lower plane 
            # e.nodes = ((ll+(nx+1), ll+1+(nx+1), ll+1+(nx+1)+((nx+1)*(ny+1)),ll+(nx+1)+((nx+1)*(ny+1)),ll, ll+1, ll+1+((nx+1)*(ny+1)),ll+((nx+1)*(ny+1))))   
            e.nodes = ((ll+nnx, ll+1+nnx, ll+1+nnx+(nnx*nny),ll+nnx+(nnx*nny),ll, ll+1, ll+1+(nnx*nny),ll+(nnx*nny))) 
              
        mesh.elements.append(e)

  mesh.bc.append(("left", range(0, (nnx * nny * nz) + (nnx * ny) + 1, nnx)))
  mesh.bc.append(("right", range(nx, (nnx * nny * nnz) + 1, nnx)))

  side = (("bottom", []))
  mesh.bc.append(side)
  for z in range(0, nnz):
    for x in range(0, nnx):
      side[1].append((z * nny) * nnx + x)

  side = (("top", []))
  mesh.bc.append(side)
  for z in range(0, nnz):
    for x in range(0, nnx):
      side[1].append((z * nny + ny) * nnx + x)

  
  # back and front as it appears with paraview
  mesh.bc.append(("back", range(0, (nx + 1) * (ny + 1))))
  mesh.bc.append(("front", range(nz * (nx + 1) * (ny + 1), (nz + 1) * (nx + 1) * (ny + 1))))


  mesh.bc.append(("left_bottom_back", [0]))
  mesh.bc.append(("right_bottom_back", [nx]))
  mesh.bc.append(("left_top_back", [nnx * ny]))
  mesh.bc.append(("right_top_back", [nnx * nny - 1]))
  mesh.bc.append(("left_bottom_front", [nnx * nny * nz]))
  mesh.bc.append(("right_bottom_front", [nnx * nny * nz + nx]))
  mesh.bc.append(("left_top_front", [nnx * nny * nz + nnx * ny]))
  mesh.bc.append(("right_top_front", [nnx * nny * nnz - 1]))
  
  if second > 0: 
    print str(second) + ' elements of secondary region (' + str(100.0 * second / (nnx * nny * nnz)) + '%)' 
   
  return mesh

# # LBM pipe_bend and two_inlet_one_outlet example as used by Pingen et al. 2007
# @param case pipe_bend or two_inlet_one_outlet 
def create_lbm2d(resolution, case, inclusion, inclusion_size): 
 
  size = 1.0 
   
  nx = resolution
  ny = nx

  mesh = Mesh(nx, ny)

  
  dx = size / nx
  
  nnx = nx+1 
  nny = ny+1 
  
  eps = 1e-4

  for y in range(nx + 1):
    for x in range(nx + 1):
      mesh.nodes.append((x * dx, y * dx))
 
  # print mesh.nodes 
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.type = QUAD4
      e.density = 1.0
      if inclusion == 'rect' and x >= nnx/2 * (1 - inclusion_size) and x < nnx/2 * (1 + inclusion_size) \
                   and y >= nny/2 * (1 - inclusion_size) and y < nny/2 * (1 + inclusion_size): 
                        e.region = 'obstacle' 
      elif x > 0 and y > 0 and x < nx-1 and y < ny -1:   
        e.region = 'design'
      else:
        e.region = 'boundary'

      # assign nodes
      ll = (nx + 1) * y + x  # lowerleft
      e.nodes = ((ll, ll + 1, ll + 1 + nx + 1, ll + nx + 1))            
      mesh.elements.append(e) 
  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx + 1) * ny]))
  mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))

  if case == 'pipe_bend': 
    mesh.ne.append(('inlet', range(int(0.1 * nx * ny + eps), int(0.3 * nx * ny + nx + eps), nx)))
    mesh.ne.append(('outlet', range(int((ny - 1) * nx + 0.7 * nx - eps), int((ny - 1) * nx + 0.9 * nx + eps))))
  elif case == 'two_inlet_one_outlet':
    mesh.ne.append(('inlet', range(int((0.25 - 1. / 16) * nx * ny + eps), int((0.25 + 1. / 16) * nx * ny + nx + eps), nx)))
    mesh.ne.append(('inlet', range(int((0.75 - 1. / 16) * nx * ny + eps), int((0.75 + 1. / 16) * nx * ny + nx + eps), nx)))
    mesh.ne.append(('outlet', range(int(0.375 * nx * ny - 1 + eps), int(0.625 * nx * ny - 1 + eps), nx)))
  elif case == 'two_inlet_two_outlet':
    inletLength = 0.15 * ny
    mesh.ne.append(('inlet',range(int(0.2*nx*ny+eps), int(0.2*nx*ny+eps + inletLength*nx),nx) ))
    mesh.ne.append(('inlet',range(int(0.8*nx*ny+eps-inletLength*nx), int(0.8*nx*ny+eps),nx) ))
    mesh.ne.append(('outlet',range(int(0.2*nx*ny+eps + nx-1), int(0.2*nx*ny+eps + inletLength*nx+ nx-1),nx) ))
    mesh.ne.append(('outlet',range(int(0.8*nx*ny+eps-inletLength*nx+ nx-1), int(0.8*nx*ny+eps+ nx-1),nx) ))
    
  elif case == 'two_inlet_two_outlet': 
    inletLength = 0.15 * ny 
    mesh.ne.append(('inlet',range(int(0.2*nx*ny+eps), int(0.2*nx*ny+eps + inletLength*nx),nx) )) 
    mesh.ne.append(('inlet',range(int(0.8*nx*ny+eps-inletLength*nx), int(0.8*nx*ny+eps),nx) )) 
    mesh.ne.append(('outlet',range(int(0.2*nx*ny+eps + nx-1), int(0.2*nx*ny+eps + inletLength*nx+ nx-1),nx) )) 
    mesh.ne.append(('outlet',range(int(0.8*nx*ny+eps-inletLength*nx+ nx-1), int(0.8*nx*ny+eps+ nx-1),nx) )) 
      
  return mesh

def create_backstep(x_res, y_res, z_res): 
  mesh = Mesh() 
 
  nx = x_res 
  ny = y_res if y_res <> None else x_res 
  nz = z_res if z_res <> None else x_res 
 
  nnx = nx+1 
  nny = ny+1 
  nnz = nz+1 
   
  width = 1.0 
  dx = width / nx 
   
  height = float(ny)/nx 
  dy = height / ny 
   
  depth = float(nz)/nx 
  dz = depth / nz 
   
  eps = 1e-4 
 
  print 'width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz) 
   
   
  # the coordinate system in Paraview is a right-hand sided coodrdinate system with z pointing to the viewer  
  # 
  #  y ^  
  #    | 
  # z (.)--> x 
  #  
  # This are the node numbers if we have only one element. The .mesh file will be transformed to 1-based               
  # x is the fastet variable, z is the slowest variable 
  # 
  #       2 --------- 3          
  #      /|          /| 
  #     / |         / | 
  #    6 --------- 7  | 
  #    |  |        |  | 
  #    |  |        |  | 
  #    |  0 -------|- 1 
  #    | /         | / 
  #    |/          |/ 
  #    4 --------- 5 
   
 
  for z in range(nnz): # slowest variable 
    for y in range(nny): 
      for x in range(nnx): # fastest variable 
        mesh.nodes.append((x * dx, y * dy, z * dz)) 
  
  second = 0 
   
  # store ymax and of zmax obstacle 
  # need this information to set inlet and outlet elements 
  ymax = 0 
 
  for z in range(nz): 
    for y in range(ny): 
      for x in range(nx): 
        e = Element() 
        e.density = 1.0 
        e.type = HEXA8 
        if z > 0 and z < nz - 1 and x > 0 and y > 0 and x < int(0.3*nx) and y < int(0.35*ny):   
          e.region = 'obstacle' 
          second += 1 
          ymax = max(ymax,y) 
        elif x > 0 and y > 0  and z > 0 and x < nx-1 and y < ny -1 and z < nz - 1:  
          e.region = 'design' 
        else: 
          e.region = 'boundary' 
     
        # assign nodes 
        ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element 
        e.nodes = ((ll+nnx, ll+1+nnx, ll+1+nnx+(nnx*nny),ll+nnx+(nnx*nny),ll, ll+1, ll+1+(nnx*nny),ll+(nnx*nny))) 
         
        mesh.elements.append(e) 
 
  print "Created " + str(second) + " obstacle elements" 
           
  mesh.bc.append(("left", range(0, (nnx*nny*z)+(nnx*ny)+1, nnx))) 
  mesh.bc.append(("right", range(nx, (nnx*nny*nnz)+1, nnx))) 
   
  side = (("bottom", [])) 
  mesh.bc.append(side) 
  for z in range(0, nnz): 
    for x in range(0, nnx): 
      side[1].append((z*nny)*nnx+x) 
 
  side = (("top", [])) 
  mesh.bc.append(side) 
  for z in range(0, nnz): 
    for x in range(0, nnx): 
      side[1].append((z*nny+ny)*nnx+x) 
 
   
  # back and front as it appears with paraview 
  mesh.bc.append(("back", range(0, (nx+1)*(ny+1)))) 
  mesh.bc.append(("front", range(nz*(nx+1)*(ny+1), (nz+1)*(nx+1)*(ny+1)))) 
 
 
  mesh.bc.append(("left_bottom_back",   [0])) 
  mesh.bc.append(("right_bottom_back",  [nx])) 
  mesh.bc.append(("left_top_back",      [nnx*ny])) 
  mesh.bc.append(("right_top_back",     [nnx*nny-1])) 
  mesh.bc.append(("left_bottom_front",  [nnx*nny*nz])) 
  mesh.bc.append(("right_bottom_front", [nnx*nny*nz+nx])) 
  mesh.bc.append(("left_top_front",     [nnx*nny*nz+nnx*ny])) 
  mesh.bc.append(("right_top_front",    [nnx*nny*nnz-1])) 
   
  for i in range(1,nz-1): 
    mesh.ne.append(('inlet',range(ymax*nx+nx+i*nx*ny,nx*ny*(i+1)-nx,nx))) 
    mesh.ne.append(('outlet',range(2*nx+i*nx*ny-1,nx*ny*(i+1)-nx,nx))) 
 
  return mesh 
 
   
def create_lbm3d(x_res, y_res, z_res, case, inclusion, inclusion_size): 
 
  nx = x_res 
  ny = y_res if y_res <> None else x_res 
  nz = z_res if z_res <> None else x_res 
 
  nnx = nx+1 
  nny = ny+1 
  nnz = nz+1 
   
  width = 1.0 
  dx = width / nx 
   
  height = float(ny)/nx 
  dy = height / ny 
   
  depth = float(nz)/nx 
  dz = depth / nz 

  mesh = Mesh(nx,ny,nz)    
  eps = 1e-4 
 
  print 'width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz) 
   
   
  # the coordinate system in Paraview is a right-hand sided coodrdinate system with z pointing to the viewer  
  # 
  #  y ^  
  #    | 
  # z (.)--> x 
  #  
  # This are the node numbers if we have only one element. The .mesh file will be transformed to 1-based               
  # x is the fastet variable, z is the slowest variable 
  # 
  #       2 --------- 3          
  #      /|          /| 
  #     / |         / | 
  #    6 --------- 7  | 
  #    |  |        |  | 
  #    |  |        |  | 
  #    |  0 -------|- 1 
  #    | /         | / 
  #    |/          |/ 
  #    4 --------- 5 
   
 
  for z in range(nnz): # slowest variable 
    for y in range(nny): 
      for x in range(nnx): # fastest variable 
        mesh.nodes.append((x * dx, y * dy, z * dz)) 
  
  second = 0 
  for z in range(nz): 
    for y in range(ny): 
      for x in range(nx): 
        e = Element() 
        e.density = 1.0 
        e.type = HEXA8 
        if inclusion == 'rect' and x >= nnx/4 * (1 - inclusion_size) and x < nnx/4 * (1 + inclusion_size) \
                   and y >= nny/2 * (1 - inclusion_size) and y < nny/2 * (1 + inclusion_size)\
                   and z >= nnz/2 * (1 - inclusion_size) and z < nnz/2 * (1 + inclusion_size): 
                        e.region = 'obstacle' 
                        second += 1 
        elif x > 0 and y > 0  and z > 0 and x < nx - 1 and y < ny -1 and z < nz - 1:  
          e.region = 'design' 
        else: 
          e.region = 'boundary' 
     
        # assign nodes 
        # ll = (nx+1)*y*(nx+1) * z + (nx+1) * y + x  # lowerleftfront 
        ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element 
        # start with upper-front-left counterclockwise in the x-z plane. Repeat in then lower plane 
        # e.nodes = ((ll+(nx+1), ll+1+(nx+1), ll+1+(nx+1)+((nx+1)*(ny+1)),ll+(nx+1)+((nx+1)*(ny+1)),ll, ll+1, ll+1+((nx+1)*(ny+1)),ll+((nx+1)*(ny+1))))   
        e.nodes = ((ll+nnx, ll+1+nnx, ll+1+nnx+(nnx*nny),ll+nnx+(nnx*nny),ll, ll+1, ll+1+(nnx*nny),ll+(nnx*nny))) 
         
        mesh.elements.append(e) 
 
  #print "Created " + str(second) + " obstacle elements" 
           
  mesh.bc.append(("left", range(0, (nnx*nny*z)+(nnx*ny)+1, nnx))) 
  mesh.bc.append(("right", range(nx, (nnx*nny*nnz)+1, nnx))) 
   
  side = (("bottom", [])) 
  mesh.bc.append(side) 
  for z in range(0, nnz): 
    for x in range(0, nnx): 
      side[1].append((z*nny)*nnx+x) 
 
  side = (("top", [])) 
  mesh.bc.append(side) 
  for z in range(0, nnz): 
    for x in range(0, nnx): 
      side[1].append((z*nny+ny)*nnx+x) 
 
   
   
  if case == 'pipe_bend': 
    area = 0.04 
    in_x = int(math.sqrt(area / (dx * dx))+eps) ## find out how many elements in one direction 
    dist_wall = 0.04 
    x_wall = int(math.sqrt(dist_wall / (dx * dx))+eps) #find out distance from wall 
    for i in range(nx*ny*nz-(x_wall+in_x)*nx*ny,nx*ny*nz-x_wall*nx*ny,nx*ny): 
      mesh.ne.append(('inlet',range(i+x_wall*nx,i+(x_wall+in_x)*nx,nx))) 
    for i in range(0,in_x*nx*ny,nx*ny): 
      mesh.ne.append(('outlet', range((x_wall+1)*nx*ny-in_x-x_wall+i,(x_wall+1)*nx*ny-x_wall+i,1))) 
  elif case == 'pipe': 
    for i in range (1,nz-1): 
      mesh.ne.append(('inlet', range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx))) 
      mesh.ne.append(('outlet', range(int(i*nx*ny + 2*nx-1),int(i*nx*ny + 2*nx-1+nx*(ny-2)),nx))) 
  elif case == 'two_inlet_one_outlet': 
    print "Not implemented yet!" 
  elif case == 'diffuser': 
    for i in range (1,nz-1): 
      mesh.ne.append(('inlet', range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx))) 
    for i in range (int(0.3*nz),int(0.7*nz)): 
      mesh.ne.append(('outlet', range(int(i*nx*ny + 0.3*nx*ny + nx-1),int(i*nx*ny + 0.7*nx*ny+ nx-1),nx))) 
  elif case == 'distributor': 
    center_x = nx / 2.0 
    center_y = ny / 2.0 
    center_z = nz / 2.0 
    width_x = 0.1 * nx  
    width_y = 0.1 * ny  
    width_z = 0.1 * nz  
    for i in range(int(round(center_z-width_z/2.0)),int(round(center_z+width_z/2.0))): 
      mesh.ne.append(('outlet',range(int(i*nx*nz-center_x+2.0*width_x),int(i*nx*nz-center_x+3.0*width_x),1))) # top face 
      mesh.ne.append(('outlet',range(int((i-1)*nx*nz+center_x+2.0*width_x),int((i-1)*nx*nz+center_x+3.0*width_x),1))) # bottom face 
    for i in range(int(center_z-int(width_z)),int(center_z+int(width_z))): 
      mesh.ne.append(('inlet',range(int(i*nx*ny+nx*int(center_y-width_y)),int(i*nx*ny+nx*int(center_y+width_y)),nx))) #left face 
    for i in range(int(round(center_y-width_y/2.0)),int(round(center_y+width_y/2.0)),1): 
      mesh.ne.append(('outlet',range(int(i*nx+center_x+2.0*width_x),int(i*nx+center_x+3.0*width_x),1))) #back face 
      mesh.ne.append(('outlet',range(int(nx*ny*(nz-1)+i*nx+center_x+2.0*width_x),int(nx*ny*(nz-1)+i*nx+center_x+3.0*width_x),1))) #front face 
   
  return mesh 

# creates a mesh from hdf5 file  
def create_mesh_from_hdf5(hdf5_f, region, bcregions, region_force=None, region_support=None, threshold=0.):
  hdf5_file = h5py.File(hdf5_f, 'r')
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'].value  # for all regions
  # assume that region[0] is design, region[1] is non-design or void 
  reg_elements_region = []
  for i in range(len(region)):
    reg_elements_region.append(hdf5_file['/Mesh/Regions/' + region[i] + '/Elements'].value)
    
  types = hdf5_file['/Mesh/Elements/Types'].value
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates'].value
  length = len(hdf5_file['/Mesh/Regions/' + region[0] + '/Nodes'].value)
  #reg_nodes = [[0 for col in range(len(region))] for row in range(length)]
  #for i in range(len(region)):
  #  reg_nodes[i][:] = hdf5_file['/Mesh/Regions/' + region[i] + '/Nodes']
  #design_var = hdf5_file['/Results/Mesh/MultiStep_1/Step_0/physicalPseudoDensity/mech/Elements/Real'].value
    
  # Create mesh  
  mesh = Mesh()
  # extract boundary force nodes from region_force if available
  if region_force <> None:
    reg_force_nodes = hdf5_file['/Mesh/Groups/' + region_force + '/Nodes']
    mesh.bc.append((region_force, reg_force_nodes[:] - 1))
  # extract boundary force nodes from region_force if available
  elif region_support <> None:
    reg_support_nodes = hdf5_file['/Mesh/Groups/' + region_support + '/Nodes']
    mesh.bc.append((region_support, reg_support_nodes[:] - 1))
  else:
    #array of boundary regions must be given, e.g. ['support','load1','load2']
    for i in range(len(bcregions)):
      bc_nodes = hdf5_file['Mesh/Groups/' + str(bcregions[i]) + '/Nodes']
      mesh.bc.append((bcregions[i], bc_nodes[:] - 1))
  # insert nodes
  for i in range(len(all_nodes)):
    mesh.nodes.append(all_nodes[i])
    
  # counter for regions
  idx = range(len(region))
  for i in range(len(region)):
    idx[i] = 0  
  for i in range(len(all_elements[:, 0])):
    e = Element()
    e.nodes = (all_elements[i, :] - 1)
    #e.density = design_var[i]
    for j in range(len(region)):
      if idx[j] < len(reg_elements_region[j]):
        if i + 1 == reg_elements_region[j][idx[j]]:
          #if e.density >= threshold:
          e.region = region[j]
          #else:
          #  e.region = 'void'
          idx[j] += 1
    e.type = mesh_type_from_hdf5(types[i])
    mesh.elements.append(e) 
  return mesh


def create_mesh_from_tetgen(meshfile, region):
  print 'read tetgenfile' + meshfile + '.1.ele'
  all_elements = numpy.loadtxt(meshfile + '.1.ele', dtype='int' , skiprows=1)
  print 'read all_elements done'
  all_nodes = numpy.loadtxt(meshfile + '.1.node', skiprows=1)
  print 'read all_nodes done'
  # all_faces = numpy.loadtxt(meshfile+'1.face',skiprows=1)
  # all_edges = numpy.loadtxt(meshfile+'1.edge',skiprows=1)
  
    
  # Create mesh 3D Tetrahedron  
  mesh = Mesh()  
  for i in range(len(all_nodes)):
    mesh.nodes.append(all_nodes[i, 1:])  
  for i in range(len(all_elements[:, 0])):
    e = Element()
    e.nodes = (all_elements[i, 1:] - 1)
    e.density = 1.
    e.region = region
    e.type = TET4
    mesh.elements.append(e) 
  return mesh


def create_validation_mesh(coords,nondes_coords, s1, s2, s3, ip_nx, grad, dir, scale,d_f,valid_position,type = "apod6",thres=0.0,csize = None,simp = None):
  centers, mi, ma = coords[0:3]  # design elements
  nondes_centers, nondes_min, nondes_max = nondes_coords[0:3]  # nondesign elements
  mesh = Mesh()
  # number of cells in one-direction per coarse cell (for validation)
  dx_f,dy_f,dz_f = d_f
  if scale <= 0:
    scale = 1.0
  
  # determine cell size of coarse grid for detection of geometry
  if csize is None:
    dx = (ma[0] - mi[0]) / ip_nx
    dy = dx
    dz = dx
  else:
    dx = csize[0]
    dy = csize[1]
    dz = csize[2]
  # calculate convex hull of non-design and design nodes
  hull = Delaunay(nondes_centers)
  if type == "robot":
    hull_des = Delaunay(centers)
  print 'calculating convex hull of non-design done'
  # choose validation for simp result or 2-scale result
  if simp is None:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1, s2, s3, dx,dy,dz)
  else:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1,None,None, dx,dy,dz)
  print 'interpolation of thicknesses done'

  # lowest density = void density
  void = min(ip_data.ravel())
  
  # add points to fine mesh including shell
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  # where we want nodes
  nx = (int(delta[0] / dx) + 1)*dx_f
  ny = (int(delta[1] / dy) + 1)*dy_f
  nz = (int(delta[2] / dz) + 1)*dz_f

  print "nx,ny,nz = " + str(nx) + ", " + str(ny) + ", " + str(nz)
  print "des: ma,mi = " + str(ma) + " " + str(mi) 
  print "nondes: ma,mi = " + str(nondes_max) + " " + str(nondes_min) 
  #thickness of shell 1mm: tx,... is thickness of non-design shell
  if type == "apod6":
    tx = 0
    ty = int(dy_f / (dy*1000))+1
    tz = 0
  elif type == "robot":
    tx = int(dx_f / dx)
    ty = int(dy_f / dy)
    tz = int(dz_f / dz)
    tx *= 5
    ty *= 5
    tz *= 5 
 
  # offset for function apod6 (valid_position), fixes a bug 
  offset = dx + 1e-6
  if ny == 0 or nz == 0 or nx == 0:
    print 'chose a higher hom_samples or smaller cell_size such that also the smallest side gets discretized'
    exit()
  print "tx, ty, tz = " + str(tx) + ", " + str(ty) + ", " + str(tz)
  # add nodes including offset for shell 
  for z in range(-tz,nz+1+tz):
    for y in range(-ty, ny + 1 + ty):
      for x in range(-tx,nx+1+tx):
        mesh.nodes.append((mi[0] + float(x) * dx/dx_f, mi[1] +  float(y) * dy/dy_f, mi[2] +float(z) * dz/dz_f))
  print 'inserting mesh.nodes done'
  nnx = nx + 2 * tx
  nny = ny + 2 * ty 
  nnz = nz + 2 * tz
  array = -1. * numpy.ones((nnx,nny,nnz))
  res = [dx_f,dy_f,dz_f]
  count = 0
  hole = -2. if type == "robot" else void
  for k in xrange(tz,nz-dz_f + 1 + tz,dz_f):
    for j in xrange(ty,ny-dy_f + 1 + ty,dy_f):
      for i in xrange(tx,nx -dx_f + 1 + tx,dx_f):    
        coord = out[count]
        if simp is None:
          s1, s2, s3 = ip_data[count][0:3]
        else:
          s1 = ip_data[count][0]
        l = [i,j,k]
        u = [i+dx_f,j+dy_f,k+dz_f]
        # if s1 < 0 point is out of the convex hull
        if s1 > 0.0 and simp is None:
          # 2scale optimization
          if s1 >= thres or s2 >= thres or s3 >= thres:
            if not valid_position(coord, coords,offset):
              create_cross_3D(array,l,u,void,void,void,void,res)
            else:
              create_cross_3D(array,l,u,s1,s2,s3,hole,res)
          else:
              create_cross_3D(array,l,u,void,void,void,hole,res)
        elif simp is None:
          create_cross_3D(array,l,u,-1.,-1.,-1.,hole,res)#void,void,void,hole,res)
        else:
          # simp
          if s1 > 0:
            if s1 >= thres:
              if not valid_position(coord, coords,offset):
                  array[l[0]:u[0],l[1]:u[1],l[2]:u[2]] = void * numpy.ones((res[0], res[1], res[2]))
              else:
                  array[l[0]:u[0],l[1]:u[1],l[2]:u[2]] = numpy.ones((res[0], res[1], res[2]))
            else:
                array[l[0]:u[0],l[1]:u[1],l[2]:u[2]] = void * numpy.ones((res[0], res[1], res[2]))
          else:
              array[l[0]:u[0],l[1]:u[1],l[2]:u[2]] = void * numpy.ones((res[0], res[1], res[2]))             
        count += 1
  print 'calculation of density array done'
  number = 0
  for z in range(nnz):
    for y in range(nny):
      for x in range(nnx):
        e = Element()
        e.type = HEXA8
        ll = (nnx + 1) * (nny + 1) * z + (nnx + 1) * y + x  # lowerleft
        e.nodes = ((ll + (nnx + 1) * (nny + 1), ll + (nnx + 1) * (nny + 1) + nnx + 1, ll + (nnx + 1) * (nny + 1) + nnx + 1 + 1, ll + (nnx + 1) * (nny + 1) + 1, ll, ll + nnx + 1, ll + nnx + 1 + 1, ll + 1))        
        condition = True #if type == "robot" else ((x < tx) or (y < ty) or (z < tz) or (x >= nx+tx) or (y >= ny+ty) or (z >= nz+tz))
        count = 0
        for i in range(len(e.nodes)):
	  node = mesh.nodes[e.nodes[i]]
          # test if element is above or below design region, bounds are given by non-design region
          if (node[1] >= -0.33403 - ((0.5*dy)/dy_f) and node[1] < -0.33303 + ((0.5*dy)/dy_f)) or (node[1] <= -0.35203 + ((0.5*dy)/dy_f) and node[1] > -0.35303 - ((0.5*dy)/dy_f)):
	    count +=1
        # calculate center of element
        center = numpy.array([0.0, 0.0, 0.0])
        len_nod = len(e.nodes)
        for n in range(len_nod):
          center += mesh.nodes[e.nodes[n]]
        center *= 1.0 / len_nod
        if count >= 5:
          # test if is in convex hull of non-design nodes
          if in_hull(center, hull):
            if not valid_position(center, coords,offset):
              e.region = 'void1'
            #elif type == "robot" and array[x,y,z] > 0.9:
            else:
		e.region = 'non-design'
 	        number += 1
            #elif type == "robot" and not in_hull(center, hull_des):
	    #  e.region = 'non-design'
	    #  number += 1.
            #else:
            #  if type == "robot":
            #    e.region = 'void3'
	    #  else:
		#e.region = "non-design"
		#number += 1
          else:
            e.region = 'void2'
        else:
          # test if is in convex hull of non-design nodes
          if in_hull(center, hull):
	    if not valid_position(center, coords,offset):
              e.region = 'void1'
            #elif type == "robot" and array[x,y,z] > 0.9:
            elif array[x,y,z] > 0.9:
              number += 1
              e.region = 'design'
            elif array[x,y,z] <= void + 1e-5:
              e.region = 'void3'
	    else:
	      e.region = "void5"
          else:
	    e.region = "void4"
        #elif array[x,y,z] <= void:
        #  e.region = 'void3'
        #else:
        #  number += 1
        #  e.region = 'design'  
        mesh.elements.append(e)
  if type == "apod6":
    # add apod6 boundary conditions to mesh      
    mesh = add_apod6_boundary_conditions(mesh)
  elif type == "robot":
    mesh = add_robot_boundary_conditions(mesh)
  print 'mesh has ' + str(number) + "design and non-design elements"
  return mesh

def in_hull(p, hull,to = None):
  # Test if points in `p` are in `hull`
  #`p` should be a `NxK` coordinates of `N` points in `K` dimensions
  #`hull` is either a scipy.spatial.Delaunay object or the `MxK` array of the 
  # coordinates of `M` points in `K`dimensions for which Delaunay triangulation
  # will be computed
  #from scipy.spatial import Delaunay
  #if not isinstance(hull,Delaunay):
  #  hull = Delaunay(hull)
  return hull.find_simplex(p,tol = to)>=0    

def create_cross_3D(array,l,u,s1,s2,s3,void,res):
  # creates 3D cross in array for fine mesh generation; validation of optimal result by FEM
  # l, u are the upper and lower bounds for the subdomain, e.g [lx,ly,lz] [ux,uy,uz]
  # s1,s2,s3 are the cross thicknesses of one cross
  # res is the discretization resolution for each cross, e.g. [resx,resy,resz]
  # array is the density array
  
  array[l[0]:u[0],l[1]:u[1],l[2]:u[2]] = void * numpy.ones((res[0], res[1], res[2]))
  offx = int((res[0] / 2.) * (1. - s1) + 0.5)
  offy = int((res[1] / 2.) * (1. - s2) + 0.5)
  offz = int((res[2] / 2.) * (1. - s3) + 0.5)
  for i in range(0, res[0]):
    for j in range(offx, res[1] - offx):
      for k in range(offx, res[2] - offx):
        array[l[0] + i,l[1] + j, l[2] + k] = 1.
  for i in range(offy, res[0] - offy):
    for j in range(0, res[1]):
      for k in range(offy, res[2] - offy):
        array[l[0] + i,l[1] + j,l[2] + k] = 1.
  for i in range(offz, res[0] - offz):
    for j in range(offz, res[1] - offz):
      for k in range(0, res[2]):
        array[l[0] + i,l[1] + j,l[2] + k] = 1.
        
def add_robot_boundary_conditions(mesh):
  # add loads and support to robot arm mesh
  # support
  m1 = [-107.44,-54.,0.]
  m2 = [-147.44,-54.,40.]
  m3 = [-187.44,-54.,0.]
  m4 = [-147.44,-54.,-40.]
  
  #loads
  m5 = [213.,0.,0.]
  m6 = [250.,0.,37.]
  m7 = [287.,0.,0.]
  m8 = [250.,0.,-37.]
  
  r = 5
  force1 = []
  support = []
  delta = 1.
  delta_y = 2.9
  for i in range(len(mesh.nodes)):
    coord = mesh.nodes[i][:]
    if abs(coord[1] + 54.) < delta_y and (coord[0] - m1[0]) ** 2 + (coord[2] - m1[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m1[0]) ** 2 + (coord[2] - m1[2]) ** 2 < (r+delta) ** 2:
      support.append(i)
    elif abs(coord[1] + 54.) < delta_y and (coord[0] - m2[0]) ** 2 + (coord[2] - m2[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m2[0]) ** 2 + (coord[2] - m2[2]) ** 2 < (r+delta) ** 2:
      support.append(i)
    elif abs(coord[1] + 54.) < delta_y and (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 < (r+delta) ** 2:
      support.append(i)
    elif abs(coord[1] + 54.) < delta_y and (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 < (r+delta) ** 2:
      support.append(i)
    elif abs(coord[1] + 0.) < delta_y and (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 < (r+delta) ** 2:
      force1.append(i)
    elif abs(coord[1] + 0.) < delta_y and (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 < (r+delta) ** 2:
      force1.append(i)
    elif abs(coord[1] + 0.) < delta_y and (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 < (r+delta) ** 2:
      force1.append(i)
    elif abs(coord[1] + 0.) < delta_y and (coord[0] - m8[0]) ** 2 + (coord[2] - m8[2]) ** 2 >= (r-delta) ** 2 and (coord[0] - m8[0]) ** 2 + (coord[2] - m8[2]) ** 2 < (r+delta) ** 2:
      force1.append(i)
      
  mesh.bc.append(('force1', force1))
  mesh.bc.append(('support', support))
  return mesh
            
def add_apod6_boundary_conditions(mesh):
  # add apod6 boundary conditions
  m1 = [33.052, -0.353, -2.474]
  m2 = [33.046, -0.353, -2.518]
  m3 = [33.131, -0.353, -2.449]
  m4 = [33.124, -0.353, -2.498]
  m5 = [32.978, -0.353, -2.436]
  m6 = [32.971, -0.353, -2.485]
  m7 = [33.023, -0.353, -2.559]
  m8 = [33.042, -0.353, -2.548]
  m9 = [33.004, -0.353, -2.443]
  r1 = 0.0195
  r2 = 0.0165
  r3 = 0.0058
  delta = 0.0002
  force1 = []
  force2 = []
  force3 = []
  support = []
  support2 = []
  support3 = []
  upper_bound = -0.35303
  dy = 0.0009
  for i in range(len(mesh.nodes)):
    coord = mesh.nodes[i][:]
    if abs(coord[1] - upper_bound) < dy and (coord[0] - m1[0]) ** 2 + (coord[2] - m1[2]) ** 2 < r1 ** 2:
      force1.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m2[0]) ** 2 + (coord[2] - m2[2]) ** 2 < r2 ** 2:
      force2.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m9[0]) ** 2 + (coord[2] - m9[2]) ** 2 < r1 ** 2:
      force3.append(i)
    elif abs(coord[1] - upper_bound) <  dy  and (coord[0] - m8[0]) ** 2 + (coord[2] - m8[2]) ** 2 < (r3+delta) ** 2:
      support2.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 < (r3+delta) ** 2:
      support3.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 < (r3+delta) ** 2:
      support3.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 < (r3+delta) ** 2:
      support3.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 < (r3+delta) ** 2:
      support3.append(i)
    elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 < (r3+delta) ** 2:
      support3.append(i)
    elif (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 < r3 ** 2:
      support.append(i)
    elif (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 < r3 ** 2:
      support.append(i)
    elif (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 < r3 ** 2:
      support.append(i)
    elif (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 < r3 ** 2:
      support.append(i)
    elif (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 < r3 ** 2:
      support.append(i)
          
  mesh.bc.append(('force1', force1))
  mesh.bc.append(('force2', force2))
  mesh.bc.append(('force3', force3))
  mesh.bc.append(('support', support))
  mesh.bc.append(('support2', support2))
  mesh.bc.append(('support3', support3))
  return mesh
  

def create_mesh_for_apod6(meshfile, all_nodes = [], elements = [], force1 = [], force2 = [],force3 = [], support = [], support2 = [], support3 =[]):
  # create element and nodes files by hand from optistruct
  if len(all_nodes) == 0:
    # load files from optistruct file
    hexa_elements = numpy.loadtxt(meshfile + '.hexa.elements', dtype='int', skiprows=1)
    wedge_elements = numpy.loadtxt(meshfile + '.wedge.elements', dtype='int', skiprows=1)
    all_nodes = numpy.loadtxt(meshfile + '.nodes', skiprows=1)
    # Rotate nodes for apod6
    ay = -0.084636333418591
    Ry = numpy.matrix(((math.cos(ay), 0., math.sin(ay)), (0., 1., 0.), (-math.sin(ay), 0., math.cos(ay))))

  # Create mesh  
  # add nodes    
  mesh = Mesh()  
  for i in range(len(all_nodes)):
    coord = numpy.matrix(((all_nodes[i, 1]), (all_nodes[i, 2]), all_nodes[i, 3])).T
    # print Rx
    # print coord  
    new_coord = coord#Ry * coord
    # new_coord = Rz * new_coord
    mesh.nodes.append([new_coord[0, 0], new_coord[1, 0], new_coord[2, 0]])
  if len(elements) == 0:  
    # hexaeder     
    for i in range(len(hexa_elements[:, 0])):
      e = Element()
      e.nodes = (hexa_elements[i, 1:] - 1)
      e.density = 1.
      shell = 0
      for j in range(8):
        coord = mesh.nodes[e.nodes[j]]
        if abs(coord[1] + 353.034) < 1.1 or abs(coord[1] + 333.034) < 1.1:
          shell += 1
      if shell > 4:
        e.region = 'non-design'
      else:
        e.region = 'design'
      e.type = HEXA8
      mesh.elements.append(e)
    
    # wedge elements  
    for i in range(len(wedge_elements[:, 0])):
      e = Element()
      e.nodes = (wedge_elements[i, 1:] - 1)
      e.density = 1.
      shell = 0
      for j in range(6):
        coord = mesh.nodes[e.nodes[j]]
        if abs(coord[1] + 353.034) < 1.1 or abs(coord[1] + 333.034) < 1.1:
          shell += 1
      if shell > 3:
        e.region = 'non-design'
      else:
        e.region = 'design'
      e.type = WEDGE6
      mesh.elements.append(e)
  else:
    for i in range(len(elements)):
      e = Element()
      e.nodes = (elements[i][1:])
      for k in range (len(e.nodes)):
        e.nodes[k] -= 1
      e.density = 1.
      shell = 0
      for j in range(len(e.nodes)):
        coord = mesh.nodes[e.nodes[j]]
        if abs(coord[1] + 353.034) < 1.1 or abs(coord[1] + 333.034) < 1.1:
          shell += 1
        if (len(e.nodes) == 6 and shell > 3) or shell > 4:
          e.region = 'non-design'
        else:
          e.region = 'design'
        if len(e.nodes) == 4:
          e.type = TET4
        elif len(e.nodes) == 6:
          e.type = WEDGE6
        elif len(e.nodes) == 8:
          e.type = HEXA8
      mesh.elements.append(e)
  # include boundary conditions for 442.mesh manually
  if len(force1) == 0 and len(force2) == 0:
    m1 = [33052., -353., -2474.]
    m2 = [33046., -353., -2518.]
    m3 = [33131., -353., -2449.]
    m4 = [33124., -353., -2498.]
    m5 = [32978., -353., -2436.]
    m6 = [32971., -353., -2485.]
    m7 = [33023., -353., -2559.]
    m8 = [33042., -353., -2548.]
    m9 = [33004., -353., -2443.]
    
    r1 = 15.5#19.5
    r2 = 12.5#16.5
    r3 = 5.8
    delta = 0.2
    force1 = []
    force2 = []
    force3 = []
    support = []
    support2 = []
    support3 = []
    for i in range(len(all_nodes)):
      coord = all_nodes[i, 1:]
      if (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 < r3 ** 2:
        support.append(i)
      elif (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 < r3 ** 2:
        support.append(i)
      elif (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 < r3 ** 2:
        support.append(i)
      elif (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 < r3 ** 2:
        support.append(i)
      elif (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 < r3 ** 2:
        support.append(i)
    upper_bound = -353.03
    dy = 0.9
    for i in range(len(mesh.elements)):
      e = mesh.elements[i]
      f1 = False
      f2 = False
      f3 = False
      sp2 = False
      sp3 = False
      if e.region == "non-design":
        for j in range(len(e.nodes)):
          coord = mesh.nodes[e.nodes[j]]
          if abs(coord[1] - upper_bound) < dy and (coord[0] - m1[0]) ** 2 + (coord[2] - m1[2]) ** 2 < r1 ** 2:
            f1 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m2[0]) ** 2 + (coord[2] - m2[2]) ** 2 < r2 ** 2:
            f2 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m9[0]) ** 2 + (coord[2] - m9[2]) ** 2 < (r1+ delta) ** 2:
            f3 = True
          elif abs(coord[1] - upper_bound) < dy and (coord[0] - m8[0]) ** 2 + (coord[2] - m8[2]) ** 2 < (r3 + delta) ** 2:
            sp2 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m3[0]) ** 2 + (coord[2] - m3[2]) ** 2 < (r3 + delta) ** 2:
            sp3 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m4[0]) ** 2 + (coord[2] - m4[2]) ** 2 < (r3 + delta) ** 2:
            sp3 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m5[0]) ** 2 + (coord[2] - m5[2]) ** 2 < (r3 + delta) ** 2:
            sp3 = True
          elif abs(coord[1] - upper_bound) < dy  and (coord[0] - m6[0]) ** 2 + (coord[2] - m6[2]) ** 2 < (r3 + delta) ** 2:
            sp3 = True
          elif abs(coord[1] - upper_bound) < dy and (coord[0] - m7[0]) ** 2 + (coord[2] - m7[2]) ** 2 < (r3 + delta) ** 2:
            sp3 = True
        for j in range(len(e.nodes)):
          coord = mesh.nodes[e.nodes[j]]
          if f1 == True and abs(coord[1] - upper_bound) < dy:
              force1.append(e.nodes[j])
          elif f2 == True:
            force2.append(e.nodes[j])
          elif f3 == True:
            force3.append(e.nodes[j])
          elif sp2 == True:
            support2.append(e.nodes[j])    
          elif sp3 == True:
            support3.append(e.nodes[j])
  mesh.bc.append(('force1', force1))
  mesh.bc.append(('force2', force2))
  mesh.bc.append(('force3', force3))
  mesh.bc.append(('support', support))
  mesh.bc.append(('support2', support2))
  mesh.bc.append(('support3', support3))

  if len(elements) == 0:
    write_gid_mesh(mesh, meshfile+".mesh")     
  return mesh

def create_mesh_from_gmsh(meshfile,type):
  # read 3D tetrahedron gmsh mesh
  inp = open(meshfile+".msh").readlines()
  nodes = []
  elem = []
  if type == "apod6":
    force1 = []
    force2 = []
    force3 = []
    support = []
    support2 = []
    support3 = []
  count = 1
  num_node = 0
  num_elem = 0
  for line in inp:
    item = str.split(line)
    # read and check header
    if count == 2:
      if float(item[0]) != 2.2:
        print 'Error: Gmsh format should be 2.2, result probably wrong'
    # read number of nodes
    elif count == 5:
      num_node = int(item[0])
    #add nodes
    elif count > 5 and count <= num_node + 5:
      nodes.append([int(item[0]),float(item[1]),float(item[2]),float(item[3])])
    elif count > num_node + 5 and count <= num_node + 7:
      #skip lines
      count += 1 
      continue
    # read number of elements
    elif count == num_node + 8:
      num_elem = int(item[0])
    # add elements
    elif count > num_node + 8 and count <= num_node + 8 + num_elem:
      # read 3D tetrahedron elements
      if int(item[1]) == 4:
        elem.append([int(item[0]),int(item[5]),int(item[6]),int(item[7]),int(item[8])])
      # read 3D hexahedron elements
      elif int(item[1]) == 5:
        elem.append([int(item[0]),int(item[5]),int(item[6]),int(item[7]),int(item[8]),int(item[9]),int(item[10]),int(item[11]),int(item[12])])
      # read 3D wedge elements
      elif int(item[1]) == 6:
        elem.append([int(item[0]),int(item[5]),int(item[6]),int(item[7]),int(item[8]),int(item[9]),int(item[10])])
      elif int(item[1]) == 15 and type == "apod6":
        # force1
        if int(item[4]) == 5:
          force1.append(int(item[5])-1)
        # force2  
        elif int(item[4]) == 6:
          force2.append(int(item[5])-1)
        # force3  
        elif int(item[4]) == 9:
          force3.append(int(item[5])-1)
        # support  
        elif int(item[4]) == 7:
          support.append(int(item[5])-1)
        elif int(item[4]) == 8:
          support2.append(int(item[5])-1)
        elif int(item[4]) == 10:
          support3.append(int(item[5])-1)
    count += 1  
  nodes = numpy.asarray(nodes)
  if type == "apod6":
    mesh = create_mesh_for_apod6(meshfile,nodes,elem)
  elif type == "aux_cells":
    mesh = create_mesh_for_aux_cells(meshfile,nodes,elem)
  else:
    print "Error: No correct type was selected! options: apod6, aux_cells"
  write_gid_mesh(mesh, meshfile+".mesh") 
  
def create_gmsh_from_cfs_hdf5(hdf5_file, region, bcregions,output):
  # force names and support name has to be set manually, default force1, force2, force3, support, support2, support3
  mesh = create_mesh_from_hdf5(hdf5_file, region, bcregions)
  write_gid_mesh(mesh, "test.mesh") 
  out = open(output, "w")
  # gmsh header
  out.write('$MeshFormat \n')
  out.write('2.2 0 8\n')
  out.write('$EndMeshFormat \n')
  out.write('$Nodes \n')
  out.write(str(len(mesh.nodes))+' \n')
  dim = len(mesh.nodes[0])
  #write nodes
  for i in range(len(mesh.nodes)):  # write one based!
    out.write(str(i + 1) + "  " + str(mesh.nodes[i][0]) + "  " + str(mesh.nodes[i][1]))
    if dim == 3:
      out.write("  " + str(mesh.nodes[i][2]) + "\n")
    else:
      out.write("  0.0\n")
  #write elements
  out.write('$EndNodes \n')
  out.write('$Elements \n')
  out.write(str(len(mesh.elements)+len(mesh.bc[0][1]) + len(mesh.bc[1][1]) + len(mesh.bc[2][1]))+ '\n') #+ len(mesh.bc[3][1]) + len(mesh.bc[4][1]) + len(mesh.bc[5][1]) + len(mesh.bc[6][1]))+ '\n')
  # 1D boundary elements support, forces
  count = 0
  for k in range(len(mesh.bc)):
    bc = mesh.bc[k]
    if bc[0] == 'force1':
      id = 5
    elif bc[0] == 'force2':
      id = 6
    elif bc[0] == 'force3':
      id = 9
    elif bc[0] == 'support':
      id = 7
    elif bc[0] == 'support2':
      id = 8
    elif bc[0] == 'support3':
      id = 10
    else:
      print 'Warning mesh.bc type not handled!'
    for l in range(len(bc[1])):
      out.write(str(count+1) + ' ' +str(15) + ' 2 0 ' + str(id) + ' ' + str(bc[1][l] + 1)+' \n')
      count +=1
  # write 3D elements
  for i in range(len(mesh.elements)):  # write one based!
    e = mesh.elements[i]
    nodes = len(e.nodes)
    out.write(str(count + 1) + ' ' + str(5 if nodes_by_type(e.type) == 8 else 6) + ' 2 0 ' + str(2 if e.region == 'design' else 3))
    count +=1
    for n in range(nodes_by_type(e.type)):
      out.write(' '+ str(e.nodes[n] + 1))
    out.write('\n')
        
  out.write('$EndElements \n')
  #write forces, support
  
  out.write(' ')
  out.close()
  
def create_nastran_mesh_from_cfs(meshfile,h5file):
  # manually select cfs hdf5 file
  print 'Set regions and boundary nodes manually, default: design, non-design, force1,force2 and support'
  mesh = create_mesh_from_hdf5(h5file, ['design','non-design'], ['force1','force2','force3','support','support2','support3'])
  out = open(meshfile, "w")
  #design nodes and non-design nodes
  #out2 = open(meshfile + '.design', "w")
  #out3 = open(meshfile + 'non-desi"w")
  # nastran header
  out.write('ENDCONTROL\n')
  out.write('SUBCASE       1\n')
  out.write('  LABEL= SUBCASE 1\n')
  out.write('LOAD =       1\n')
  out.write('SUBCASE       2\n')
  out.write('  LABEL= SUBCASE 2\n')
  out.write('  LOAD =       2\n')
  out.write('BEGIN BULK\n')
  # write nodes
  for i in range(len(mesh.nodes)):
    n = mesh.nodes[i]
    out.write('GRID%12d%8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
    #out.write('GRID    ' + '%-8d%-8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
  # Hexaeder elements
  for i in range(len(mesh.elements)):
    e = mesh.elements[i]
    n = mesh.elements[i].nodes
    if e.type == HEXA8 and e.region == 'design':
      out.write('CHEXA%11d%8d%8d%8d%8d%8d%8d%8d+\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      out.write('+       %8d%8d\n'%(n[6]+1,n[7]+1))
      #out.write('CHEXA   ' + '%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+E%-6d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1,i+1))
      #out.write('+E%-6d%-8d%-8d\n'%(i+1,n[6]+1,n[7]+1))
    elif e.type == HEXA8 and e.region == 'non-design':
      out.write('CHEXA%11d%8d%8d%8d%8d%8d%8d%8d+\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      out.write('+       %8d%8d\n'%(n[6]+1,n[7]+1))
      #out.write('CHEXA   ' + '%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+E%-6d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1,i+1))
      #out.write('+E%-6d%-8d%-8d\n'%(i+1,n[6]+1,n[7]+1))
  # Wedge elements
  for i in range(len(mesh.elements)):
    e = mesh.elements[i]
    n = mesh.elements[i].nodes
    if e.type == WEDGE6 and e.region == 'design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == WEDGE6 and e.region == 'non-design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
  # write forces1
  for i in range(len(mesh.bc[0][1])):
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(1,mesh.bc[0][1][i]+1,0,5000./len(mesh.bc[0][1])))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(1,mesh.bc[0][1][i]+1,0,5000./len(mesh.bc[0][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

    # write forces2
  for i in range(len(mesh.bc[1][1])):
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,mesh.bc[1][1][i]+1,0,5000./len(mesh.bc[1][1])))
    
      # write forces3
  for i in range(len(mesh.bc[2][1])):
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,mesh.bc[2][1][i]+1,0,5000./len(mesh.bc[2][1])))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(2,mesh.bc[1][1][i]+1,0,5000./len(mesh.bc[1][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

  for i in range(len(mesh.bc[3][1])):
    out.write('SPC%13d%8d  13     \n'%(1,mesh.bc[3][1][i]+1))
  for i in range(len(mesh.bc[4][1])):
    out.write('SPC%13d%8d  2     \n'%(1,mesh.bc[4][1][i]+1))
  for i in range(len(mesh.bc[5][1])):
    out.write('SPC%13d%8d  2     \n'%(1,mesh.bc[5][1][i]+1))
    #out.write('SPC     ' + '%-8d%-8d%-8d%-8d%-8d%-8f\n'%(1,mesh.bc[2][1][i]+1,1,2,3,0.))
    
  #out.write('PSOLID         1       1\n')          
  #out.write('PSOLID         2       1\n')
  out.write('PSOLID  1       1       \n')          
  out.write('PSOLID  2       1       \n')            
  #out.write('MAT1    1       1.00E0  0.34     0.0     0.785E-5  12.E-6                +M1\n')     
  #out.write('+M1     100.    -100.   100.\n')    
  #out.write('MAT1    2       7.00E4  0.34     0.0     0.785E-5  12.E-6                +M1\n')     
  #out.write('+M1     100.    -100.   100.\n')
  # Ti6Al4
  out.write('MAT1    %-8d%-.2e        %-8f%-8f%-8f%-8f%-8f\n'%(1,1.2E11,0.342,1.0,1.0,0.,1.))
  # Aluminum
  out.write('MAT1    %-8d%-.2e        %-8f%-8f%-8f%-8f%-8f\n'%(2,6.8E10,0.36,1.0,1.0,0.,1.)) 
  out.write('ENDDATA\n')
  
  out.close()  
  
def create_mesh_apod6_from_optistruct(meshfile,scale):
  # read 3D optistruct mesh with hexa and wedge elements for apod6 got by M. Muir (12/2015)
  inp = open(meshfile).readlines()
  nodes = []
  elem = []
  force1 = []
  force2 = []
  force3 = []
  support = []
  support2 = []
  support3 = []
  design = []
  nondesign = []
  count = 1
  num_node = 0
  num_elem = 0
  des = False
  nondes = False
  for i in range(len(inp)):
    item = str.split(inp[i])
    if i < len(inp)-1:
      item_n = str.split(inp[i+1])
    if len(item) == 0:
      continue
    # read and check header
    if item[0] == 'GRID':
      #add nodes
      if len(item) > 3:
        x = float(item[2][:])
        y = float(item[3][0:8])
        z = float(item[3][8:16])
      else:
        x = float(item[2][0:8])
        y = float(item[2][8:16])
        z = float(item[2][16:24])
      nodes.append([int(item[1]),x,y,z])
    elif item[0] == 'CPENTA':
      # read 3D wedge elements
      elem.append([int(item[1]),int(item[3]),int(item[4]),int(item[5]),int(item[6]),int(item[7]),int(item[8])])
    elif item[0] == 'CHEXA':
      # read 3D hexahedron elements
      elem.append([int(item[1]),int(item[3]),int(item[4]),int(item[5]),int(item[6]),int(item[7]),int(item[8]),int(item_n[1]),int(item_n[2])])
    elif item[0] == '$HMMOVE' and item[1] == '6':
      # set flag for design material
      des = True
    elif item[0] == '$HMMOVE' and item[1] == '7':
      # set flag for nondesign material
      des = False
      nondes = True
    elif des == True and len(item) > 1:
      # nondesign elements
      for i in range(len(item)-1):
        design.append(item[i+1])
    elif nondes == True and len(item) > 1:
      # add nondesign elements
      for i in range(len(item)-1):
        nondesign.append(item[i+1])
    else:
      des = False
      nondes = False
    count += 1  
  nodes = numpy.asarray(nodes)
  mesh = create_mesh_for_apod6(meshfile,nodes,elem)
  write_gid_mesh(mesh, meshfile+".mesh",scale) 
  
def create_mesh_for_aux_cells(meshfile, all_nodes = [], elements = []):
  mesh = Mesh()
  
  # insert nodes  
  for i in range(len(all_nodes)):
    mesh.nodes.append([all_nodes[i, 1], all_nodes[i, 2], all_nodes[i, 3]])
  min_diam_x = 1000000. 
  min_diam_y = 1000000. 
  min_diam_z = 1000000. 
  # insert elements
  for i in range(len(elements)):
      e = Element()
      e.nodes = (elements[i][1:])
      for k in range(len(e.nodes)):
        e.nodes[k] -= 1
      count = 0
      for k in range (len(e.nodes)):
        if count + 1 == len(e.nodes):
          min_diam_x = min(min_diam_x,abs(mesh.nodes[e.nodes[count]][0] - mesh.nodes[e.nodes[0]][0]))
          min_diam_y = min(min_diam_y,abs(mesh.nodes[e.nodes[count]][1] - mesh.nodes[e.nodes[0]][1]))
          min_diam_z = min(min_diam_z,abs(mesh.nodes[e.nodes[count]][2] - mesh.nodes[e.nodes[0]][2]))
        else:
          min_diam_x = min(min_diam_x,abs(mesh.nodes[e.nodes[count]][0] - mesh.nodes[e.nodes[count+1]][0]))
          min_diam_y = min(min_diam_y,abs(mesh.nodes[e.nodes[count]][1] - mesh.nodes[e.nodes[count+1]][1]))
          min_diam_z = min(min_diam_z,abs(mesh.nodes[e.nodes[count]][2] - mesh.nodes[e.nodes[count+1]][2]))
        count += 1
      print "min_diam = [" + str(min_diam_x) + ", " +  str(min_diam_y) + ", " +  str(min_diam_z)
      e.density = 1.
      e.region = 'mech'
      if len(e.nodes) == 4:
        e.type = TET4
      elif len(e.nodes) == 6:
        e.type = WEDGE6
      elif len(e.nodes) == 8:
        e.type = HEXA8
      mesh.elements.append(e)
       
  # define boundary nodes
  print "min_diam = [" + str(min_diam_x) + ", " +  str(min_diam_y) + ", " +  str(min_diam_z)
  # calculate extrem values along each coordinate
  ma_x = -100000.
  mi_x = 100000.
  ma_y = -100000.
  mi_y = 100000.
  ma_z = -100000.
  mi_z = 100000.
  for i in range(len(mesh.nodes)):
    coord = mesh.nodes[i]
    ma_x = ma_x if ma_x > coord[0] else coord[0]
    mi_x = mi_x if mi_x < coord[0] else coord[0]
    ma_y = ma_y if ma_y > coord[1] else coord[1]
    mi_y = mi_y if mi_y < coord[1] else coord[1]
    ma_z = ma_z if ma_z > coord[2] else coord[2]   
    mi_z = mi_z if mi_z < coord[2] else coord[2]

  delta = 1e-4
  top = []
  bottom = []
  left = []
  right = []
  front = []
  back = []
  for i in range(len(mesh.nodes)):
    if abs(mesh.nodes[i][0] - mi_x) < min_diam_x - delta:
      left.append(i)
    elif abs(mesh.nodes[i][0] - ma_x) < min_diam_x - delta:
      right.append(i)
    elif abs(mesh.nodes[i][1] - mi_y) < min_diam_y - delta:
      bottom.append(i) 
    elif abs(mesh.nodes[i][1] - ma_y) < min_diam_y - delta:
      top.append(i) 
    elif abs(mesh.nodes[i][2] - mi_z) < min_diam_z - delta:
      back.append(i)
    elif abs(mesh.nodes[i][2] - ma_z) < min_diam_z - delta:
      front.append(i)
      
  #add boundary nodes    
  mesh.bc.append(('top', top))
  mesh.bc.append(('bottom', bottom))
  mesh.bc.append(('left', left))
  mesh.bc.append(('right', right))
  mesh.bc.append(('front', front))
  mesh.bc.append(('back', back))
  
  return mesh
  
