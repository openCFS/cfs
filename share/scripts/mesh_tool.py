#!/usr/bin/env python

import Image, sys, os, copy, numpy, math
from hdf5_tools import *


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
  def __init__(self):
   self.nodes = []  # list 2d tupels (float, float) or 3d tuples
   self.elements = []  # list of Element
   # list of boundary conditon nodes
   self.bc = []  # list of tupel (name, <list of zero based nodes>)
   # list of named nodes (save element in gd)
   self.ne = []  # list of tupel (name, <list of zero based elements>)

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


def create_dense_mesh_img(input_img, mesh, threshold, scale, rhomin, shearAngle, pressure=False, color_mode="random"):
  input_pix = input_img.load()
  nx, ny = input_img.size
  create_dense_mesh(input_pix, nx, ny, mesh, threshold, scale, rhomin, 1, shearAngle, pressure, color_mode)

def create_dense_mesh_density(numpy_array, mesh, threshold, scale, rhomin, multi_d=1):
  if multi_d == 1:
    nx, ny = numpy_array.shape
  else:
    nx, ny, nz, m = numpy_array.shape
  create_dense_mesh(numpy_array, nx, ny, mesh, threshold, scale, rhomin, multi_design=multi_d, shearAngle=0.0)
  
def create_dense_mesh(input_array, nx, ny, mesh, threshold, scale, rhomin, multi_design=1, shearAngle=0, pressure=False, color_mode="random"):
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
        if is_gray or (input_array[x, y][0] == input_array[x, y][1] and input_array[x, y][1] == input_array[x, y][2]):  
          if e.density >= threshold:
            e.region = 'mech'
            mech_count += 1
          else:
            e.region = 'void'
        else:
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
  if pressure:
    y = ny - 1
    for x in range(nx):
      if (x >= int(0.8 * nx) and y == ny - 1):
        b = Element()
        b.type = LINE
        ll = x
        b.nodes = ((ll, ll + 1))
        if e.region == 'mech' or e.region == 'colorful' or e.region == 'red': 
          b.region = 'pressure2'
        else:
          b.region = 'void'
        mesh.elements.append(b)
    print 'Warning: pressure area has to be set manually in method create_dense_mesh.'
  mesh.bc.append(("south", range(0, nx + 1)))
  mesh.bc.append(("north", range((nx + 1) * ny, (nx + 1) * (ny + 1))))
  mesh.bc.append(("west", range(0, (nx + 1) * ny + 1, nx + 1)))
  mesh.bc.append(("east", range(nx, (nx + 1) * (ny + 1), nx + 1)))
  mesh.bc.append(("pressure2", range(int(0.8 * nx), nx + 1)))



  
  # print mesh.bc

  print "dense resolution: " + str(nx) + " x " + str(ny) + " elements (" + str(scale) + "m x " + str(float(ny) / nx * scale) + "m)",
  print " -> " + str(mech_count) + " mech elements out of " + str(nx * ny) + " (" + str(float(mech_count) / (nx * ny) * 100.0) + " %)",
  print " with threshold " + str(threshold) 
  if colorful_count > 0:
    print 'plus ' + str(colorful_count) + ' non gray elements'

# @param mesh dense mesh (input)
# @return sparse mesh
def convert_to_sparse_mesh(dense):
  sparse = Mesh()

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
      out.write(str(i + 1) + ' ' + str(e.type) + ' ' + str(nodes) + ' ' + e.region + "\n")
    
      # prepare for second order elements
      for n in range(nodes):
        out.write(str(e.nodes[n] + 1) + ("\n" if n == len(e.nodes) - 1 else " "))  # write one based node numbers


def write_gid_mesh(mesh, filename):
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
    out.write(str(i + 1) + "  " + str(mesh.nodes[i][0]) + "  " + str(mesh.nodes[i][1]))
    if dim == 3:
      out.write("  " + str(mesh.nodes[i][2]) + "\n")
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
  out.write('#ElemNr Level\n\n')
  for e in range(len(mesh.ne)):
    ne = mesh.ne[e]
    for n in range(len(ne[1])):
      out.write(str(ne[1][n] + 1) + " " + ne[0] + "\n")

  out.close()
  

# # creates a 2D mesh of predefined geometry
def create_2d_mesh(type, x_res, y_res):
  mesh = Mesh()
  
  # mbb reinforced does not work!!! check and fix the mbb stuff!
  nx = x_res
  # buld2d case
  ny = y_res if y_res <> None else x_res
  width = 1.0
  height = float(ny) / nx 
  # offset for node coordinates
  offx = 0.
  offy = 0.
  if type.startswith('cantilever2d'):
    width = 3.0
    height = 2.0
    ny = int(nx * (2. / 3.))
  if type.startswith('mbb'):
    width = 2.0
    height = 1.0
    ny = int(nx * 0.5)
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
     
  dx = width / nx
  dy = height / ny

  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((offx + x * dx, offy + y * dy))
  
  scale = 2 if type == 'triangle_msfem' else 1
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
      else:
        e.region = 'mech'
      # assign nodes
      ll = (nx + 1) * y + x  # lowerleft
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
  # if type == 'msfem_test':
  #  for y in range(ny):
  #    b = Element()
  #    b.type = LINE
  #    ll = (nx + 1) * y + nx - 1
  #    b.nodes = ((ll + 1, ll + 1 + nx + 1)) 
  #    b.region = 'pressure'
  #   mesh.elements.append(b)
  if type == 'msfem_test2':
    for x in range(int(0.8 * nx), nx):
      b = Element()
      b.type = LINE
      ll = x
      b.nodes = ((ll, ll + 1)) 
      b.region = 'pressure2'
      mesh.elements.append(b)
      
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

# # creates a mesh of predefined geometry
def create_3d_mesh(x_res, y_res, z_res):
  mesh = Mesh()

  nx = x_res
  ny = y_res if y_res <> None else x_res
  nz = z_res if z_res <> None else x_res

  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1
  
  width = 1.0
  dx = width / nx
  
  height = float(ny) / nx
  dy = height / ny
  
  depth = float(nz) / nx
  dz = depth / nz

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
 
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.density = 1.0
        e.type = HEXA8
        e.region = 'mech'
  
        # assign nodes
        # ll = (nx+1)*y*(nx+1) * z + (nx+1) * y + x  # lowerleftfront
        ll = nnx * nny * z + nnx * y + x  # lower-left-front of current element
        # start with upper-front-left counterclockwise in the x-z plane. Repeat in then lower plane
        # e.nodes = ((ll+(nx+1), ll+1+(nx+1), ll+1+(nx+1)+((nx+1)*(ny+1)),ll+(nx+1)+((nx+1)*(ny+1)),ll, ll+1, ll+1+((nx+1)*(ny+1)),ll+((nx+1)*(ny+1))))  
        e.nodes = ((ll + nnx, ll + 1 + nnx, ll + 1 + nnx + (nnx * nny), ll + nnx + (nnx * nny), ll, ll + 1, ll + 1 + (nnx * nny), ll + (nnx * nny)))
        
        mesh.elements.append(e)

  mesh.bc.append(("left", range(0, (nnx * nny * z) + (nnx * ny) + 1, nnx)))
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
  
  return mesh

# # LBM pipe_bend and two_inlet_one_outlet example as used by Pingen et al. 2007
# @param case pipe_bend or two_inlet_one_outlet 
def create_lbm(resolution, case):
  mesh = Mesh()
 
  size = 1.0 
   
  nx = resolution
  ny = nx
  
  dx = size / nx
  
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
      if x > 0 and y > 0 and x < nx - 1 and y < nx - 1: 
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
      
  return mesh

# creates a mesh from hdf5 file  
def create_mesh_from_hdf5(hdf5_file, region, bcregions, region_force=None, region_support=None, threshold=0.):
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'].value  # for all regions
  reg_elements_des = hdf5_file['/Mesh/Regions/' + region[0] + '/Elements'].value
  if len(region) > 1:
    reg_elements_nondes = hdf5_file['/Mesh/Regions/' + region[1] + '/Elements'].value
    reg_elements_void = hdf5_file['/Mesh/Regions/' + region[2] + '/Elements'].value
  types = hdf5_file['/Mesh/Elements/Types'].value
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates'].value
  length = len(hdf5_file['/Mesh/Regions/' + region[0] + '/Nodes'].value)
  reg_nodes = [[0 for col in range(len(region))] for row in range(length)]
  for i in range(len(region)):
    reg_nodes[i][:] = hdf5_file['/Mesh/Regions/' + region[i] + '/Nodes']
  design_var = hdf5_file['/Results/Mesh/MultiStep_1/Step_0/physicalPseudoDensity/mech/Elements/Real'].value
    
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
    for i in range(len(bcregions)):
      bc_nodes = hdf5_file['Mesh/Groups/' + str(bcregions[i]) + '/Nodes']
      mesh.bc.append((bcregions[i], bc_nodes[:] - 1))
  
  for i in range(len(all_nodes)):
    mesh.nodes.append(all_nodes[i])
  idx = 0
  idx2 = 0
  idx3 = 0  
  for i in range(len(all_elements[:, 0])):
    e = Element()
    e.nodes = (all_elements[i, :] - 1)
    e.density = design_var[i]
    if idx < len(reg_elements_des):
      if i + 1 == reg_elements_des[idx]:
        if e.density >= threshold:
          e.region = region[0]
        else:
          e.region = 'void'
        idx += 1
    if len(region) > 1:
      if idx2 < len(reg_elements_nondes):
        if i + 1 == reg_elements_nondes[idx2]:
          e.region = 'nondesign'
          idx2 += 1
      if idx3 < len(reg_elements_void):
        if i + 1 == reg_elements_void[idx3]:
          e.region = 'void'
          idx3 += 1
    e.type = mesh_type_from_hdf5(types[i])
    mesh.elements.append(e) 
  return mesh


def create_mesh_from_tetgen(meshfile, region):
  print meshfile + '.1.ele'
  all_elements = numpy.loadtxt(meshfile + '.1.ele', dtype='int' , skiprows=1)
  print 'read all_elements done'
  all_nodes = numpy.loadtxt(meshfile + '.1.node', skiprows=1)
  print 'read all_nodes done'
  # all_faces = numpy.loadtxt(meshfile+'1.face',skiprows=1)
  # all_edges = numpy.loadtxt(meshfile+'1.edge',skiprows=1)
  
    
  # Create mesh  
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

def create_mesh_from_optistruct(meshfile):
  # create element and nodes files by hand from optistruct
  
  # load files
  hexa_elements = numpy.loadtxt(meshfile + '.hexa.elements', dtype='int', skiprows=1)
  wedge_elements = numpy.loadtxt(meshfile + '.wedge.elements', dtype='int', skiprows=1)
  all_nodes = numpy.loadtxt(meshfile + '.nodes', skiprows=1)
  # ax = 0.0849
  ay = -0.084636333418591
  # az = 0.1358
  # Rx = numpy.matrix(((1., 0., 0.), (0., math.cos(ax), -math.sin(ax)), (0., math.sin(ax), math.cos(ax))))
  Ry = numpy.matrix(((math.cos(ay), 0., math.sin(ay)), (0., 1., 0.), (-math.sin(ay), 0., math.cos(ay))))
  # Rz = numpy.matrix(((math.cos(az), -math.sin(az), 0.), (math.sin(az), math.cos(az), 0.), (0., 0., 1.)))
  # Create mesh  
  # add nodes    
  mesh = Mesh()  
  for i in range(len(all_nodes)):
    coord = numpy.matrix(((all_nodes[i, 1]), (all_nodes[i, 2]), all_nodes[i, 3])).T
    # print Rx
    # print coord  
    new_coord = Ry * coord
    # new_coord = Rz * new_coord
    mesh.nodes.append([new_coord[0, 0], new_coord[1, 0], new_coord[2, 0]])
    
  # hexaeder     
  for i in range(len(hexa_elements[:, 0])):
    e = Element()
    e.nodes = (hexa_elements[i, 1:] - 1)
    e.density = 1.
    shell = 0
    for j in range(8):
      coord = mesh.nodes[e.nodes[j]]
      if (coord[1] + 353.) < 1. or abs(coord[1] + 333.) < 1.:
        shell += 1
    if shell > 3:
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
      if (coord[1] + 353.) < 1. or abs(coord[1] + 333.) < 1.:
        shell += 1
    if shell > 2:
      e.region = 'non-design'
    else:
      e.region = 'design'
    e.type = WEDGE6
    mesh.elements.append(e)
  
  # include boundary conditions for 442.mesh manually
  m1 = [33052., -353., -2474.]
  m2 = [33046., -353., -2518.]
  m3 = [33131., -353., -2449.]
  m4 = [33124., -353., -2498.]
  m5 = [32978., -353., -2436.]
  m6 = [32971., -353., -2485.]
  m7 = [33023., -353., -2559.]
  r1 = 19.5
  r2 = 16.5
  r3 = 5.8
  force1 = []
  force2 = []
  support = []
  for i in range(len(all_nodes)):
    coord = all_nodes[i, 1:]
    if abs(coord[1] + 353.) < 1. and (coord[0] - m1[0]) ** 2 + (coord[2] - m1[2]) ** 2 < r1 ** 2:
      force1.append(i)
    elif abs(coord[1] + 353.) < 1.  and (coord[0] - m2[0]) ** 2 + (coord[2] - m2[2]) ** 2 < r2 ** 2:
      force2.append(i)
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
  mesh.bc.append(('support', support))    
  return mesh
  
