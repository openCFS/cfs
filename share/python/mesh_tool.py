#!/usr/bin/env python
import platform
# from PIL import Image # load only in the function we need it!
import sys, os, numpy, math
from copy import deepcopy
try:
  from hdf5_tools import *
except:
  print("failed to import hdf5_tools, hopefully we don't need it")    
import scipy.interpolate as ip
from numpy import ceil
import scipy.spatial
from special_mesh_tools import *


# writes a dense two region mesh
# def write_dense_mesh(pixels, size, file, threshold):

# element types as in gid (simInputMESH.cc -> AnsysType2ElemType)
TRIA3 = 4
QUAD4 = 6
TET4 = 8
HEXA8 = 10
WEDGE6 = 14
LINE = 100

def set_array_point(array,point,hx,hy,hz,minx,miny,minz,val):
  i = int((point[0] - minx)/hx)
  j = int((point[1] - miny)/hy)
  k = int((point[2] - minz)/hz)
  array[i,j,k] = val

# checks whether point p is inside cube
# @param tetra: list with 4 vertices of tetrahedron
def inside_cube(p,tetra):
  d0 = numpy.linalg.det(numpy.array([[tetra[0][0],tetra[0][1],tetra[0][2],1], [tetra[1][0],tetra[1][1],tetra[1][2],1], [tetra[2][0],tetra[2][1],tetra[2][2],1], [tetra[3][0],tetra[3][1],tetra[3][2],1]]))
  d1 = numpy.linalg.det(numpy.array([[p[0],p[1],p[2],1], [tetra[1][0],tetra[1][1],tetra[1][2],1], [tetra[2][0],tetra[2][1],tetra[2][2],1], [tetra[3][0],tetra[3][1],tetra[3][2],1]]))
  d2 = numpy.linalg.det(numpy.array([[tetra[0][0],tetra[0][1],tetra[0][2],1], [p[0],p[1],p[2],1], [tetra[2][0],tetra[2][1],tetra[2][2],1], [tetra[3][0],tetra[3][1],tetra[3][2],1]]))
  d3 = numpy.linalg.det(numpy.array([[tetra[0][0],tetra[0][1],tetra[0][2],1], [tetra[1][0],tetra[1][1],tetra[1][2],1], [p[0],p[1],p[2],1], [tetra[3][0],tetra[3][1],tetra[3][2],1]]))
  d4 = numpy.linalg.det(numpy.array([[tetra[0][0],tetra[0][1],tetra[0][2],1], [tetra[1][0],tetra[1][1],tetra[1][2],1], [tetra[2][0],tetra[2][1],tetra[2][2],1], [p[0],p[1],p[2],1]]))
  
  det = [d0,d1,d2,d3,d4]
  # point is inside cube if all determinats have same sign
  return all(item >= 0 for item in det) or all(item < 0 for item in det)
# calculates barycenter of element e
def calc_barycenter(mesh,e):
  center = numpy.array([0.0, 0.0, 0.0])
  len_nod = len(e.nodes)
  for n in range(len_nod):
    center += mesh.nodes[e.nodes[n]]
  center *= 1.0 / len_nod
  
  return center

# node1/2: coordinates of respective node
def calc_edge_length(mesh,node1,node2):
  return numpy.linalg.norm(numpy.array(node1)-numpy.array(node2))

# calculates longest edge of element e
def calc_longest_edge(mesh,e):
  # so far only tested for 3D elements
  assert(elem_dim(e.type) == 3) 
  result = 0
  
  nodes = e.nodes
  
  for i,node in enumerate(nodes):
    for j,other in enumerate(nodes):
      if i == j:
        continue
      # euklidian distance
      distance = calc_edge_length(mesh, mesh.nodes[node], mesh.nodes[other])
    
      if distance > result:
        result = distance  
    
  return result  

def calc_min_max_coords(mesh):
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
    

  return mi_x, mi_y, mi_z, ma_x, ma_y, ma_z

# for given node ide (e.g. from optistruct), return id in gid mesh
#def node_id_to_mesh_id(nodeId):
  

def nodes_by_type(type):
  if type == QUAD4:
    return 4
  if type == HEXA8:
    return 8
  if type == WEDGE6:
    return 6
  if type == TRIA3:
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
    return TRIA3
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
    print(self.nodes)
    print(self.region)
    print(self.density)     


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
  
# adds same number of boundary nodes on adjacent sides to assure periodic b.c
# @ min_diam
def add_nodes_for_periodic_bc(mesh,min_diam_x=1e-3,min_diam_y=1e-3,min_diam_z=1e-3,delta=1e-3):
  left_c = 0
  right_c = 0
  top_c = 0
  bottom_c = 0
  back_c = 0
  front_c = 0
  
  top = []
  bottom = []
  left = []
  right = []
  front = []
  back = []
  
  mi_x, mi_y, mi_z, ma_x, ma_y, ma_z = calc_min_max_coords(mesh)
  
  # count number of boundary nodes per region
  for i in range(len(mesh.nodes)):
    if abs(mesh.nodes[i][0] - mi_x) < min_diam_x + delta:
      left_c += 1
    elif abs(mesh.nodes[i][0] - ma_x) < min_diam_x + delta:
      right_c += 1
    elif abs(mesh.nodes[i][1] - mi_y) < min_diam_y + delta:
      bottom_c += 1 
    elif abs(mesh.nodes[i][1] - ma_y) < min_diam_y + delta:
      top_c += 1
    elif abs(mesh.nodes[i][2] - mi_z) < min_diam_z + delta:
      back_c += 1
    elif abs(mesh.nodes[i][2] - ma_z) < min_diam_z + delta:
      front_c += 1
  
  lr_counter = min(left_c,right_c)
  bt_counter = min(bottom_c,top_c)
  bf_counter = min(back_c,front_c)
  
  print("left_c: ", left_c, " right_c: ", right_c, " bottom_c:", bottom_c, " top_c: ", top_c, " back_c: ", back_c, " front_c: ",front_c) 
  
  left_c = 0
  right_c = 0
  top_c = 0
  bottom_c = 0
  back_c = 0
  front_c = 0
  
  for i in range(len(mesh.nodes)):
    if abs(mesh.nodes[i][0] - mi_x) < min_diam_x + delta and left_c < lr_counter:
      left.append(i)
      left_c +=1
    elif abs(mesh.nodes[i][0] - ma_x) < min_diam_x + delta and right_c < lr_counter:
      right.append(i)
      right_c +=1
    elif abs(mesh.nodes[i][1] - mi_y) < min_diam_y + delta and bottom_c < bt_counter:
      bottom.append(i)
      bottom_c +=1 
    elif abs(mesh.nodes[i][1] - ma_y) < min_diam_y + delta and top_c < bt_counter:
      top.append(i)
      top_c +=1 
    elif abs(mesh.nodes[i][2] - mi_z) < min_diam_z + delta and back_c < bf_counter:
      back.append(i)
      back_c +=1
    elif abs(mesh.nodes[i][2] - ma_z) < min_diam_z + delta and front_c < bf_counter:
      front.append(i)
      front_c +=1 
      
  print("left_c: ", left_c, " right_c: ", right_c, " bottom_c:", bottom_c, " top_c: ", top_c, " back_c: ", back_c, " front_c: ",front_c)    
  
  mesh.bc = []
  #add boundary nodes    
  mesh.bc.append(('top', top))
  mesh.bc.append(('bottom', bottom))
  mesh.bc.append(('left', left))
  mesh.bc.append(('right', right))
  mesh.bc.append(('front', front))
  mesh.bc.append(('back', back))
  
  return mesh     
  
def show_dense_mesh_image(mesh, shape, binary, size):
  from PIL import Image  
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
      
  check_img = check_img.resize((size, int(ny * size / nx)))    
  check_img.show()     


def create_dense_mesh_img(input_img, mesh, threshold, scale, rhomin, shearAngle, type = 1, color_mode="random"):
  input_pix = input_img.load()
  nx, ny = input_img.size
  create_dense_mesh(input_pix, nx, ny, mesh, threshold, scale, rhomin, 1, shearAngle, type, color_mode)

def create_dense_mesh_density(numpy_array, mesh, threshold, scale, rhomin, multi_d=1):
  # handle different types of numpy_array, further description in create_dense_mesh
  # only one design variable
  if multi_d == 1 and len(numpy_array.shape) == 3: 
    nx, ny, nz = numpy_array.shape   
    create_3d_mesh('bulk3d', x_res = nx, y_res = ny, z_res = nz, data = numpy_array, threshold = threshold, ext_mesh = mesh, scale = scale)
  else:
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
    print('angle has to be between -pi/2 + 1e-6 and pi/2 - 1e-6')
    return 0 
  dx = scale / nx / math.cos(angle)
  # from daniel ?! dy = scale/ny
  dy = dx
  
  mesh.nx = nx
  mesh.ny = ny
  
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
  mesh.bc.append(("south", list(range(0, nx + 1))))
  mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
  mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
  mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
  if type == 2:
    # pressure boundary for a test case
    print('Warning: pressure area has to be set manually in method create_dense_mesh.')
    mesh.bc.append(("pressure2", list(range(int(0.8 * nx), nx + 1))))
  elif type == 3:
    # other boundary conditions for msfem test
    # lower/upper loads
    mid = int((nx+1.)/2.)
    off_x = int(0.05 * nx)
    off_y = int(0.1 * ny)
    if (nx+1)/2. % 2 == 0:
      mesh.bc.append(("load1", list(range(mid-off_x,mid + off_x))))
      mesh.bc.append(("load2", list(range((nx+1)*ny + mid-off_x, (nx+1)*ny + mid + off_x))))
    else:
      mesh.bc.append(("load1", list(range(mid-off_x,mid + 1 + off_x))))
      mesh.bc.append(("load2", list(range((nx+1)*ny + mid-off_x,(nx+1)*ny + mid + off_x + 1))))
    # support lower left
    mesh.bc.append(("support", list(range(0,off_x+1))))
    mesh.bc.append(("support", list(range(0,(nx+1)*off_y+1,nx+1))))
    # support lower right
    mesh.bc.append(("support", list(range(nx - off_x,nx+1))))
    mesh.bc.append(("support", list(range(nx,(nx+1)*(off_y+1),nx+1))))
    # support upper left
    mesh.bc.append(("support", list(range((nx+1)*ny, (nx+1)*ny+off_x+1))))
    mesh.bc.append(("support", list(range((nx+1)*ny-(nx+1)*off_y,(nx+1)*ny+1,nx+1))))

    # support upper right
    mesh.bc.append(("support", list(range((nx+1)*(ny+1) -1 - off_x,(nx+1)*(ny+1)))))
    mesh.bc.append(("support", list(range((nx+1)*(ny+1)-(nx+1)*off_y-1,(nx+1)*(ny+1),nx+1))))


  
  # print mesh.bc
  msg = "dense resolution: " + str(nx) + " x " + str(ny) + " elements (" + str(scale) + "m x " + str(float(ny) / nx * scale) + "m)"
  msg += " -> " + str(mech_count) + " mech elements out of " + str(nx * ny) + " (" + str(float(mech_count) / (nx * ny) * 100.0) + " %)"
  msg += " with threshold " + str(threshold)
  if colorful_count > 0:
    msg += 'plus ' + str(colorful_count) + ' non gray elements'
  print(msg)  
    
# @param mesh dense mesh (input)
# @return sparse mesh
def convert_to_sparse_mesh(dense):
  sparse = Mesh(dense.nx, dense.ny, dense.nz)

  # necessary 0-based nodes as unique set
  nns = set()
  
  # copy element, the indices of the nodes will be replaced later
  for i in range(len(dense.elements)):
    e = dense.elements[i]
    if e.region != 'void':
      sparse.elements.append(deepcopy(e))
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
       assert(el.nodes[n] != -1)
     el.nodes = newnodes

  # finally handle the boundary conditions
  sparse.bc = []
  for b in range(len(dense.bc)):
    bc = dense.bc[b]
    dnn = bc[1]  # dense nodes
    nodes = []
    for n in range(len(dnn)):
      print('old number '+str(dnn[n]) + ' new number '+str(map[dnn[n]]))
      if map[dnn[n]] != -1:
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
  quad4 = count_elements(mesh.elements, QUAD4)
  hexa8 = count_elements(mesh.elements, HEXA8)
  wedge6 = count_elements(mesh.elements, WEDGE6)
  line = count_elements(mesh.elements, LINE)
  tet4 = count_elements(mesh.elements, TET4)
  tri3 = count_elements(mesh.elements, TRIA3)
  num_1d = line
  num_2d = quad4 + tri3
  num_3d = hexa8 + wedge6 + tet4
  assert(num_1d + num_2d + num_3d == len(mesh.elements))
  print('number of elements ' + str(num_1d + num_2d + num_3d))
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
      out.write(str(bc[1][n] + 1) + " " + bc[0] + "\n") # bc[1][n]+1 is node number and bc[0] label
  
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

# names all nodes on edges, faces and corners of cubic domain
def name_bc_nodes(mesh):
  assert(mesh != None)
  assert(mesh.nz > 1)
  nx = mesh.nx
  ny = mesh.ny
  nz = mesh.nz
  
  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1
  
  # naming faces of cube
  mesh.bc.append(("left", list(range(0, (nnx * nny * nz) + (nnx * ny) + 1, nnx))))
  mesh.bc.append(("right", list(range(nx, (nnx * nny * nnz) + 1, nnx))))

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
  mesh.bc.append(("back", list(range(0, (nx + 1) * (ny + 1)))))
  mesh.bc.append(("front", list(range(nz * (nx + 1) * (ny + 1), (nz + 1) * (nx + 1) * (ny + 1)))))

  # naming cube corners
  mesh.bc.append(("left_bottom_back", [0]))
  mesh.bc.append(("right_bottom_back", [nx]))
  mesh.bc.append(("left_top_back", [nnx * ny]))
  mesh.bc.append(("right_top_back", [nnx * nny - 1]))
  mesh.bc.append(("left_bottom_front", [nnx * nny * nz]))
  mesh.bc.append(("right_bottom_front", [nnx * nny * nz + nx]))
  mesh.bc.append(("left_top_front", [nnx * nny * nz + nnx * ny]))
  mesh.bc.append(("right_top_front", [nnx * nny * nnz - 1]))
  
  # naming cube edges
  mesh.bc.append(("bottom_back",list(range(nnx))))
  mesh.bc.append(("bottom_front",list(range(nnx*nny*(nnz-1),nnx*nny*(nnz-1)+nnx))))
  mesh.bc.append(("bottom_left",list(range(0,nnx*nny*nnz-nnx-1,nnx*nny))))
  mesh.bc.append(("bottom_right",list(range(nnx-1,nnx*nny*nnz-1-1,nnx*nny))))
  mesh.bc.append(("top_back",list(range(nnx*nny-nnx,nnx*nny))))
  mesh.bc.append(("top_front",list(range(nnx*nny*nny-nnx,nnx*nny*nny))))
  mesh.bc.append(("top_left",list(range(nnx*nny-nnx,nnx*nny*nnz,nnx*nny))))
  mesh.bc.append(("top_right",list(range(nnx*nny-1,nnx*nny*nnz,nnx*nny))))
  mesh.bc.append(("back_left",list(range(0,nnx*nny-nnx+1,nnx))))
  mesh.bc.append(("back_right",list(range(nnx-1,nnx*nny,nnx))))
  mesh.bc.append(("front_left",list(range(nnx * nny * nz,nnx * nny * nz + nnx * ny+1,nnx))))
  mesh.bc.append(("front_right",list(range(nnx * nny * nz + nx,nnx * nny * nnz,nnx))))
  
  
  return mesh  
 
def validate_periodicity(mesh):
#   assert(mesh.nz > 1)
  countLeft = len([x for x in mesh.bc if x[0] == 'left'][0][1]);
  countRight = len([x for x in mesh.bc if x[0] == 'right'][0][1]);
  countFront = len([x for x in mesh.bc if x[0] == 'front'][0][1]);
  countBack = len([x for x in mesh.bc if x[0] == 'back'][0][1]);
  countTop = len([x for x in mesh.bc if x[0] == 'top'][0][1]);
  countBottom = len([x for x in mesh.bc if x[0] == 'bottom'][0][1]);
  
  if countLeft != countRight:
    print("left: ", countLeft, " right: ", countRight)
   
  if countFront != countBack:
    print("front: ", countFront, " back: ", countBack)
   
  if countTop != countBottom:
    print("top: ", countTop, " bottom: ", countBottom)
  
  return countLeft, countRight, countFront, countBack, countTop, countBottom
  
## creates a 2D mesh of predefined geometry
def create_2d_mesh(type, x_res, y_res, width, opt_height = None, inclusion = None, inclusion_size = None, patch = None):
  
  assert(type == 'bulk2d' or type == 'cantilever2d' or type == 'cantilever2d_reinforced' or type == 'msfem_two_load' or type == 'two_load' or type.startswith('force_inverter') or type.startswith('gripper'))
  assert(inclusion == None or inclusion == "rect" or inclusion == "ball")
  assert(inclusion_size == None or inclusion_size <= 2.0)
  
  nx = x_res
  ny = x_res if y_res is None else y_res

  mesh = Mesh(nx, ny)
  
  # buld2d case
  ny = y_res if y_res != None else x_res
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
  if type == 'msfem_two_load' or type == 'two_load':
    width= 2.
    height = 1.
  if type == 'gripper_half' or type == 'force_inverter_half':
    height = 0.5
    ny = int(nx/2)
    
     
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
      e.type = TRIA3 if type == 'triangle_msfem' else QUAD4
      if type == 'cantilever2d_reinforced' and float(x) >= (28. / 30. * nx):
        e.region = 'reinforce'
      # strange: assure that x is meant to be 2.0 and y is meant to be 1.0 ?!  
      elif type == 'mbb_reinforced' and (x + 1 <= .015 * nx + 1e-5 or x >= 0.985 * nx - 1e-5 or y + 1 <= 0.03 * ny + 1e-5 or y >= 0.97 * ny - 1e-5):
        e.region = 'reinforce'
      elif type == 'force_inverter' and (((x == 0 or x < int(round(nx/100))) and (y <= 2*(ny-1)/25+1e-15 or y >= (ny-1)-2*(ny-1)/25-1e-15 or (y >= (ny-1)/2-(ny-1)/50-1e-15 and y <= (ny-1)/2+(ny-1)/50+1e-15))) or (x >= nx-int(round(nx/100)) and (y >= (ny-1)/2-(ny-1)/50-1e-15 and y <= (ny-1)/2+(ny-1)/50+1e-15))):
        e.region = 'reinforce'
      elif type == 'force_inverter_half' and (((x == 0 or x < int(round(nx/100))) and (y <= 4*(ny-1)/25+1e-15 or y >= (ny-1)-2*(ny-1)/50-1e-15)) or (x >= nx-int(round(nx/100)) and y >= (ny-1)-2*(ny-1)/50-1e-15)):
        e.region = 'reinforce'
      elif type == 'gripper' and (((x == 0 or x < int(round(nx/100))) and (y <= 2*(ny-1)/25+1e-15 or y >= (ny-1)-2*(ny-1)/25-1e-15 or (y >= (ny-1)/2-(ny-1)/50-1e-15 and y <= (ny-1)/2+(ny-1)/50+1e-15))) \
                                  or (x >= int(round(.7*(nx-1))) and ((y >= .35*(ny-1)-(ny-1)/50-1e-15 and y <= .35*(ny-1)+1e-15) or (y >= .65*(ny-1) and y <= .65*(ny-1)+(ny-1)/50+1e-15)))):
        e.region = 'reinforce'
      elif type == 'gripper_half' and (((x == 0 or x < int(round(nx/40))) and (y <= 4*(ny-1)/25+1e-15 or y >= (ny-1)-2*(ny-1)/50-1e-15 )) \
                                  or (x >= int(round(.7*(nx-1))) and ((y >= .7*(ny-1)-2*(ny-1)/40-1 and y <= .7*(ny-1)+1e-15)))):
        e.region = 'reinforce'
      elif type == 'gripper' and (x >= int(round(.7*(nx-1))) and y >= .35*(ny-1) and y <= .65*(ny-1)):
        e.region = 'void'
      elif type == 'gripper_half' and (x >= int(round(.7*(nx-1))) and y >= .7*(ny-1)):
        e.region = 'void'
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
      if inclusion != None:
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
    mesh.bc.append(("south", list(range(0, nx))))
    mesh.bc.append(("north", list(range((nx + 1) * (ny - 1), (nx + 1) * (ny - 1) + nx))))
    mesh.bc.append(("west", list(range(0, (nx + 1) * (ny - 1) + 1, nx + 1))))
    mesh.bc.append(("east", list(range(nx - 1, (nx + 1) * (ny - 1) + nx, nx + 1))))

    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx - 1]))
    mesh.bc.append(("left_upper", [(nx - 1) * (ny)]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny - 1) + nx - 1]))
    
    mesh.bc.append(("south_ghost", list(range(0, nx + 1))))
    mesh.bc.append(("north_ghost", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("west_ghost", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("east_ghost", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
  
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
    mesh.bc.append(("south", list(range(0, nx + 1))))
    mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
    mesh.bc.append(("pressure2", list(range(int(0.8 * nx), nx + 1))))
  
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh
  elif type == 'triangle_msfem':
    mesh.bc.append(("bottom", list(range(0, nx + 1))))
    mesh.bc.append(("top", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("left", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("right", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
  
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx + 1) * ny]))
    mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))
    return mesh
  elif type == 'two_load':
    mid = int((nx+1.)/2.)
    mesh.bc.append(("load1", list(range(mid,mid+1))))
    mesh.bc.append(("load2", list(range((nx+1)*ny + mid, (nx+1)*ny + mid+1)))) 
    mesh.bc.append(("support", list(range(0,1))))
    mesh.bc.append(("support", list(range(nx,nx+1))))
    mesh.bc.append(("support", list(range((nx+1)*ny,(nx+1)*ny+1,nx+1))))
    mesh.bc.append(("support", list(range((nx+1)*(ny+1)-1,(nx+1)*(ny+1)))))

 
  elif type == 'msfem_two_load':
    # lower/upper loads
    mid = int((nx+1.)/2.)
    off_x = int(0.05 * nx)
    off_y = int(0.1 * ny)
    if (nx+1)/2. % 2 == 0:
      mesh.bc.append(("load1", list(range(mid-off_x,mid + off_x))))
      mesh.bc.append(("load2", list(range((nx+1)*ny + mid-off_x, (nx+1)*ny + mid + off_x))))
    else:
      mesh.bc.append(("load1", list(range(mid-off_x,mid + 1 + off_x))))
      mesh.bc.append(("load2", list(range((nx+1)*ny + mid-off_x,(nx+1)*ny + mid + off_x + 1))))
    # support lower left
    mesh.bc.append(("support", list(range(0,off_x+1))))
    mesh.bc.append(("support", list(range(0,(nx+1)*off_y+1,nx+1))))
    # support lower right
    mesh.bc.append(("support", list(range(nx - off_x,nx+1))))
    mesh.bc.append(("support", list(range(nx,(nx+1)*(off_y+1),nx+1))))
    # support upper left
    mesh.bc.append(("support", list(range((nx+1)*ny, (nx+1)*ny+off_x+1))))
    mesh.bc.append(("support", list(range((nx+1)*ny-(nx+1)*off_y,(nx+1)*ny+1,nx+1))))
  elif type.startswith('force_inverter') or type.startswith('gripper'): 
    mesh.bc.append(("south", list(range(0, nx + 1))))
    mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
    if type.endswith('half'):
      factor = 2
      mesh.bc.append(("left_upper", [(nx+1)*ny]))
    else:
      factor = 1
      mesh.bc.append(("left_upper", numpy.arange(int(round((nx+1)*ny-2*(ny-1)/25*(nx+1))),(nx+1)*ny,nx+1)))
    mesh.bc.append(("left_lower", numpy.arange(0,int(round(factor*2*(ny-1)/25*(nx+1))),nx+1)))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))
  else: 
    mesh.bc.append(("south", list(range(0, nx + 1))))
    mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
    mesh.bc.append(("left_lower", [0]))
    mesh.bc.append(("right_lower", [nx]))
    mesh.bc.append(("left_upper", [(nx+1)*ny]))
    mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))

  print('width=' + str(width) + ' height=' + str(height) + ' dx=' + str(dx) + ' dy=' + str(dy))
  
  liu = numpy.unique(left_iface)
  
  if len(left_iface) > 0:
    mesh.bc.append(("left_iface", liu))
    print(str(len(liu)) + ' nodes (' + str(len(liu) / float(ny) * 100.) + '%) at the left interface for radius ' + str(inclusion_size))
  if len(right_iface) > 0:
    mesh.bc.append(("right_iface", numpy.unique(right_iface)))
  if len(lower_iface) > 0:
    mesh.bc.append(("lower_iface", numpy.unique(lower_iface)))
  if len(upper_iface) > 0:
    mesh.bc.append(("upper_iface", numpy.unique(upper_iface)))
  
  if second > 0:
    print(str(second) + ' elements of secondary region (' + str(100.0 * second / (nx * ny)) + '%)')

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

## creates a mesh of predefined geometry# inclusion is optional  
# data and threshold for sparse mesh from create_density. data is a numpy.array in 3D! 
# @param ext_mesh if given use it 
# @return a mesh, either ext_mesh or a newly created 
def create_3d_mesh(type, x_res, y_res = None, z_res = None, inclusion = None, inclusion_size = None, data = None, threshold = None, ext_mesh = None, scale = 1.0): 
  assert(type == "bulk3d" or type == "validation_test")

  nx = x_res  
   
  if type == "bulk3d": 
      ny = y_res if y_res != None else x_res 
      nz = z_res if z_res != None else x_res 
      width = scale
      height = scale*float(ny)/nx 
      depth = scale*float(nz)/nx 
  elif type == "cantilever3d": 
    ny = int(nx * (2./3.))
    nz = int(nx * (2./3.)) 
    width = 3.0 
    height = 2.0 
    depth = 2.0 
  elif type == "validation_test":
    ny = nx 
    nz = nx 
    width = 1.
    height = 1.
    depth = 1. 
    
     
  dx = width / nx 
  dy = height / ny 
  dz = depth / nz 

  assert(data is None or (len(data.shape) == 3 and data.shape[0] == nx and data.shape[1] == ny and data.shape[2] == nz)) 
  assert(data is None or threshold is not None) # theshold is mandatory when data is set 
  assert(not (data is None and threshold is not None)) # set threshold only when data is not set 
      
  mesh = Mesh(nx, ny, nz) if ext_mesh is None else ext_mesh
  
  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1
    
  print('width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz))
  
  
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
  mech_count = 0
  
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.density = 1.0 if data is None else data[x,y,z]
        e.type = HEXA8
        if inclusion == 'rect' and x >= nnx/2 * (1 - inclusion_size) and x < nnx/2 * (1 + inclusion_size)\
        and y >= nny/2 * (1 - inclusion_size) and y < nny/2 * (1 + inclusion_size) \
        and z >= nnz/2 * (1 - inclusion_size) and z < nnz/2 * (1 + inclusion_size) : 
          e.region = 'inner' 
          second += 1   
        elif inclusion == 'ball' and numpy.sqrt((x-nnx/2)**2 + (y-nny/2)**2 + (z-nnz/2)**2) <= nnx*inclusion_size: 
          e.region = 'inner' if not threshold or e.density > threshold else 'void'  
          second += 1
        elif type == "validation_test" and (y < int(0.1*ny) or y >= int(0.9*ny)):
          second += 1
          e.region = "non-design" if not threshold or e.density > threshold else 'void'       
        else: 
          e.region = 'mech' if not threshold or e.density > threshold else 'void' 
          # assign nodes 
          # ll = (nx+1)*y*(nx+1) * z + (nx+1) * y + x  # lowerleftfront 
        ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element 
        # start with upper-front-left counterclockwise in the x-z plane. Repeat in then lower plane 
        # e.nodes = ((ll+(nx+1), ll+1+(nx+1), ll+1+(nx+1)+((nx+1)*(ny+1)),ll+(nx+1)+((nx+1)*(ny+1)),ll, ll+1, ll+1+((nx+1)*(ny+1)),ll+((nx+1)*(ny+1))))   
        e.nodes = ((ll+nnx, ll+1+nnx, ll+1+nnx+(nnx*nny),ll+nnx+(nnx*nny),ll, ll+1, ll+1+(nnx*nny),ll+(nnx*nny))) 
        mesh.elements.append(e)

  if type == "validation_test":
    # create four support pins on bottom face
    side = (("support", []))
    mesh.bc.append(side)
    for z in range(0, int(0.2*nnz)+1):
      for x in range(0, int(0.2*nnx)+1):
        side[1].append((z * nny) * nnx + x)
    for z in range(0, int(0.2*nnz)+1):
      for x in range(int(0.8*nnx), nnx):
        side[1].append((z * nny) * nnx + x)
    for z in range(int(0.8*nnz), nnz):
      for x in range(int(0.8*nnx), nnx):
        side[1].append((z * nny) * nnx + x)
    for z in range(int(0.8*nnz), nnz):
      for x in range(0, int(0.2*nnx)+1):
        side[1].append((z * nny) * nnx + x)
    # create force region on top surface
    side = (("force", []))
    mesh.bc.append(side)
    for z in range(int(0.4*nnz), int(0.6*nnz)+1):
      for x in range(int(0.4*nnx), int(0.6*nnx)+1):
        side[1].append((z * nny + ny) * nnx + x)
  
  msg =  "dense resolution: " + str(nx) + " x " + str(ny) + " x " + str(nz) + " elements "
  msg += " -> " + str(mech_count) + " mech elements out of " + str(nx * ny * nz) + " (" + str(float(mech_count) / (nx * ny *nz) * 100.0) + " %)"
  msg += " with threshold " + str(threshold) 
  print(msg)
  
  if second > 0: 
    print(str(second) + ' elements of secondary region (' + str(100.0 * second / (nnx * nny * nnz)) + '%)') 
   
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
  
  mesh.bc.append(("west", list(range(nx+2, (nx+1) * ny, nx+1))))
  mesh.bc.append(("north", list(range(nx*ny, (nx+1) * ny-1, 1))))
  mesh.bc.append(("south", list(range(nx+2, 2*nx+1, 1))))
  #mesh.bc.append(("south",[nx,2*nx]))
  #mesh.bc.append(("south", [nx]))
  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx + 1) * ny]))
  mesh.bc.append(("right_upper", [(nx + 1) * (ny + 1) - 1]))

  if case == 'pipe_bend': 
    mesh.ne.append(('inlet', list(range(int(0.1 * nx * ny + eps), int(0.3 * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int((ny - 1) * nx + 0.7 * nx - eps), int((ny - 1) * nx + 0.9 * nx + eps)))))
  elif case == 'two_inlet_one_outlet':
    mesh.ne.append(('inlet', list(range(int((0.25 - 1. / 16) * nx * ny + eps), int((0.25 + 1. / 16) * nx * ny + nx + eps), nx))))
    mesh.ne.append(('inlet', list(range(int((0.75 - 1. / 16) * nx * ny + eps), int((0.75 + 1. / 16) * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int(0.375 * nx * ny - 1 + eps), int(0.625 * nx * ny - 1 + eps), nx))))
  elif case == 'two_inlet_two_outlet': 
    inletLength = 0.15 * ny 
    mesh.ne.append(('inlet',list(range(int(0.2*nx*ny+eps), int(0.2*nx*ny+eps + inletLength*nx),nx)) )) 
    mesh.ne.append(('inlet',list(range(int(0.8*nx*ny+eps-inletLength*nx), int(0.8*nx*ny+eps),nx)) )) 
    mesh.ne.append(('outlet',list(range(int(0.2*nx*ny+eps + nx-1), int(0.2*nx*ny+eps + inletLength*nx+ nx-1),nx)) )) 
    mesh.ne.append(('outlet',list(range(int(0.8*nx*ny+eps-inletLength*nx+ nx-1), int(0.8*nx*ny+eps+ nx-1),nx)) ))
  elif case == "pipe":
    mesh.ne.append(('inlet',list(range(nx,nx*(ny-1),nx))))   
    mesh.ne.append(('outlet',list(range(2*nx-1,nx*ny-1,nx))))  
  elif case == "diffuser":
    mesh.ne.append(('inlet',list(range(nx,nx*(ny-1),nx))))
    mesh.ne.append(('outlet',list(range(int(0.3*nx*ny+nx-1),int(0.7*nx*ny+nx-1),nx))))
  elif case == "low_in_high_out":
    mesh.ne.append(('inlet', list(range(int(0.1 * nx * ny + eps), int(0.3 * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int(0.7 * nx * ny - 1 + eps), int(0.9 * nx * ny - 1 + eps), nx))))
  else:
    print("unkwnon lbm case '" + case + "'")
    sys.exit(-1)
      
  return mesh

def create_backstep(x_res, y_res, z_res): 
  mesh = Mesh() 
 
  nx = x_res 
  ny = y_res if y_res != None else x_res 
  nz = z_res if z_res != None else x_res 
 
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
 
  print('width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz)) 
   
   
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
 
  print("Created " + str(second) + " obstacle elements") 
           
  mesh.bc.append(("left", list(range(0, (nnx*nny*z)+(nnx*ny)+1, nnx)))) 
  mesh.bc.append(("right", list(range(nx, (nnx*nny*nnz)+1, nnx)))) 
   
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
  mesh.bc.append(("back", list(range(0, (nx+1)*(ny+1))))) 
  mesh.bc.append(("front", list(range(nz*(nx+1)*(ny+1), (nz+1)*(nx+1)*(ny+1))))) 
 
 
  mesh.bc.append(("left_bottom_back",   [0])) 
  mesh.bc.append(("right_bottom_back",  [nx])) 
  mesh.bc.append(("left_top_back",      [nnx*ny])) 
  mesh.bc.append(("right_top_back",     [nnx*nny-1])) 
  mesh.bc.append(("left_bottom_front",  [nnx*nny*nz])) 
  mesh.bc.append(("right_bottom_front", [nnx*nny*nz+nx])) 
  mesh.bc.append(("left_top_front",     [nnx*nny*nz+nnx*ny])) 
  mesh.bc.append(("right_top_front",    [nnx*nny*nnz-1])) 
   
  for i in range(1,nz-1): 
    mesh.ne.append(('inlet',list(range(ymax*nx+nx+i*nx*ny,nx*ny*(i+1)-nx,nx)))) 
    mesh.ne.append(('outlet',list(range(2*nx+i*nx*ny-1,nx*ny*(i+1)-nx,nx)))) 
 
  return mesh 
 
   
def create_lbm3d(x_res, y_res, z_res, case, inclusion, inclusion_size): 
 
  nx = x_res 
  ny = y_res if y_res != None else x_res 
  nz = z_res if z_res != None else x_res 
 
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
 
  print('width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(dx) + ' dy=' + str(dy) + ' dz=' + str(dz)) 
   
   
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
           
  mesh.bc.append(("left", list(range(0, (nnx*nny*z)+(nnx*ny)+1, nnx)))) 
  mesh.bc.append(("right", list(range(nx, (nnx*nny*nnz)+1, nnx))))
  
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
      
  side = (("heat_bottom", []))
  mesh.bc.append(side)
  for z in range(1, nnz-1):
    for x in range(1,nnx-1):
      side[1].append((z*nny+1)*nnx+x) 
  
  side = (("heat_top", []))
  mesh.bc.append(side)
  for z in range(1, nnz-1):
    for x in range(1,nnx-1):
      side[1].append((z*nny+ny-1)*nnx+x)
      
  side = (("heat_back", []))
  mesh.bc.append(side)
  for y in range(1,nny-1):
    for x in range(nnx*nnz + 1 , nnx*nnz + nnx-1,1):
      side[1].append(x+nnx*y) 
  
  side = (("heat_front", []))
  mesh.bc.append(side)
  for y in range(1,nny-1):
    for x in range(nnx*nny*(nnz-2) + 1 , nnx*nny*(nnz-2) + nnx-1,1):
      side[1].append(x+nnx*y)   
   
  if case == 'pipe_bend': 
    area = 0.04 
    in_x = int(math.sqrt(area / (dx * dx))+eps) ## find out how many elements in one direction 
    dist_wall = 0.04 
    x_wall = int(math.sqrt(dist_wall / (dx * dx))+eps) #find out distance from wall 
    for i in range(nx*ny*nz-(x_wall+in_x)*nx*ny,nx*ny*nz-x_wall*nx*ny,nx*ny): 
      mesh.ne.append(('inlet',list(range(i+x_wall*nx,i+(x_wall+in_x)*nx,nx)))) 
    for i in range(0,in_x*nx*ny,nx*ny): 
      mesh.ne.append(('outlet', list(range((x_wall+1)*nx*ny-in_x-x_wall+i,(x_wall+1)*nx*ny-x_wall+i,1)))) 
  elif case == 'pipe': 
    for i in range (1,nz-1): 
      mesh.ne.append(('inlet', list(range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx)))) 
      mesh.ne.append(('outlet', list(range(int(i*nx*ny + 2*nx-1),int(i*nx*ny + 2*nx-1+nx*(ny-2)),nx)))) 
  elif case == 'two_inlet_one_outlet': 
    print("Not implemented yet!") 
  elif case == 'diffuser': 
    for i in range (1,nz-1): 
      mesh.ne.append(('inlet', list(range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx)))) 
    for i in range (int(0.3*nz),int(0.7*nz)): 
      mesh.ne.append(('outlet', list(range(int(i*nx*ny + 0.3*nx*ny + nx-1),int(i*nx*ny + 0.7*nx*ny+ nx-1),nx)))) 
  elif case == 'distributor': 
    center_x = nx / 2.0 
    center_y = ny / 2.0 
    center_z = nz / 2.0 
    width_x = 0.1 * nx  
    width_y = 0.1 * ny  
    width_z = 0.1 * nz  
    for i in range(int(round(center_z-width_z/2.0)),int(round(center_z+width_z/2.0))): 
      mesh.ne.append(('outlet',list(range(int(i*nx*nz-center_x+2.0*width_x),int(i*nx*nz-center_x+3.0*width_x),1)))) # top face 
      mesh.ne.append(('outlet',list(range(int((i-1)*nx*nz+center_x+2.0*width_x),int((i-1)*nx*nz+center_x+3.0*width_x),1)))) # bottom face 
    for i in range(int(center_z-int(width_z)),int(center_z+int(width_z))): 
      mesh.ne.append(('inlet',list(range(int(i*nx*ny+nx*int(center_y-width_y)),int(i*nx*ny+nx*int(center_y+width_y)),nx)))) #left face 
    for i in range(int(round(center_y-width_y/2.0)),int(round(center_y+width_y/2.0)),1): 
      mesh.ne.append(('outlet',list(range(int(i*nx+center_x+2.0*width_x),int(i*nx+center_x+3.0*width_x),1)))) #back face 
      mesh.ne.append(('outlet',list(range(int(nx*ny*(nz-1)+i*nx+center_x+2.0*width_x),int(nx*ny*(nz-1)+i*nx+center_x+3.0*width_x),1)))) #front face 
   
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
  if region_force != None:
    reg_force_nodes = hdf5_file['/Mesh/Groups/' + region_force + '/Nodes']
    mesh.bc.append((region_force, reg_force_nodes[:] - 1))
  # extract boundary force nodes from region_force if available
  elif region_support != None:
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
  idx = list(range(len(region)))
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
  print('read tetgenfile' + meshfile + '.1.ele')
  all_elements = numpy.loadtxt(meshfile + '.1.ele', dtype='int' , skiprows=1)
  print('read all_elements done')
  all_nodes = numpy.loadtxt(meshfile + '.1.node', skiprows=1)
  print('read all_nodes done')
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
  
def create_mesh_from_gmsh(meshfile,type):
  #from two_scale_tools import create_mesh_for_aux_cells, create_mesh_for_apod6
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
        print('Error: Gmsh format should be 2.2, result probably wrong')
    # read number of nodes
    elif count == 5:
      num_node = int(item[0])
    #add nodes
    elif count > 5 and count <= num_node + 5:
      nodes.append([float(item[1]),float(item[2]),float(item[3])])
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
  elif type == "aux_cells" or type == "base_cell":
    mesh = create_mesh_for_aux_cells(nodes,elem,1)
  else:
    print("Error: No correct type was selected! options: apod6, aux_cells")
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
      print('Warning mesh.bc type not handled!')
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
  print('Set regions and boundary nodes manually, default: design, non-design, force1,force2 and support')
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
  
def create_optistruct_mesh_from_cfs(meshfile,h5file):
  # manually select cfs hdf5 file
  print('Set regions and boundary nodes manually, default: design, non-design, non-design2, force1,force2, force3, support, support2, support3')
  mesh = create_mesh_from_hdf5(h5file, ['design','non-design','non-design2'], ['force1','force2','force3','support','support2','support3'])
  out = open(meshfile, "w")
  #design nodes and non-design nodes
  #out2 = open(meshfile + '.design', "w")
  #out3 = open(meshfile + 'non-desi"w")
  
  # optistruct header
  out.write('$\n')
  out.write('DESOBJ(MIN)=2\n')
  out.write('$\n')
  out.write('$$--------------------------------------------------------------\n')
  out.write('$$ HYPERMESH TAGS\n')
  out.write('$$--------------------------------------------------------------\n')
  out.write('$$BEGIN TAGS\n')
  out.write('$$END TAGS\n')
  out.write('$\n')
  out.write('BEGIN BULK\n')
  out.write('$$\n')
  out.write('$$  Stacking Information for Ply-Based Composite Definition\n')
  out.write('$$\n')
  out.write('$\n')
  out.write('$HMNAME OPTICONTROLS       1"optistruct_opticontrol"\n')
  out.write('$\n')
  out.write('DOPTPRM DESMAX  80      MINDIM  1.5     DISCRETE3.0     CHECKER 0\n')
  out.write('\n')
  out.write('\n')
  out.write('$HMNAME DESVARS        1Optimal\n')
  out.write('DTPL    1       PSOLID  1\n')
  out.write('+       STRESS  450.0\n')
  out.write('+       FATIGUE LIFE    80000.0\n')
  out.write('$$\n')
  out.write('$$  OPTIRESPONSES Data\n')
  out.write('$$\n')
  out.write('DRESP1  1       Vol_FracVOLFRAC PSOLID                                 1\n')
  out.write('DRESP1  2       COMPLIANCOMP\n')
  out.write('$$\n')
  out.write('$$  OPTICONSTRAINTS Data\n')
  out.write('$$\n')
  out.write('$\n')
  out.write('$HMNAME OPTICONSTRAINTS       1VOLUME\n')
  out.write('$\n')
  out.write('DCONSTR        1       1        0.35\n')
  out.write('\n')
  out.write('DCONADD        2       1\n')
  out.write('$$\n')
  out.write('$$  GRID Data\n')
  out.write('$$\n')
  # write nodes
  for i in range(len(mesh.nodes)):
    n = mesh.nodes[i]
    out.write('GRID%12d        '% (i+1) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
    #out.write('GRID    ' + '%-8d%-8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
  # Hexaeder elements
  out.write('$$\n')
  out.write('$$  SPOINT Data\n')
  out.write('$$\n')
  out.write('$\n')
  out.write('$  CHEXA Elements: First Order\n')
  out.write('$\n')
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
    elif e.type == HEXA8 and e.region == 'non-design2':
      out.write('CHEXA%11d%8d%8d%8d%8d%8d%8d%8d+\n'%(i+1,3,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      out.write('+       %8d%8d\n'%(n[6]+1,n[7]+1))
      #out.write('CHEXA   ' + '%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+E%-6d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1,i+1))
      #out.write('+E%-6d%-8d%-8d\n'%(i+1,n[6]+1,n[7]+1))
  out.write('$\n')
  out.write('$  CPENTA Elements 6-noded\n')
  out.write('$\n')
  # Wedge elements
  for i in range(len(mesh.elements)):
    e = mesh.elements[i]
    n = mesh.elements[i].nodes
    if e.type == WEDGE6 and e.region == 'design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == WEDGE6 and e.region == 'non-design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == WEDGE6 and e.region == 'non-design2':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,3,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
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

  out.write('$$------------------------------------------------------------------------------$\n')
  out.write('$$    HyperMesh name and color information for generic components               $\n')
  out.write('$$------------------------------------------------------------------------------$\n')
  out.write('$HMNAME COMP                   1"Lattice"\n')
  out.write('$HWCOLOR COMP                  2      56\n')
  out.write('$\n')
  out.write('$HMNAME COMP                   2"Wall"\n')
  out.write('$HWCOLOR COMP                  2      49\n')
  out.write('$\n')
  out.write('$HMNAME COMP                   3"Force Wall"\n')
  out.write('$HWCOLOR COMP                  3      52\n')
  out.write('$\n')
      
  #out.write('PSOLID         1       1\n')          
  #out.write('PSOLID         2       1\n')
  out.write('PSOLID  1       1       \n')          
  out.write('PSOLID  2       2       \n')
  out.write('PSOLID  3       1       \n')            
  #out.write('MAT1    1       1.00E0  0.34     0.0     0.785E-5  12.E-6                +M1\n')     
  #out.write('+M1     100.    -100.   100.\n')    
  #out.write('MAT1    2       7.00E4  0.34     0.0     0.785E-5  12.E-6                +M1\n')     
  #out.write('+M1     100.    -100.   100.\n')
  # Ti6Al4
  out.write('$ Ti6Al4\n')
  out.write('MAT1    %-8d%-.2e        %-8f%-8f%-8f%-8f%-8f\n'%(1,1.145E11,0.32,1.0,1.0,0.,1.))
  # Aluminum
  out.write('$ Aluminum\n')
  out.write('MAT2    %-8d%-.2e        %-8f%-8f%-8f%-8f%-8f\n'%(2,6.2E10,0.31,1.0,1.0,0.,1.)) 
  out.write('ENDDATA\n')
  
  out.close()    
  
def create_mesh_from_optistruct(meshfile,scale,type,offset = 0):
  # currently only used for apod6
  # read 3D optistruct mesh with hexa and wedge elements for apod6 got by M. Muir (12/2015)
  
  file = open(meshfile)
  inp = file.readlines()
  elem = []
  force1 = []
  force2 = []
  force3 = []
  forces = []
  support = []
  support2 = []
  support3 = []
  supports = []
  design = []
  nondesign = []
  count = 1
  num_node = 0
  num_elem = 0
  des = False
  nondes = False
  
  last_node_id = 0
  
  for i in range(len(inp)):
    if inp[i][0:8].strip() == 'GRID':
      last_node_id = int(inp[i][8:16].strip())
      
  nodes = [None] * (last_node_id+1)
  
  #rewind file
  file.seek(0)
  
  nodes_last_idx = 0 # count current last index of list 'nodes'
  
  for i in range(len(inp)):
    #item = str.split(inp[i])
    #if i < len(inp)-1:
    #  item_n = str.split(inp[i+1])
        
    if len(inp[i]) == 0:
      continue
    # read and check header
    if inp[i][0:8].strip() == 'GRID':
      if offset == -1: # set offset automatically
        offset = int(inp[i][8:16].strip()) 
      assert(offset >= 0)  
      #add nodes
      #x_str = 0
      #y_str = 0
      #z_str = 0
      #if len(item) > 3:
      #  # remove weird optistruct exponential function writing  
      #  for i in range(1,len(item[2][:])-1):
      #    if item[2][i] == '-':
      #      x_str = item[2][0:i-1] + str('e-') + item[2][i+1:8]
      #  for i in range(1,7):       
      #    if item[3][i] == '-':
      #      y_str = item[3][0:i-1] + str('e-') + item[3][i+1:8] 
      #  for i in range(9,15):
      #    if item[3][i] == '-':
      #      z_str = item[3][8:i-1] + str('e-') + item[3][i+1:16]     
      #else:
      #  for i in range(1,7):
      #    if item[2][i] == '-':
      #      x_str = item[2][0:i-1] + str('e-') + item[2][i+1:8]
      #  for i in range(9,15):
      #    if item[2][i] == '-':
      #      y_str = item[2][8:i-1] + str('e-') + item[2][i+1:16]
      #  for i in range(17,23):
      #    if item[2][i] == '-':
      #      z_str = item[2][16:i-1] + str('e-') + item[2][i+1:24]
      x = float(convert_optistruct_notation(inp[i][24:32],[0,8]))
      y = float(convert_optistruct_notation(inp[i][32:40],[0,8]))
      z = float(convert_optistruct_notation(inp[i][40:48],[0,8]))
      
#      map_mesh_nodeId[]
      nodes[int(inp[i][8:16].strip())-offset] = [x,y,z]
#       nodes.append([int(inp[i][8:16].strip()),x,y,z])
      
    elif inp[i][0:8].strip() == 'CTETRA':
      # read 3D tetra elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip())])
    elif inp[i][0:8].strip() == 'CPENTA':
      # read 3D wedge elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip()),int(inp[i][56:64].strip()),int(inp[i][64:72].strip())])
    elif inp[i][0:8].strip() == 'CHEXA':
      # read 3D hexahedron elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip()),int(inp[i][56:64].strip()),int(inp[i][64:72].strip()),int(inp[i][8:16].strip()),int(inp[i][16:24].strip())])
      i += 1     
    elif inp[i][0:8].strip() == 'CTRIA3':
      # read 2D triangle elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip())])
    elif inp[i][0:8].strip() == 'RBE2':
      support = []
      n = 0
      for k in range(5):
        support.append(int(inp[i][32+n:40+n].strip()))
        n += 8
      n = 0
      while inp[i+1][0] == '+':
        i += 1
        end = True
        k = 0
        n = 0
        while k < 8 and end:
          if inp[i][8+n:16+n].strip() == '':
            end = False
          else:
            support.append(int(inp[i][8+n:16+n].strip()))
          n += 8
          k +=1
      supports.append(support)
    elif inp[i][0:8].strip() == 'RBE3':
      force = []
      n = 0
      for k in range(2):
        force.append(int(inp[i][56+n:64+n].strip()))
        n += 8
      n = 0
      while inp[i+1][0] == '+':
        i += 1
        end = True
        k = 0
        n = 0
        while k < 8 and end:
          if inp[i][8+n:16+n].strip() == '':
            end = False
          else:
            force.append(int(inp[i][8+n:16+n].strip()))
          n += 8
          k += 1
      forces.append(force)

    #elif inp[i][0:8].strip() == '$HMMOVE' and inp[i][8:16].strip() == '6':
      # set flag for design material
    #  des = True
    #elif inp[i][0:8].strip() == '$HMMOVE' and inp[i][8:16].strip() == '7':
      # set flag for nondesign material
    #  des = False
    #  nondes = True
    #elif des == True and len(inp[i]) > 1:
      # nondesign elements
    #  item = str.split(inp[i])
    #  for j in range(1,9):
    #    design.append(item[j])
    #elif nondes == True and len(inp[i]) > 1:
      # add nondesign elements
    #  item = str.split(inp[i])
    #  for j in range(1,9):
    #    nondesign.append(item[j])
    else:
      des = False
      nondes = False
    count += 1  
  nodes = numpy.asarray(nodes)
  if type == "apod6":
    mesh = create_mesh_for_apod6(meshfile,nodes,elem)
  elif type == "cell_opt":
    # TUHH cell optimization
    mesh = create_mesh_for_aux_cells(nodes,elem,offset)
  elif type == "lufo_bracket":
    print('len support '+str(len(supports)))
    print(' len forces '+str(len(forces)))
    mesh = create_mesh_for_lufo_bracket(meshfile,nodes,elem,offset,forces,[], supports)
  else:
    print("Error: No correct type was selected! options: apod6, cell_opt, lufo_bracket")
#   write_gid_mesh(mesh, meshfile+".mesh",scale) # moved to create_mesh.py
  
  return mesh

# @profile
def voxelize_mesh_from_optistruct(filename,res):
  eps = 1e-3
  mesh = create_mesh_from_optistruct(filename, 1.0, 'cell_opt')
  
  array = numpy.zeros((res,res,res))
  minx, miny, minz, maxx, maxy, maxz = calc_min_max_coords(mesh)
  widthx = maxx-minx
  widthy = maxy-miny
  widthz = maxz-minz
  
  hx = widthx / res
  hy = widthy / res
  hz = widthz / res
  
  elems = mesh.elements
  
#   for e in elems:
#     coords = calc_barycenter(mesh,e)
#     i = int((coords[0] - minx)/hx - eps)
#     j = int((coords[1] - miny)/hy - eps)
#     k = int((coords[2] - minz)/hz - eps)
#     array[i,j,k] = 1
    
  for e in elems:
    barycenter = calc_barycenter(mesh,e)
    set_array_point(array,barycenter, hx, hy, hz, minx, miny, minz, 1)
    
    # calc longest side of triangle
    long_edge = calc_longest_edge(mesh,e)
#     print "longest edge: ",long_edge
#     print "barycenter: ", barycenter
#     print "hx: ", hx

    # if original mesh is too coarse for new resolution res**3, sample more points inside TETRA elem
    if long_edge > 0.9*hx: # assume hx = hy = hz
      points = []
      for idx,node in enumerate(e.nodes):
        points.append(mesh.nodes[node])
      
      tri = scipy.spatial.Delaunay(points)
      
      # create virtual cube around barycenter
      for x in numpy.arange(barycenter[0]-0.5*long_edge,barycenter[0]+0.51*long_edge,0.5*hx):
        for y in numpy.arange(barycenter[1]-0.5*long_edge,barycenter[1]+0.51*long_edge,0.5*hy):
          for z in numpy.arange(barycenter[2]-0.5*long_edge,barycenter[2]+0.51*long_edge,0.5*hz):
            if tri.find_simplex((x,y,z)) >= 0:
#               if inside_cube((x,y,z), points):
              set_array_point(array,(x,y,z), hx, hy, hz, minx, miny, minz, 2)
#     else:
#       print "long_edge:",long_edge
      
  minDim = [minx,miny,minz]
  maxDim = [maxx,maxy,maxz]
  meshNew = create_3d_mesh_from_array(array, True, widthx, widthy, widthz, minDim, maxDim)
  
  meshNew.nx = res
  meshNew.ny = res
  meshNew.nz = res
  
  meshNew = convert_to_sparse_mesh(meshNew)
  
  meshNew = add_nodes_for_periodic_bc(meshNew)
  
  validate_periodicity(meshNew)
  
  # moved to create_mesh.py
  # write_gid_mesh(meshNew, filename[:-4]+"_voxelized_res_" + str(res) + ".mesh")
  
  return meshNew
  

def convert_optistruct_notation(s,indexes):
  #remove weird optistruct exponential function writing  
  for i in range(indexes[0]+1,indexes[1]-1):
    if s[i] == '-':
      return s[indexes[0]:i-1] + str('e-') + s[i+1:indexes[1]]
  return s[indexes[0]:indexes[1]]
  


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
        

            



 
## create_mesh from unstructured mesh and assure that number of periodic boundary nodes are equal
#@param all_nodes can be unstructered
#@param elements can be unstructured
#@param offset optional offset for node numbering 
def create_mesh_for_aux_cells(all_nodes = [], elements = [],offset = 0.):
  mesh = Mesh()
  mesh.nodes = all_nodes
  
  min_diam_x = 1000000. 
  min_diam_y = 1000000. 
  min_diam_z = 1000000. 
  # insert elements
  for i in range(len(elements)):
      e = Element()
      e.nodes = (elements[i][1:]) 
      for k in range(len(e.nodes)):
        e.nodes[k] -= offset
      count = 0
      for k in range (len(e.nodes)):
        # determine the min_diam and max_diam of an element  
        if count + 1 >= len(e.nodes):
          min_diam_x = min(min_diam_x,abs(mesh.nodes[e.nodes[count]][0] - mesh.nodes[e.nodes[0]][0]))
          min_diam_y = min(min_diam_y,abs(mesh.nodes[e.nodes[count]][1] - mesh.nodes[e.nodes[0]][1]))
          min_diam_z = min(min_diam_z,abs(mesh.nodes[e.nodes[count]][2] - mesh.nodes[e.nodes[0]][2]))
        else:
          assert(count + 1 < len(e.nodes))
          min_diam_x = min(min_diam_x,abs(mesh.nodes[e.nodes[count]][0] - mesh.nodes[e.nodes[count+1]][0]))
          min_diam_y = min(min_diam_y,abs(mesh.nodes[e.nodes[count]][1] - mesh.nodes[e.nodes[count+1]][1]))
          min_diam_z = min(min_diam_z,abs(mesh.nodes[e.nodes[count]][2] - mesh.nodes[e.nodes[count+1]][2]))
        count += 1
      e.density = 1.
      e.region = 'mech'
      if len(e.nodes) == 4:
        e.type = TET4
      elif len(e.nodes) == 6:
        e.type = WEDGE6
      elif len(e.nodes) == 8:
        e.type = HEXA8
      mesh.elements.append(e)
      
  mesh = convert_to_sparse_mesh(mesh)
  
  mi_x, mi_y, mi_z, ma_x, ma_y, ma_z = calc_min_max_coords(mesh)
  
  delta = 1e-3
  
  mesh = add_nodes_for_periodic_bc(mesh, min_diam_x, min_diam_y, min_diam_z,delta)
  
  return mesh

# @param array to be written out
# @param singRegion do we want a mesh with only one region?
# @param minDim and maxDim contain for x,y,z direction the minimum/maximum coordinate
def create_3d_mesh_from_array(array,singRegion,widthx=1.0,widthy=1.0,widthz=1.0,minDim=[0.0,0.0,0.0],maxDim=[1.0,1.0,1.0]):
  nx, ny, nz = array.shape
  mesh = Mesh(nx,ny,nz)
  
  count = 0
  
  dx = widthx / nx
  dy = widthy / ny
  dz = widthz / nz
  
  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1

  for k in range(nnz):
    for j in range(nny):
      for i in range(nnx):
        mesh.nodes.append((minDim[0] + i * dx, minDim[1] + j * dy, minDim[2] + k * dz))
    
  for z in range(nz):    
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.type = HEXA8
        if (array[x][y][z] >= 0.0 and not singRegion):
          e.region = "mech" + str(int(array[x][y][z]))
          count += 1
        elif (array[x][y][z] >= 0.0 and singRegion):
          e.region = "mech"
          count += 1
        else:
          e.region = "void"
          
        ll = nnx*nny*z + nnx*y + x
        e.nodes = ((ll+nnx, ll+1+nnx, ll+1+nnx+(nnx*nny),ll+nnx+(nnx*nny),ll, ll+1, ll+1+(nnx*nny),ll+(nnx*nny)))
        mesh.elements.append(e)
  
  mesh = convert_to_sparse_mesh(mesh)  
  mesh = add_nodes_for_periodic_bc(mesh)
  
  return mesh

def create_2d_mesh_from_array(array):
  nx, ny, dummy = array.shape
  
  mesh = Mesh(nx,ny)
  
  dx = 1.0 / nx
  dy = 1.0 / ny
  
  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((x * dx, y * dy))
     
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.type = QUAD4
      if (array[x][y] > 0.0):
        e.region = "mech" + str(int(array[x][y]))
      else:
        e.region = "void"
       
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))
      mesh.elements.append(e)
  
  mesh.bc.append(("south", list(range(0, nx + 1))))
  mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
  mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
  mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx+1)*ny]))
  mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))
  
  return mesh
def create_validation_mesh(coords,nondes_coords, s1, s2, s3, ip_nx, grad, dir, scale,d_f,valid_position, valid_ring_position, type = "apod6",thres=0.0,csize = None,simp = None):
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
  print('calculating convex hull of non-design done')
  # choose validation for simp result or 2-scale result
  if simp is None:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1, s2, s3, dx,dy,dz)
  else:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1,None,None, dx,dy,dz)
  print('interpolation of thicknesses done')

  # lowest density = void density
  void = min(ip_data.ravel())
  
  # add points to fine mesh including shell
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  # where we want nodes
  nx = (int(delta[0] / dx) + 1)*dx_f
  ny = (int(delta[1] / dy) + 1)*dy_f
  nz = (int(delta[2] / dz) + 1)*dz_f

  print("nx,ny,nz = " + str(nx) + ", " + str(ny) + ", " + str(nz))
  print("des: ma,mi = " + str(ma) + " " + str(mi)) 
  print("nondes: ma,mi = " + str(nondes_max) + " " + str(nondes_min)) 
  #thickness of shell 1mm: tx,... is thickness of non-design shell
  if type == "apod6":
    tx = 0
    #ty = int(dy_f / (dy*1000))+1
    ty = dy_f / int(dy*1000) + 1
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
    print('chose a higher hom_samples or smaller cell_size such that also the smallest side gets discretized')
    exit()
  print("tx, ty, tz = " + str(tx) + ", " + str(ty) + ", " + str(tz))
  # add nodes including offset for shell 
  for z in range(-tz,nz+1+tz):
    for y in range(-ty, ny + 1 + ty):
      for x in range(-tx,nx+1+tx):
        mesh.nodes.append((mi[0] + float(x) * dx/dx_f, mi[1] +  float(y) * dy/dy_f, mi[2] +float(z) * dz/dz_f))
  print('inserting mesh.nodes done')
  nnx = nx + 2 * tx
  nny = ny + 2 * ty 
  nnz = nz + 2 * tz
  array = -1 * numpy.ones((nnx,nny,nnz))
  res = [dx_f,dy_f,dz_f]
  count = 0
  hole = -2. if type == "robot" else void
  for k in range(tz,nz-dz_f + 1 + tz,dz_f):
    for j in range(ty,ny-dy_f + 1 + ty,dy_f):
      for i in range(tx,nx -dx_f + 1 + tx,dx_f):    
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
  print('calculation of density array done')
  number = 0
  void3_count = 0
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
        center = calc_barycenter(mesh, e)
        
        if count >= 5:
          # test if is in convex hull of non-design nodes
          if in_hull(center, hull):
            if not valid_position(center, coords,offset):
              e.region = 'void1'
              #elif type == "robot" and array[x,y,z] > 0.9:
            else:
              e.region = 'non-design'
              number += 1
          else:
            e.region = 'void2'
        else:
          # test if is in convex hull of non-design nodes
          if in_hull(center, hull):
            if not valid_position(center, coords,offset):
              e.region = 'void1'
            #elif type == "robot" and array[x,y,z] > 0.9:
            elif valid_ring_position(center, coords,offset): # add material ring around holes 
              e.region = 'non-design2'  
            elif array[x,y,z] > 0.9:
              number += 1
              e.region = 'design'
            elif array[x,y,z] <= void + 1e-5:
              e.region = 'void3'
              void3_count += 1
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
  print('mesh has ' + str(number) + "design and non-design elements")
  print('volume = ' +str(float(number)/float(number + void3_count)))
  return mesh
