#!/usr/bin/env python

import Image, sys, os, copy, numpy, math

# writes a dense two region mesh
# def write_dense_mesh(pixels, size, file, threshold):

# element types as in gid (simInputMESH.cc -> AnsysType2ElemType)
QUAD4 = 6
HEXA8= 10
WEDGE6 = 14

def nodes_by_type(type):
  if type == QUAD4:
    return 4
  if type == HEXA8:
    return 8
  if type == WEDGE6:
    return 6
  assert(False)

# gid element
class Element: 
  def __init__(self):
    self.nodes = [] # list of zero based node indices. counter-clock wise
    self.region = None # region name
    self.density = 1 # from lower_bound to 1, not necessarily used
    self.stiff1 = 0
    self.stiff2 = 0
    self.stiff3 = 0
    self.rotAngle = 0
    self.type = -1
    
  def dump(self):
    print self.nodes
    print self.region
    print self.density     


# gid Mesh
class Mesh:
  def __init__(self):
   self.nodes = []    # list 2d tupels (float, float) or 3d tuples
   self.elements2d = [] # list of Element
   self.elements3d = [] # list of Element
   self.bc = []       # list of tupel (name, <list of zero based nodes>)  

def show_dense_mesh_image(mesh, shape, binary, size):
  check_img = Image.new("RGB", shape, "white")
  check_pix = check_img.load()

  nx, ny = shape
  
  for x in range(nx):
    for y in range(ny):
      #print input_pix[x,y]
      e = mesh.elements2d[x * ny + y]
      val = 1-e.density # black is 0 in the image but 1 as density
      # print str(val) + " - " + str(barrier)
      show = (200,10,10) if binary else (int(val*255),int(val*255),int(val*255)) 
      check_pix[x,ny-y-1] = show if e.region == 'mech' else (10,10,200)
      
  check_img = check_img.resize((size, ny * size/nx))    
  check_img.show()     


def create_dense_mesh_img(input_img, mesh, threshold, scale, rhomin, shearAngle):
  input_pix = input_img.load()
  nx, ny = input_img.size
  create_dense_mesh(input_pix, nx, ny, mesh, threshold, scale, rhomin,True,1,shearAngle)

def create_dense_mesh_density(numpy_array, mesh, threshold, scale, rhomin,multi_d=1):
  if multi_d == 1:
    nx, ny = numpy_array.shape
  else:
    nx,ny,nz,m = numpy_array.shape
  create_dense_mesh(numpy_array, nx, ny, mesh, threshold, scale, rhomin,False,multi_design = multi_d,shearAngle = 0.0)
  
def create_dense_mesh(input_array, nx, ny,  mesh, threshold, scale, rhomin,img = True,multi_design=1,shearAngle=0):
  # convert angle to rad and check for feasibility
  angle = shearAngle/180 * math.pi
  if (abs(angle) > math.pi/2 - 1e-6):
    print 'angle has to be between -pi/2 + 1e-6 and pi/2 - 1e-6'
    return 0 
  input_pix = input_array
  dx = scale/nx/math.cos(angle)
  dy = scale/ny
    # create mesh.nodes
  for y in range(ny + 1):
    for x in range(nx + 1):
      if angle == 0.0:
        mesh.nodes.append((x * dx, y * dy))
      else:
        x_Coord = round(x * dx - y * dy * math.tan(angle), 8)
        if abs(x_Coord) < 1e-8:
          x_Coord = 0.0
        mesh.nodes.append((x_Coord,y * dy))
  # print mesh.nodes 
  mech_count = 0
  for x in range(nx):
    for y in range(ny):
      e = Element()
      e.type = QUAD4
      if img:
        # convert to black is one and white = 0
        e.density = 1.0 - (input_pix[x,y] / 255.0)
      else:
        if multi_design == 1:
          e.density = input_pix[x,y]
        else:
          e.stiff1 = input_pix[x,y,0,0]
          e.stiff2 = input_pix[x,y,0,1]
          if multi_design == 3:
            e.rotAngle = input_pix[x,y,0,2]
      if multi_design == 1:
        if e.density < rhomin:
          e.density = rhomin
      else:
        if e.stiff1 < rhomin:
          e.stiff1 = rhomin
        elif e.stiff2 < rhomin:
          e.stiff2 = rhomin
      if multi_design == 1:    
        if float(e.density) >= float(threshold):
          e.region = 'mech'
          mech_count += 1
        else:
          e.region = 'void'
      else:
        if (float(e.stiff1) >= float(threshold)) | (float(e.stiff2) >= float(threshold)):  
          e.region = 'mech'
          mech_count += 1
        else:
          e.region = 'void'
      # assign nodes
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))
      mesh.elements2d.append(e)
      # e.dump()
  
  mesh.bc.append(("bottom", range(0, nx+1)))
  mesh.bc.append(("top", range((nx+1)*ny, (nx+1)*(ny+1))))
  mesh.bc.append(("left", range(0, (nx+1)*ny+1, nx+1)))
  mesh.bc.append(("right", range(nx, (nx+1)*(ny+1), nx+1)))
  
  #print mesh.bc

  print "dense resolution: " + str(nx) + " x " + str(ny) + " elements (" + str(scale) + "m x " + str(float(ny)/nx * scale) + "m)",
  print " -> " + str(mech_count) + " mech elements out of " + str(nx*ny) + " (" + str(float(mech_count) / (nx*ny) * 100.0) + " %)",
  print " with threshold " + str(threshold) 

# @param mesh dense mesh (input)
# @return sparse mesh
def convert_to_sparse_mesh(dense):
  sparse = Mesh()

  # necessary 0-based nodes as unique set
  nns = set()
  
  # copy element, the indices of the nodes will be replaced later
  for i in range(len(dense.elements)):
    e = dense.elements[i]
    if e.region == 'mech':
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
     newnodes = [] # el.nodes is a tuple, values cannot be replaces
     for n in range(len(el.nodes)):
       newnodes.append(map[el.nodes[n]])
       assert(el.nodes[n] <> -1)
     el.nodes = newnodes

  # finally handle the boundary conditions
  sparse.bc = []
  for b in range(len(dense.bc)):
    bc = dense.bc[b]
    dnn = bc[1] # dense nodes
    nodes = []
    for n in range(len(dnn)):
      if map[dnn[n]] <> -1:
        nodes.append(map[dnn[n]])
    sparse.bc.append((bc[0], nodes))
    print "sparse boundary condition '" + bc[0] + "' with " + str(len(nodes)) + " nodes"

  return sparse

def count_elements(elements, type):
  count = 0
  for e in elements:
    if e.type == type:
      count += 1
  return count

def write_gid_elements(out, elements):
  for i in range(len(elements)): # write one based!
    e = elements[i]
    nodes = len(e.nodes)
    out.write(str(i+1) + ' ' + str(e.type) + ' ' + str(nodes) + ' ' + e.region  + "\n")
    
    # prepare for second order elements
    for n in range(nodes):
      out.write(str(e.nodes[n]+1) + ("\n" if n == len(e.nodes) - 1 else " ")) # write one based node numbers


def write_gid_mesh(mesh, filename):
  dim = 3 if len(mesh.elements3d) > 0 else 2
  quad4  = count_elements(mesh.elements2d, QUAD4)
  hexa8  = count_elements(mesh.elements3d, HEXA8)
  wedge6 = count_elements(mesh.elements3d, WEDGE6)
  assert(quad4 == len(mesh.elements2d))
  assert(hexa8 + wedge6 == len(mesh.elements3d))
  
  out = open(filename, "w")
  
  out.write('[Info]\n')
  out.write('Version 1\n')
  out.write('Dimension ' + str(dim) + '\n')
  out.write('NumNodes ' + str(len(mesh.nodes)) + '\n')
  out.write('Num3DElements ' + str(len(mesh.elements3d)) + '\n')
  out.write('Num2DElements ' + str(len(mesh.elements2d)) + '\n')
  out.write('Num1DElements 0\n')
  bcn = 0
  for i in range(len(mesh.bc)):
    bcn += len(mesh.bc[i][1])
  out.write('NumNodeBC ' + str(bcn) + '\n')
  out.write('NumSaveNodes 0\n')
  out.write('NumSaveElements 0\n')
  out.write('Num 2d-line      : 0\n')
  out.write('Num 2d-line,quad : 0\n')
  out.write('Num 3d-line      : 0\n')
  out.write('Num 3d-line,quad : 0\n')
  out.write('Num triangle     : 0\n')
  out.write('Num triangle,quad: 0\n')
  out.write('Num quadr        : ' + str(quad4) + '\n')
  out.write('Num quadr,quad   : 0\n')
  out.write('Num tetra        : 0\n')
  out.write('Num tetra,quad   : 0\n')
  out.write('Num brick        : ' + str(hexa8) + '\n')
  out.write('Num brick,quad   : 0\n')
  out.write('Num pyramid      : 0\n')
  out.write('Num pyramid,quad : 0\n')
  out.write('Num wedge        : ' + str(wedge6) + '\n')
  out.write('Num wedge,quad   : 0\n')
  
  out.write('\n[Nodes]\n')
  out.write('#NodeNr x-coord y-coord z-coord\n')
  for i in range(len(mesh.nodes)): # write one based!
    out.write(str(i+1) + "  " + str(mesh.nodes[i][0]) + "  " + str(mesh.nodes[i][1]))
    if dim == 3:
      out.write("  "  + str(mesh.nodes[i][2]) + "\n")
    else:
      out.write("  0.0\n")

  out.write('\n[1D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  
  out.write('\n[2D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_gid_elements(out, mesh.elements2d)
 
  out.write('\n[3D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_gid_elements(out, mesh.elements3d)
  
  out.write('\n[Node BC]\n')
  out.write('#NodeNr Level\n')
  for b in range(len(mesh.bc)):
    bc = mesh.bc[b]
    for n in range(len(bc[1])):
      out.write(str(bc[1][n]+1) + " " + bc[0] + "\n")
      
  out.write('\n[Save Nodes]\n')
  out.write('#NodeNr Level\n')
  out.write('\n[Save Elements]\n')
  out.write('#ElemNr Level\n\n')

  out.close()
  
## creates a mesh of predefined geometry
def create_cantilever2d_mesh(type, resolution):
  mesh = Mesh()
  
  width = 3.0
  height = 2.0
  
  nx = resolution
  ny = int(nx * (2./3.))
  
  dx = width / nx
  dy = height / ny 

  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((x * dx, y * dy))
 
  # print mesh.nodes 
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.density = 1.0
      e.type = QUAD4
      if type == 'cantilever2d_reinforced' and float(x) >= (28./30. * nx):
        e.region = 'reinforce'
      else:
        e.region = 'mech'

      # assign nodes
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))
            
      mesh.elements2d.append(e)
  
  mesh.bc.append(("south", range(0, nx+1)))
  mesh.bc.append(("north", range((nx+1)*ny, (nx+1)*(ny+1))))
  mesh.bc.append(("west", range(0, (nx+1)*ny+1, nx+1)))
  mesh.bc.append(("east", range(nx, (nx+1)*(ny+1), nx+1)))

  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx+1)*ny]))
  mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))
  
  return mesh

## creates a mesh of predefined geometry
def create_mbb_mesh(type, resolution):
  mesh = Mesh()
  
  width = 2.0
  height = 1.0
  
  nx = resolution
  ny = int(nx * 0.5)
  
  dx = width / nx
  dy = height / ny 
  
  e = 1e-4


  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((x * dx, y * dy))
 
  # print mesh.nodes 
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.type = QUAD4
      e.density = 1.0
      # if type == 'mbb_reinforced' and (float(x) <= (.03 * nx + e) or float(x) >= (.97 * nx - e) or float(y) <= (.03 * ny + e) or float(y) >= (.97 * ny - e)):
      if type == 'mbb_reinforced' and (x+1 <= .015 * nx + 1e-5 or x >= 0.985 * nx - 1e-5 or y+1 <= 0.03 * ny + 1e-5 or y >= 0.97 * ny - 1e-5):
          e.region = 'reinforce'
      else:
        e.region = 'mech'

      #if y == 0:
      #  print "x=" + str(x) + " -> " + str(.015 * nx) + " r=" + e.region


      # assign nodes
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))
            
      mesh.elements2d.append(e)
  
  mesh.bc.append(("south", range(0, nx+1)))
  mesh.bc.append(("north", range((nx+1)*ny, (nx+1)*(ny+1))))
  mesh.bc.append(("west", range(0, (nx+1)*ny+1, nx+1)))
  mesh.bc.append(("east", range(nx, (nx+1)*(ny+1), nx+1)))

  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx+1)*ny]))
  mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))
  
  return mesh
