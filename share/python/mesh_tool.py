 #!/usr/bin/env python
import platform
import itertools
from enum import IntEnum
# from PIL import Image # load only in the function we need it!
import sys
import os
import numpy as np
import pickle

# cfs Elem::FEType for set_cfs_mesh()
class cfsFE(IntEnum):
  UNDEF   =  0
  POINT   =  1
  LINE2   =  2
  LINE3   =  3
  TRIA3   =  4
  TRIA6   =  5
  QUAD4   =  6
  QUAD8   =  7
  QUAD9   =  8
  TET4    =  9
  TET10   = 10
  HEXA8   = 11
  HEXA20  = 12
  HEXA27  = 13
  PYRA5   = 14
  PYRA13  = 15
  PYRA14  = 19
  WEDGE6  = 16
  WEDGE15 = 17
  WEDGE18 = 18
  POLYGON = 20
  POLYHEDRON = 21
  
# element types as in Ansys/gid (simInputMESH.cc -> AnsysType2ElemType)
# only partially coinciding with cfs's Elem::FEType
class Ansys(IntEnum):
  UNDEFINED = -1
  TRIA3 = 4
  QUAD4 = 6
  TET4 = 8
  TET10 = 9
  HEXA8 = 10
  WEDGE6 = 14
  LINE = 100

def cfs_type_to_ansys_type(cfs_type):
  if cfs_type == cfsFE.TRIA3:
    return Ansys.TRIA3
  if cfs_type == cfsFE.QUAD4:
    return Ansys.QUAD4
  if cfs_type == cfsFE.TET4:
    return Ansys.TET4
  if cfs_type == cfsFE.TET10:
    return Ansys.TET10
  if cfs_type == cfsFE.HEXA8:
    return Ansys.HEXA8
  if cfs_type == cfsFE.WEDGE6:
    return Ansys.WEDGE6
  if cfs_type == cfsFE.LINE2:
    return Ansys.LINE
  assert(False)

def ansys_type_to_cfs_type(ansys_type):
  if ansys_type == Ansys.TRIA3:
    return cfsFE.TRIA3
  if ansys_type == Ansys.QUAD4:
    return cfsFE.QUAD4
  if ansys_type == Ansys.TET4:
    return cfsFE.TET4
  if ansys_type == Ansys.TET10:
    return cfsFE.TET10
  if ansys_type == Ansys.HEXA8:
    return cfsFE.HEXA8
  if ansys_type == Ansys.WEDGE6:
    return cfsFE.WEDGE6
  if ansys_type == Ansys.LINE:
    return cfsFE.LINE2
  assert(False)
  
  
# type is Ansys type
def node_by_ansys_type(type):
  if type == Ansys.QUAD4:
    return 4
  if type == Ansys.HEXA8:
    return 8
  if type == Ansys.WEDGE6:
    return 6
  if type == Ansys.TRIA3:
    return 3
  if type == Ansys.LINE:
    return 2
  if type == Ansys.TET4:
    return 4
  if type == Ansys.TET10:
    return 10
  assert(False)


# dimension by Ansys-Type
def elem_dim_ansys(type):
  if type == Ansys.HEXA8 or type == Ansys.WEDGE6 or type == Ansys.TET4 or type == Ansys.TET10:
    return 3
  elif type == Ansys.LINE:
    return 1
  else:
    return 2



# gid element
class Element:
  def __init__(self, region = None, type = Ansys.UNDEFINED, density = 1):
    self.nodes = []  # list of zero based node indices. counter-clock wise
    self.region = region  # region name
    self.density = density  # from lower_bound to 1, not necessarily used, can be used to set the region name. see Mesh.set_density()
    self.type = type # this is for historical reasons the Ansys type!

  def dump(self):
    print(self.nodes)
    print(self.region)
    print(self.density)

# gid Mesh
class Mesh:
  # provides the structure of the mesh, does not fill it. Use create_2d_mesh or create_3d_mesh for it
  # @row_major is standard. Fastest variable is x, slowest is z (3d) or y (2d). If not created like this, some methos fail
  def __init__(self, nx = -1, ny = -1, nz = -1, row_major = True):
   self.nodes = [] # the nodes shall be an numpy array with 2 or three columns - some tools first usa a list and convert thend
   if nx > 0 and ny > 0 and nz <= 0:
     self.nodes = np.zeros(((nx+1) * (ny+1), 2)) 
   if nx * ny * nz > 0:
     self.nodes = np.zeros(((nx+1) * (ny+1) * (nz+1), 3))
           
   self.elements = []  # list of Element
   # list of boundary conditon nodes
   self.bc = []  # list of tupel (name, <list of zero based nodes>)
   # list of named element (save element in gd)
   self.ne = []  # list of tupel (name, <list of zero based elements>)
   self.nx = nx
   self.ny = ny
   self.nz = nz
   self.dx = None # with by nx. Set in create_[2/3]d_mesh
   self.dy = None
   self.dz = None
   
   self.row_major = row_major

  # returns the element index, respects row_major
  def element_idx(self, x, y, z = -1):
    assert self.nx > -1 # needs to be set
    if z == -1:
      if self.row_major:
        return y * self.nx + x
      else:
        return x * self.ny + y
    else: # no col_major in 3D yet
      return z * (self.ny*self.nx) + y * self.nx + x

  # determines mesh dimension by number of node coordinates
  def is2d(self):
    return len(self.nodes[0]) == 2 

  # return debug string for element by index
  def dump_elem(self, elem_idx):
    e = self.elements[elem_idx]
    nodes = ''
    for n in e.nodes:
      nodes += ' n_' + str(n) + '=' + str(self.nodes[n]) 
    return 'e=' + str(elem_idx) + ' r=' + e.region + ' d=' + str(e.density) + nodes
    
  # returns the x, y (,z) coordinates for the position, respect row_major.
  # this is the inverse of element_idx()
  # @param z_dummy if set do also return in 2D three values with the given z_dummy as value
  def element_idx_by_pos(self,pos, z_dummy = None):
    assert self.nx > -1 # needs to be set
    if self.nz <= 0:
      x = None; y = None
      if self.row_major:
        y = int(pos / self.nx)
        x = pos - y*self.nx
      else:
        x = int(pos / self.ny)
        y = pos - x * self.ny
      if z_dummy is None:
        return x, y
      else:
        return x, y, z_dummy 
    else:
      z = int(pos / (self.ny * self.nx))
      znynx = z*self.ny*self.nx 
      y = int((pos - znynx)/self.nx)
      x = pos - znynx - y*self.nx
      return x, y, z
        

  # return the elment at the given index position, Respects row_major
  def element(self, x, y, z = -1):
    return self.elements[self.element_idx(x,y,z)]

  # return the node index, respects row_major
  def node_idx(self,x,y,z=-1):
    assert self.nx > -1 # needs to be set
    if z == -1:
      if self.row_major:
        return y * (self.nx+1) + x
      else:
        return x * (self.ny+1) + y
    else: # no col_major in 3D yet
      return z * ((self.ny+1)*(self.nx+1)) + y * (self.nx+1) + x

  # return the node at the given index position. Respects row_major. 
  def node(self,x,y,z=-1):
    return self.nodes[self.node_idx(x,y,z)] 

  # compute the barycenter of the given element.
  # this is not in Element as we need the nodes
  def calc_barycenter(self,e):
    # there is no usage yet of nodes being a numpy array .. go ahead :)
    center = np.zeros(len(self.nodes[0]))
    for n in e.nodes:
      center += np.array(self.nodes[n])
    center /= len(e.nodes)
    return center
 
  # takes a 2d or 3d numpy array of density and set the Element.density accordingly. Also sets the region name
  # Element.region name is 'mech' if data value > threshold and 'void' for <= threshold
  # @multiregion set non-void to 'mech<int value>' of the data
  def set_density(self, data, threshold = 1e-6, multiregion = False):
    dx = None; dy = None; dz = -1 # data dimensions vs. self.nx dimension

    if len(data.shape) == 2:
      dx, dy = data.shape
    else:
      dx, dy, dz = data.shape
    # we access the density array as x,y or x,y,z which is not how we store the mesh (z,y,x)
    if data.size != self.nx * self.ny * abs(self.nz): 
      raise RuntimeError('incompatible data ' + str(data.shape) + ' for mesh (' + str(self.nx) + ', ' + str(self.ny) + ', ' + str(self.nz) + ')')
    
    # in 2d a numpy array has shape = (nx,ny). 
    # print(array) prints x as rows and y as columns: (6,4) from cantilever2d is 6 rows of 4 columns. 
    # bottom left image point (0,0) with show_density is data[0,0] and is the first element
    # top right image point is data[5,3] and is the last element
    # print(data) gives array([[1, 7, 13, 19],[2,8,14,20],[3,9,15,21],...,[6,12,18,24]]) als element numbers
    # you can test with with m = read_density('mech.density.xml','nr'); write_density_file('killme.density.xml', m/m.max()) 
    #print(data.shape)
    #print(data)
    assert self.row_major
    
    for z in range(1 if dz == -1 else dz):
      for y in range(dy):
        for x in range(dx):
          idx = self.element_idx(x,y,z)
          #print(x,y,z,idx)
          e = self.elements[idx]
          
          #e = self.element(z,y,x)
          #idx = dx * y + x
          #print(x,y,idx,data[x,y])
          #e = self.elements[idx]
          
          v = data[x,y] if dz == -1 else data[x,y,z] # for historical reaasons we access the matrix like this
          e.density = float(v)
          if v <= threshold:
            e.region = 'void'
          else:
            e.region = 'mech' + (str(int(v)) if multiregion else '')  
            
  # extract the density values from the elements. Note that Element.density is not always meaningfully set!
  # see get_density for ordering conventions
  # @return a 2D or 3D numpy array
  def get_density(self):
    assert self.nx * self.ny > 0
    assert self.row_major # otherwise check setting of data
     
    data = np.zeros((self.nx, self.ny)) if self.nz <= 0 else np.zeros((self.nx, self.ny, self.nz))
     
    for i, e in enumerate(self.elements):
      idx = self.element_idx_by_pos(i)         
      if self.nz <= 0:
        data[idx[0],idx[1]] = e.density
      else:          
        data[idx[0],idx[1], idx[2]] = e.density
         
    return data     
  # scale the nodal coordinates which are in m by default. 1e-2 makes them in cm
  # nodes shall be an numpy array
  def scale(self, scale):
    assert type(self.nodes) == np.ndarray # could easily be converted
    self.nodes *= scale 
            
  # translatate (add offset) the nodal coordinates 
  # nodes shall be a numpy array
  # @param offset needs to be a tuple, list or array of size dim 
  def move(self, offset):
    assert type(self.nodes) == np.ndarray # could easily be converted
    self.nodes += offset 
      
  # count number for given ansys type
  def count_elements(self, type):
    count = 0
    for e in self.elements:
      if e.type == type:
        count += 1
    return count
  
  # determine regions and their number of elements
  # @return list of tuples (name, num)
  def count_regions(self):
    dic = {}
    for e in self.elements:
      if e.region in dic:
        dic[e.region] += 1
      else:
        dic[e.region] = 1
        
    return [(k,v) for k,v in dic.items()] # note that order will be arbitrary   
        
  # rename region from old name to new name
  # @minimum_elements if less elements are replaced, an error is raised
  def rename_region(self, old_name, new_name, minimum_elements = 0):      
    cnt = 0
    for e in self.elements:   
      if e.region == old_name:  
        e.region = new_name
        cnt += 1
    if cnt < minimum_elements:
      raise ValueError("replaced " + str(cnt) + " element.region from '" + old_name + "' to '" + new_name + "' but minimum_elements is " + str(minimum_elements))
        
# node1/2: coordinates of respective node
def calc_edge_length(mesh,node1,node2):
  return np.linalg.norm(np.array(node1)-np.array(node2))

# calculates longest edge of element e
def calc_longest_edge(mesh,e):
  # so far only tested for 3D elements
  assert(elem_dim_ansys(e.type) == 3)
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

# return min_x, min_y, max_x, max_y for 2D and 6 values for 3D
# this is equivalent to a calc_bounding_box() method
# @dummy_z if given return this value for the 2D case to return also 6 values
def calc_min_max_coords(mesh, dummy_z = None):
  # make use of mesh.nodes being a numpy array
  min = np.amin(mesh.nodes,axis=0) 
  assert len(min) == 2 or len(min) == 3
  max = np.amax(mesh.nodes,axis=0)
  
  if mesh.nz <= 0:
    if dummy_z is None:
      return min[0], min[1], max[0], max[1]
    else:
      return min[0], min[1], dummy_z, max[0], max[1], dummy_z
  else:
    return min[0], min[1], min[2], max[0], max[1], max[2]
        
# finds the mesh barycenter coordinate out of the nodal coordinates
# the Mesh.calc_barycenter is for elements - not nice, next refactoring :)
# returns an arry of dim size
def calc_barycenter(mesh):        
  min = np.amin(mesh.nodes,0)
  max = np.amax(mesh.nodes,0)

  # m = a + .5 * (b-a) = .5 * (a+b)
  return .5 * (max + min)        
        
## creates a 2D mesh 
def create_2d_mesh(x_res, y_res = None, width = 1.0, height = None, pfem=False, row_major=True, triangles = False, silent=False):

  assert x_res is not None
  if width is None: # e.g. for the PythonMesher the external parameters might be None
    width = 1.0
   
  nx = x_res
  ny = y_res # might be None
  # if only one of height and y_res is set, the other is deduced.
  if y_res is None and height is None:
    ny = nx
    height = width
  elif y_res is None:
    ny = round(nx * (height/width)) # is int
  elif height is None: 
    height = round(width * (ny / nx),14)    
  assert ny is not None and height is not None  
    
  mesh = Mesh(nx, ny, row_major=row_major)

  mesh.dx = round(width / nx,15) # usually we round with 14 digits, but then the base shall be finer!
  mesh.dy = round(height / ny,15)

  # set nodes
  for y in range(ny+1):
    for x in range(nx+1):
      #print('x',x,'y',y,mesh.nodeindex(x, y))
      mesh.nodes[mesh.node_idx(x, y)] = [round(x * mesh.dx,14), round(y * mesh.dy,14)] # nodeindex is row_major sensitive

  # create volume (2d) elements
  for l1 in range(ny if row_major else nx):
    for l2 in range(nx if row_major else ny):
      x = l2 if row_major else l1
      y = l1 if row_major else l2
      e = Element()
      e.density = 1.0
      e.type = Ansys.TRIA3 if triangles else Ansys.QUAD4
      e.region = 'mech'
      # assign nodes
      ll = (nx+1) * y + x if row_major else (ny+1) * x + y # lowerleft
      if not triangles:
        if row_major:
          e.nodes = ((ll, ll + 1, ll + 1 + nx + 1, ll + nx + 1))
        else:
          e.nodes = ((ll, ll + ny + 1, ll + 1 + ny + 1, ll + 1))
        mesh.elements.append(e)
      else:          
        e.nodes = ((ll, ll + 1, ll + 1 + nx + 1))
        e2 = Element()
        e2.density = e.density
        e2.region = e.region
        e2.type = e.type
        e2.nodes = ((ll + 1 + nx + 1, ll + nx + 1, ll))
        mesh.elements.append(e)
        mesh.elements.append(e2)

  # create surface (1d) elements for pfem bounary conditions
  if pfem:
    for y in range(ny):
      for x in [0, nx]:
        e = Element()
        e.type = Ansys.LINE
        e.density = 1
        # assign nodes
        ll = (nx+1) * y + x if row_major else (ny+1) * x + y  # lowerleft
        e.nodes = ((ll,ll+nx+1 if row_major else ll+1))
        # set region for appropriate load case
        e.region = "left_surf" if x == 0 else "right_surf"
        mesh.elements.append(e)
  
    # create lines on top and bottom boundaries for pfem
    for x in range(nx):
      for y in [0, ny]:
        e = Element()
        e.type = Ansys.LINE
        e.density = 1
        # assign nodes
        ll = (nx+1) * y + x if row_major else (ny+1) * x + y  # lowerleft
        e.nodes = ((ll,ll+1 if row_major else ll+ny+1))
        # set region for appropriate load case
        e.region = "bottom_surf" if y == 0 else "top_surf"
        mesh.elements.append(e)

  # add named nodes for boundary conitions  
  if row_major:
    mesh.bc.append(("bottom", list(range(0, nx + 1))))
    mesh.bc.append(("top", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
    mesh.bc.append(("left", list(range(0, (nx + 1) * ny + 1, nx + 1))))
    mesh.bc.append(("right", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
    mesh.bc.append(("bottom_left", [0]))
    mesh.bc.append(("bottom_right", [nx]))
    mesh.bc.append(("top_left", [(nx+1)*ny]))
    mesh.bc.append(("top_right", [(nx+1)*(ny+1)-1]))
  else:
    mesh.bc.append(("bottom", list(range(0, (ny + 1) * nx + 1, ny + 1))))
    mesh.bc.append(("top", list(range(ny, (ny + 1) * (nx + 1), ny + 1))))
    mesh.bc.append(("left", list(range(0, ny + 1))))
    mesh.bc.append(("right", list(range((ny + 1) * nx, (ny + 1) * (nx + 1)))))
    mesh.bc.append(("bottom_left", [0]))
    mesh.bc.append(("bottom_right", [(ny+1)*nx]))
    mesh.bc.append(("top_left", [ny]))
    mesh.bc.append(("top_right", [(ny+1)*(nx+1)-1]))

  if not silent:
    print('width=' + str(width) + ' height=' + str(height) + ' dx=' + str(mesh.dx) + ' dy=' + str(mesh.dy))

  return mesh

# convenience function to create a mesh from an x,y np.array.
# for special options, call create_2d_mesh() and Mesh.set_density() manually
# threshold all > is mech, all below is void 
def create_2d_mesh_from_array(data, threshold = 1e-6):
  nx, ny = data.shape
  mesh = create_2d_mesh(nx, ny)
  mesh.set_density(data, threshold)
  return mesh

# create mesh from PIL image object (color or grayscale)
# all <= threshold becomes void
def create_2d_mesh_from_image(img, threshold):
  color = img.mode == 'RGB'
  pa = img.load() # pixel access
  nx, ny = img.size
  data = np.zeros((nx,ny)) # for some strange reasons np.array(pa) is still pa ?! Do it manually
  for y in range(ny): 
    for x in range(nx):  
      val = sum(pa[x,ny-1-y][0:3]) / 3.0 if color else pa[x,ny-1-y] # our y goes upwards (origin is left-lower). The image origin is left-upper
      data[x,y] = 1.0 - (val / 255.0)

  return create_2d_mesh_from_array(data, threshold)

## creates a mesh
# default is a unit cube with edge 1m, y/z_res and height/depth influence each other if only one is given
# elment.density is set to 1. Mesh.set_density() afterwards if you have density to set  
# @param offset is the nodal offset (x,y,z)
# note that you can use Mesh.scale() and Mesh.move() for an offset afterwards - just to keep the list of options clean
def create_3d_mesh(x_res, y_res = None, z_res = None, width = 1.0, height = None, depth = None, pfem=False):
  # we have x_res and width not None
  nx = x_res
  ny = y_res # might be None
  nz = z_res # might be None
  # if nothing is set, we do default
  if y_res is None and height is None:
    ny = nx
    height = width
  elif y_res is None:
    ny = round(nx * (height/width)) # is int
  elif height is None:
    height = round(width * (ny / nx),14)

  if z_res is None and depth is None:  
    nz = nx
    depth = width
  elif z_res is None:
    nz = round(nx * (depth/width))
  elif depth is None:
    depth = round(width * (nz / nx), 14)        
  
  assert ny is not None and nz is not None and height is not None and depth is not None  

  mesh = Mesh(nx, ny, nz)

  mesh.dx = round(width  / nx,15) # see create_2d_mesh for comment on 15 vs. 14 digits
  mesh.dy = round(height / ny,15)
  mesh.dz = round(depth  / nz,15)

  nnx = nx + 1
  nny = ny + 1
  nnz = nz + 1

  print('width=' + str(width) + ' height=' + str(height) + ' depth=' + str(depth) + ' dx=' + str(mesh.dx) + ' dy=' + str(mesh.dy) + ' dz=' + str(mesh.dz))

  # the coordinate system in Paraview is a right-hand sided coodrdinate system with z pointing to the viewer
  #
  #  y ^
  #    |
  # z (.)--> x
  #
  # These are the node numbers if we have only one element. The .mesh file will be transformed to 1-based
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
  #
  # definition of cube faces:
  # left: x=0, right x=1
  # bottom: y=0, top: y=1
  # back z=0, front: z=1
  for z in range(nnz):  # slowest variable
    for y in range(nny):
      for x in range(nnx):  # fastest variable
        mesh.nodes[mesh.node_idx(x, y, z)] = [round(x * mesh.dx,14), round(y * mesh.dy,14), round(z * mesh.dz,14)] # round coarser than dx, ...

  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        e = Element()
        e.density = 1.0
        e.type = Ansys.HEXA8
        e.region = 'mech' 
        # assign nodes
        # left: x=0, right x=1
        # bottom: y=0, top: y=1
        # back z=0, front: z=1
        ll = nnx*nny*z + nnx*y + x  #   left-bottom-back of current element
        # Local nodes definifing element: 0-1-3-2-4-5-7-6
        e.nodes = ( (ll, ll+1, ll+nnx+1,  ll+nnx, ll+nnx*nny, ll+(nnx*nny)+1, ll+nnx+(nnx*nny)+1,ll+nnx+(nnx*nny)) )
        mesh.elements.append(e)
  
  if pfem:
    for z in range(nz):
        for y in range(ny):
          for x in [0,nx]:
            e = Element()
            e.type = Ansys.QUAD4
            e.density = 1
            ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element
            e.nodes = ((ll,ll+nnx,ll+nnx*nny+nnx,ll+nnx*nny))
            mesh.elements.append(e)
            e.region = "left_surf" if x == 0 else "right_surf"
  
    for z in range(nz):
      for x in range(nx):
        for y in [0,ny]:
          e = Element()
          e.type = Ansys.QUAD4
          e.density = 1
          ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element
          e.nodes = ((ll,ll+nnx*nny,ll+nnx*nny+1,ll+1))
          e.region = "bottom_surf" if y == 0 else "top_surf"
          mesh.elements.append(e)
  
    for y in range(ny):
      for x in range(nx):
        for z in [0,nz]:
          e = Element()
          e.type = Ansys.QUAD4
          e.density = 1
          ll = nnx*nny*z + nnx*y + x  # lower-left-front of current element
          e.nodes = ((ll,ll+1,ll+nnx+1,ll+nnx))
          e.region = "back_surf" if z == 0 else "front_surf"
          mesh.elements.append(e)

  # naming faces of nodes of cube
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
  mesh.bc.append(("top_front",list(range(nnx*nny*nnz-nnx,nnx*nny*nnz))))
  mesh.bc.append(("top_left",list(range(nnx*nny-nnx,nnx*nny*nnz,nnx*nny))))
  mesh.bc.append(("top_right",list(range(nnx*nny-1,nnx*nny*nnz,nnx*nny))))
  mesh.bc.append(("back_left",list(range(0,nnx*nny-nnx+1,nnx))))
  mesh.bc.append(("back_right",list(range(nnx-1,nnx*nny,nnx))))
  mesh.bc.append(("front_left",list(range(nnx * nny * nz,nnx * nny * nz + nnx * ny+1,nnx))))
  mesh.bc.append(("front_right",list(range(nnx * nny * nz + nx,nnx * nny * nnz,nnx))))

  return mesh


# @see create_2d_mesh_from array. array with [x,y,z] indizces
def create_3d_mesh_from_array(array, threshold = 1e-6):
  nx, ny, nz = array.shape
  mesh = Mesh(nx,ny,nz)
  mesh.set_density(array, threshold)
  return mesh

      
     
# adds same number of boundary nodes on adjacent sides to assure periodic b.c
# same numbeer is minimal number. Works e.g. for meshes from stl
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
  for node in mesh.nodes:
    if abs(node[0] - mi_x) < min_diam_x + delta:
      left_c += 1
    elif abs(node[0] - ma_x) < min_diam_x + delta:
      right_c += 1
    elif abs(node[1] - mi_y) < min_diam_y + delta:
      bottom_c += 1
    elif abs(node[1] - ma_y) < min_diam_y + delta:
      top_c += 1
    elif abs(node[2] - mi_z) < min_diam_z + delta:
      back_c += 1
    elif abs(node[2] - ma_z) < min_diam_z + delta:
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
  
  # now add on eaach side dones up the the minum of both sides -> possible some will be ignored

  for i,node in enumerate(mesh.nodes):
    if abs(node[0] - mi_x) < min_diam_x + delta and left_c < lr_counter:
      left.append(i)
      left_c +=1
    elif abs(node[0] - ma_x) < min_diam_x + delta and right_c < lr_counter:
      right.append(i)
      right_c +=1
    elif abs(node[1] - mi_y) < min_diam_y + delta and bottom_c < bt_counter:
      bottom.append(i)
      bottom_c +=1
    elif abs(node[1] - ma_y) < min_diam_y + delta and top_c < bt_counter:
      top.append(i)
      top_c +=1
    elif abs(node[2] - mi_z) < min_diam_z + delta and back_c < bf_counter:
      back.append(i)
      back_c +=1
    elif abs(node[2] - ma_z) < min_diam_z + delta and front_c < bf_counter:
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


# name nodes on faces of bbox ("top","bottom",...)
# works for tet mesh and other meshees e.g. greated for stl
def name_bb_faces(mesh,bounds):
  nodes = mesh.nodes
  top = []
  bottom = []
  left = []
  right = []
  back = []
  front = []
  xmin,ymin,zmin,xmax,ymax,zmax = bounds
  for i, node in enumerate(nodes):
    if np.isclose(node[0],xmin,1e-4):
      left.append(i)
    elif np.isclose(node[0],xmax,1e-4):
      right.append(i)
    elif np.isclose(node[1],ymin,1e-4):
      front.append(i)
    elif np.isclose(node[1],ymax,1e-4):
      back.append(i)
    elif np.isclose(node[2],zmin,1e-4):
      bottom.append(i)
    elif np.isclose(node[2],zmax,1e-4):
      top.append(i)

  mesh.bc.append(("top",top))
  mesh.bc.append(("bottom",bottom))
  mesh.bc.append(("left",left))
  mesh.bc.append(("right",right))
  mesh.bc.append(("front",front))
  mesh.bc.append(("back",back))

  return mesh

# count prefined bc name and if pairs do not match, print  
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




# visualize a mesh for debugging purpose
def show_dense_mesh_image(mesh, binary, size):
  assert mesh.nx > 0
  from PIL import Image
  shape = (mesh.nx, mesh.ny)
  
  check_img = Image.new("RGB", shape, "white")
  check_pix = check_img.load()

  for x in range(mesh.nx):
    for y in range(mesh.ny):
      # print input_pix[x,y]
      e = mesh.element(x,y)
      val = 1 - e.density  # black is 0 in the image but 1 as density
      # print str(val) + " - " + str(barrier)
      show = (200, 10, 10) if binary else (int(val * 255), int(val * 255), int(val * 255))
      check_pix[x, mesh.ny - y - 1] = show if e.region == 'mech' else (10, 10, 200) if e.region == 'void' else (200, 10, 10)

  check_img = check_img.resize((size, int(mesh.ny * size / mesh.nx)), Image.NEAREST)
  check_img.show()


# @param mesh dense mesh (input)
# @return sparse mesh
def convert_to_sparse_mesh(dense):
  sparse = Mesh()

  # necessary 0-based nodes as unique set
  nns = set()

  # copy element, the indices of the nodes will be replaced later
  for e in dense.elements:
    if e.region != 'void':
      sparse.elements.append(pickle.loads(pickle.dumps(e))) # sending through pickle is faster than deepcopy
      for node in e.nodes:
        nns.add(node)

  # nns contains the required nodes uniquely and ordered
  # convert to nnl to access the indices
  nnl = list(nns)
  for i in nnl:
    sparse.nodes.append(dense.nodes[i])
    
  # we can now convert the nodes to the expected numpy array
  sparse.nodes = np.array(sparse.nodes)
  

  # next we need to correct the element node indices
  # the map is indexed by the dense numbering and contains the sparse indices or -1
  map = len(dense.nodes) * [-1]
  for i,idx in enumerate(nnl):
    map[idx] = i

  # now correct the element nodes
  for el in sparse.elements:
     newnodes = []  # el.nodes is a tuple, values cannot be replaces
     for node in el.nodes:
       newnodes.append(map[node])
       assert(node != -1)
     el.nodes = newnodes

  # finally handle the boundary conditions
  sparse.bc = []
  for bc in dense.bc:
    dnn = bc[1]  # dense nodes
    nodes = []
    for n in dnn:
#       print('old number '+str(dnn[n]) + ' new number '+str(map[dnn[n]]))
      if map[n] != -1:
        nodes.append(map[n])
    sparse.bc.append((bc[0], nodes))
    
    
    
  return sparse



# helper for write_ansys_mesh()
def write_ansys_elements(out, elements, dim):
  for i,e in enumerate(elements):  # write one based!
    if elem_dim_ansys(e.type) == dim:
      nodes = len(e.nodes)
      # we need the int value of the element type, not the enum name
      out.write(str(i + 1) + ' ' + str(int(e.type)) + ' ' + str(node_by_ansys_type(e.type)) + ' ' + e.region + "\n")

      # prepare for second order elements
      for n in range(node_by_ansys_type(e.type)):
        out.write(str(e.nodes[n] + 1) + ("\n" if n == node_by_ansys_type(e.type) - 1 else " "))  # write one based node numbers


# the ansys format is also used by the gid mesher export (not used any more) 
def write_ansys_mesh(mesh, filename, scale = 1):
  # Warning: mesh dimensions should be in [m]
  quad4  = mesh.count_elements(Ansys.QUAD4)
  hexa8  = mesh.count_elements(Ansys.HEXA8)
  wedge6 = mesh.count_elements(Ansys.WEDGE6)
  line   = mesh.count_elements(Ansys.LINE)
  tet4   = mesh.count_elements(Ansys.TET4)
  tet10  = mesh.count_elements(Ansys.TET10)
  tri3   = mesh.count_elements(Ansys.TRIA3)
  num_1d = line
  num_2d = quad4 + tri3
  num_3d = hexa8 + wedge6 + tet4 + tet10
  assert(num_1d + num_2d + num_3d == len(mesh.elements))
  #print('number of elements ' + str(num_1d + num_2d + num_3d))
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
  for bc in mesh.bc:
    bcn += len(bc[1])
  out.write('NumNodeBC ' + str(bcn) + '\n')
  nen = 0
  for ne in mesh.ne:
    nen += len(ne[1])
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
  out.write('Num tetra,quad   : ' + str(tet10) + '\n')
  out.write('Num brick        : ' + str(hexa8) + '\n')
  out.write('Num brick,quad   : 0\n')
  out.write('Num pyramid      : 0\n')
  out.write('Num pyramid,quad : 0\n')
  out.write('Num wedge        : ' + str(wedge6) + '\n')
  out.write('Num wedge,quad   : 0\n')

  out.write('\n[Nodes]\n')
  out.write('#NodeNr x-coord y-coord z-coord\n')
  for i,node in enumerate(mesh.nodes):  # write one based!
    out.write(str(i + 1) + "  " + str(node[0]/scale) + "  " + str(node[1]/scale))
    if dim == 3:
      out.write("  " + str(node[2]/scale) + "\n")
    else:
      out.write("  0.0\n")

  out.write('\n[1D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_ansys_elements(out, mesh.elements, 1)

  out.write('\n[2D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_ansys_elements(out, mesh.elements, 2)

  out.write('\n[3D Elements]\n')
  out.write('#ElemNr  ElemType  NrOfNodes  Level\n')
  out.write('#Node1 Node2 ... NodeNrOfNodes\n')
  write_ansys_elements(out, mesh.elements, 3)

  out.write('\n[Node BC]\n')
  out.write('#NodeNr Level\n')
  for bc in mesh.bc:
    for v in bc[1]:
      if v >= len(mesh.nodes):
        print(bc[0], v, ' larger', len(mesh.nodes)-1)
      assert(v >= 0)
      assert(v < len(mesh.nodes))
      out.write(str(v + 1) + " " + bc[0] + "\n") # bc[1][n]+1 is node number and bc[0] label

  out.write('\n[Save Nodes]\n')
  out.write('#NodeNr Level\n')
  out.write('\n[Save Elements]\n')
  out.write('#ElemNr Level\n')
  for ne in mesh.ne:
    for n in ne[1]:
      out.write(str(n + 1) + " " + ne[0] + "\n")

  out.write("\n \n")
  out.close()

# this is for being called via embedded python from cfs via the python input reader.
# cfs calls set_cfs_mesh() which creates a mesh and then calls this function to call the SimInputPython interface functions
# @see set_cfs_mesh()
def call_cfs_input_reader(mesh):
  #print('call_cfs_input_reader()')
  import cfs
  assert type(mesh.nodes) == np.ndarray # if still list, it could be easily converted
  # print(nodes)
  cfs.set_nodes(mesh.nodes);
  #print('nodes',nodes)

  # find regions and identify them by an id
  regions = list(set([e.region for e in mesh.elements]))
  #print(regions)
  cfs.set_regions(regions)
  
  # for add_elements we need to sort by type 
  total = len(mesh.elements)
  curr  = 0 # counter of provided elements to asign element number
  ansys_types = list(set([e.type for e in mesh.elements]))
  for at in ansys_types:
    cfs_type = ansys_type_to_cfs_type(at)
    #  array with the columns element number (1-based), region id, node numbers (1-based) corresponding to fe type */
    data = []
    for e in mesh.elements:
      if e.type == at:
        data.append(list(itertools.chain([curr+1,regions.index(e.region)],[n+1 for n in e.nodes]))) # 1-based element and node number
        curr += 1
    npd = np.array(data,dtype=int)     
    #print(total,cfs_type,npd)
    cfs.add_elements(total, cfs_type, npd)

  # named nodes
  for name, nodes in mesh.bc:
    if len(nodes) > 0:
      npn = np.array(nodes, dtype=np.uintc) # unsigned int 
      npn += 1 # make 1-based
      #print(name, npn) 
      cfs.add_named_nodes(name, npn)

  # named elements
  for name, elems in mesh.ne:
    if len(elems) > 0:
      npe = np.array(elems, dtype=np.uintc) # unsigned int 
      npe += 1 # make 1-based
      #print(name, npe) 
      cfs.add_named_elements(name, npe)

  
  # ignore named elemeents
  assert len(mesh.ne) == 0  


## this provides simple reference implementation for basic mesh generation with create_2d_mesh and create_3d_mesh
# to be called from cfs input python with options set with parameters type, x_res and optionally y_res, z_res, width, height, depth
# type is a subset of the create_mesh.py type attributes
def set_cfs_mesh(opt):
  print('\nset_cfs_mesh() from mesh_tool.py with', str(opt))
  if not ('type' in opt and 'x_res' in opt):
    raise "set at lright the options type and x_res" 

  type = opt['type']
  x_res  = int(opt['x_res'])
  y_res  = int(opt['y_res'])  if 'y_res' in opt else None
  z_res  = int(opt['z_res'])  if 'z_res' in opt else None
  width  = float(opt['width'])  if 'width' in opt else None
  height = float(opt['height']) if 'height' in opt else None
  depth  = float(opt['depth'])  if 'depth' in opt else None
  
  #print(x_res, y_res, z_res, width, height, depth)
  
  # 'bulk2d', 'bulk3d'
  if type == 'bulk2d':   
    mesh = create_2d_mesh(x_res, y_res, width=width, height=height)
  elif type == 'bulk3d':   
    mesh = create_3d_mesh(x_res, y_res, z_res, width=width, height=height, depth=depth)
  else:
    raise "unknown mesh type " + str(type) + ", only 'bulk2d' and 'bulk3d' known"  
      
  call_cfs_input_reader(mesh)


# creates a mesh from hdf5 file
def create_mesh_from_hdf5(hdf5_f, region, bcregions, region_force=None, region_support=None, threshold=0.):
  hdf5_file = h5py.File(hdf5_f, 'r')
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'].value  # for all regions
  # assume that region[0] is design, region[1] is non-design or void
  reg_elements_region = []
  for reg in region:
    reg_elements_region.append(hdf5_file['/Mesh/Regions/' + reg + '/Elements'].value)

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
    for bcreg in bcregions:
      bc_nodes = hdf5_file['Mesh/Groups/' + str(bcreg) + '/Nodes']
      mesh.bc.append((bcregions[i], bc_nodes[:] - 1))
  # insert nodes, usually nodes is a numpy array - transform it later
  for node in all_nodes:
    mesh.nodes.append(node)
  mesh.nodes = np.array(mesh.nodes)  

  # counter for regions
  idx = np.zeros(len(region))
  for el in all_elements:
    e = Element()
    e.nodes = el - 1
    #e.density = design_var[i]
    for j, regcount in idx:
      if regcount < len(reg_elements_region[j]):
        if i + 1 == reg_elements_region[j][regcount]:
          #if e.density >= threshold:
          e.region = region[j]
          #else:
          #  e.region = 'void'
          regcount += 1
    e.type = cfs_type_to_ansys_type(types[i])
    mesh.elements.append(e)
  return mesh

def create_mesh_from_tetgen(meshfile, region):
  print('read tetgenfile' + meshfile + '.1.ele')
  all_elements = np.loadtxt(meshfile + '.1.ele', dtype='int' , skiprows=1)
  print('read all_elements done')
  all_nodes = np.loadtxt(meshfile + '.1.node', skiprows=1)
  print('read all_nodes done')
  # all_faces = np.loadtxt(meshfile+'1.face',skiprows=1)
  # all_edges = np.loadtxt(meshfile+'1.edge',skiprows=1)


  # Create mesh 3D Tetrahedron
  mesh = Mesh()
  assert all_nodes.shape[1] == 3
  mesh.nodes = all_nodes

  for el in all_elements:
    e = Element()
    e.nodes = el[1:] - 1
    e.density = 1.
    e.region = region
    if len(e.nodes) == 4:
      e.type = Ansys.TET4
    else:
      assert(len(e.nodes) == 10)
      e.type = Ansys.TET10
    mesh.elements.append(e)
  return mesh


def create_mesh_from_gmsh(meshfile,regionnumbers=None,surfaceBCnumbers=[]):
  #from two_scale_tools import create_mesh_for_aux_cells, create_mesh_for_apod6
  # read 3D tetrahedron gmsh mesh
  if not meshfile.endswith(".msh"):
    meshfile = meshfile + ".msh"
  inp = open(meshfile).readlines()
  nodes = []
  if regionnumbers != None:
    regions = [[] for nums in regionnumbers]
  else:
    regions = []
  if (surfaceBCnumbers != None):
    bcs = [[] for nums in surfaceBCnumbers]
  count = 1
  num_node = 0
  num_elem = 0
  nodeListStart = -1
  elemListStart = -1
  for line in inp:
    item = str.split(line)
    # read and check header
    if count == 2:
      if float(item[0]) != 2.2:
        print('Error: Gmsh format should be 2.2, result probably wrong')
    # read starting line of nodes
    if '$Nodes' in line:
      nodeListStart = count+1
    # read ending line of nodes
    if '$Elements' in line:
      elemListStart = count+1
    # read number of nodes
    if count == nodeListStart:
      num_node = int(item[0])
    # read nodes
    elif count > nodeListStart and count <= (nodeListStart+num_node):
      nodes.append([float(item[1]),float(item[2]),float(item[3])])
    # read number of elements
    elif count == elemListStart:
      num_elem = int(item[0])
    # add elements
    elif count > elemListStart and count <= elemListStart+num_elem:
      el = [int(item[0])]
      for i in range(3+int(item[2]),len(item)):
        el.append(int(item[i]))
      # read 2D surface triangles (WARNING: unreliable, check results)
      if int(item[1]) == 2:
        for idx, bcnum in surfaceBCnumbers:
          for tag in range(int(item[2])):
            if int(item[3+tag]) == bcnum:
              for j in range(1,len(el)):
                bcs[idx].append(el[j])
              break
      # read 3D elements
      if int(item[1]) != 2:
        if regionnumbers != None:
          for idx, region in enumerate(regionnumbers):
            for tag in range(int(item[2])):
              if int(item[3+tag]) == region:
                regions[idx].append(el)
                break
        else:
          regions.append(el)
    count += 1

  # Create mesh
  # add nodes
  mesh = Mesh()
  mesh.nodes = np.array(nodes)

#   if regionnumbers == None:
#     regionnumbers = 'mech'
  # seems not to work if no region was specified
  if regionnumbers != None:
    print("regionnumbers:",regionnumbers)
    for j, region in enumerate(regionnumbers):
      for reg in regions[j]:
        e = Element()
        e.nodes = (reg[1:])
        for node in e.nodes:
          node -= 1
        e.density = 1.
        e.region = str(region)
        if len(e.nodes) == 4:
          e.type = Ansys.TET4
        elif len(e.nodes) == 6:
          e.type = Ansys.WEDGE6
        elif len(e.nodes) == 8:
          e.type = Ansys.HEXA8
        mesh.elements.append(e)
  else: # workaround
    for reg in regions:
      e = Element()
      e.nodes = (reg[1:])
      for node in e.nodes:
        node -= 1
      e.density = 1.
      e.region = "mech"
      if len(e.nodes) == 4:
        e.type = Ansys.TET4
      elif len(e.nodes) == 6:
        e.type = Ansys.WEDGE6
      elif len(e.nodes) == 8:
        e.type = Ansys.HEXA8
      mesh.elements.append(e)
  for idx, bcnum in enumerate(surfaceBCnumbers):
    mesh.bc.append((str(bcnum), list(set(bcs[idx]))))

## Manually add simple boundary conditions
  load = []
  support = []

  bounds = calc_min_max_coords(mesh)
#   print("xmin,ymin,zmin,xmax,ymax,zmax:",xmin,ymin,zmin,xmax,ymax,zmax)

  #mesh = add_bc_for_ppbox(mesh,(xmin,xmax,ymin,ymax,zmin,zmax))
  #mesh = add_bc_for_box_varel(mesh,bounds)
  #mesh = name_bb_faces(mesh,bounds)

  mesh = add_nodes_for_periodic_bc(mesh)

  write_ansys_mesh(mesh, meshfile[:-4]+".mesh")



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
  for el in elements:
      e = Element()
      e.nodes = el[2:]
      for node in e.nodes:
        node -= offset
      count = 0
      for _ in e.nodes:
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
        e.type = Ansys.TET4
      elif len(e.nodes) == 6:
        e.type = Ansys.WEDGE6
      elif len(e.nodes) == 8:
        e.type = Ansys.HEXA8
      mesh.elements.append(e)

  mesh = convert_to_sparse_mesh(mesh)

  mi_x, mi_y, mi_z, ma_x, ma_y, ma_z = calc_min_max_coords(mesh)

  delta = 1e-3

  mesh = add_nodes_for_periodic_bc(mesh, min_diam_x, min_diam_y, min_diam_z,delta)

  return mesh



def create_volume_mesh_from_stl(stlName,type=None,write_vtk=True):
  assert(stlName.endswith(".stl"))
  # -p Tetrahedralizes a piececwise linear complex
  # -k Outputs mesh to .vtk file for viewing by Paraview
  # -O3 optimization level 3
  command = "tetgen -pk -O3" if write_vtk else "tetgen -p"
  cfs_utils.execute(command + " " + stlName)
  mesh = create_mesh_from_tetgen(stlName[:-4],"mech")
  if type == "box_varel":
    bounds = calc_min_max_coords(mesh)
    mesh = add_bc_for_box_varel(mesh,bounds)
    mesh = name_bb_faces(mesh,bounds)
  else:
    add_nodes_for_periodic_bc(mesh)
    validate_periodicity(mesh)

  return mesh

def create_volume_mesh_with_gmsh(stlName):
  baseName = stlName[:-4]
  # write .geo file for gmsh
  geoName = baseName + ".geo"
  out = open(geoName,"w")
  out.write("Merge '" + stlName + "';\n")
  out.write("// add a geometrical volume \n")
  out.write("Surface Loop(1) = {1};\n")
  out.write("Volume(1) = {1};\n")
  out.flush()
  out.close()

  # -3: tetrahedralize (3D)
  # -optimize: use netgen's mesh optimization tool
  command = "gmsh -3 " + geoName
  cfs_utils.execute(command)
  create_mesh_from_gmsh(baseName)

