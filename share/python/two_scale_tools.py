#!/usr/bin/env python
from mesh_tool import *

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
  r2 = 0.0173
  r3 = 0.0058
  delta = 0.0012
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
    
    r1 = 19.51
    r2 = 17.31
    r3 = 7.1
    delta = 0.
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
          elif f2 == True and abs(coord[1] - upper_bound) < dy:
            force2.append(e.nodes[j])
          elif f3 == True and abs(coord[1] - upper_bound) < dy:
            force3.append(e.nodes[j])
          elif sp2 == True and abs(coord[1] - upper_bound) < dy:
            support2.append(e.nodes[j])    
          elif sp3 == True and abs(coord[1] - upper_bound) < dy:
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
        # determine the min_diam and max_diam of an element  
        if count + 1 == len(e.nodes):
          min_diam_x = min(min_diam_x,abs(mesh.nodes[e.nodes[count]][0] - mesh.nodes[e.nodes[0]][0]))
          min_diam_y = min(min_diam_y,abs(mesh.nodes[e.nodes[count]][1] - mesh.nodes[e.nodes[0]][1]))
          min_diam_z = min(min_diam_z,abs(mesh.nodes[e.nodes[count]][2] - mesh.nodes[e.nodes[0]][2]))
        else:
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
  
  print "min = [" + str(mi_x) + ", " +  str(mi_y) + ", " +  str(mi_z) + "]"
  print "max = [" + str(ma_x) + ", " +  str(ma_y) + ", " +  str(ma_z) + "]"
  delta = 1e-4
  top = []
  bottom = []
  left = []
  right = []
  front = []
  back = []
  # counter necessary to have the same number of left and right nodes, ...
  left_c = 0
  right_c = 0
  bottom_c = 0
  top_c = 0
  back_c = 0
  front_c = 0
  min_diam_x = 1e-3
  min_diam_y = 1e-3
  min_diam_z = 1e-3
  # count number of boundary nodes per region
  for i in range(len(mesh.nodes)):
    if abs(mesh.nodes[i][0] - mi_x) < min_diam_x - delta:
      left_c += 1
    elif abs(mesh.nodes[i][0] - ma_x) < min_diam_x - delta:
      right_c += 1
    elif abs(mesh.nodes[i][1] - mi_y) < min_diam_y - delta:
      bottom_c += 1 
    elif abs(mesh.nodes[i][1] - ma_y) < min_diam_y - delta:
      top_c += 1
    elif abs(mesh.nodes[i][2] - mi_z) < min_diam_z - delta:
      back_c += 1
    elif abs(mesh.nodes[i][2] - ma_z) < min_diam_z - delta:
      front_c += 1
      
  lr_counter = min(left_c,right_c)
  bt_counter = min(bottom_c,top_c)
  bf_counter = min(back_c,front_c)
  left_c = 0
  right_c = 0
  top_c = 0
  bottom_c = 0
  back_c = 0
  front_c = 0
  print 'lr_counter = ' + str(lr_counter) + ', bt_counter = ' + str(bt_counter) + ', bf_counter = ' + str(bf_counter)
  # insert the same number of nodes to each region pair
  for i in range(len(mesh.nodes)):
    if abs(mesh.nodes[i][0] - mi_x) < min_diam_x - delta and left_c < lr_counter:
      left.append(i)
    elif abs(mesh.nodes[i][0] - ma_x) < min_diam_x - delta and right_c < lr_counter:
      right.append(i)
    elif abs(mesh.nodes[i][1] - mi_y) < min_diam_y - delta and bottom_c < bt_counter:
      bottom.append(i) 
    elif abs(mesh.nodes[i][1] - ma_y) < min_diam_y - delta and top_c < bt_counter:
      top.append(i) 
    elif abs(mesh.nodes[i][2] - mi_z) < min_diam_z - delta and back_c < bf_counter:
      back.append(i)
    elif abs(mesh.nodes[i][2] - ma_z) < min_diam_z - delta and front_c < bf_counter:
      front.append(i)  
      
  #add boundary nodes    
  mesh.bc.append(('top', top))
  mesh.bc.append(('bottom', bottom))
  mesh.bc.append(('left', left))
  mesh.bc.append(('right', right))
  mesh.bc.append(('front', front))
  mesh.bc.append(('back', back))
  
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
  array = -1 * numpy.ones((nnx,nny,nnz))
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
          else:
            e.region = 'void2'
        else:
          # test if is in convex hull of non-design nodes
          if in_hull(center, hull):
            if not valid_position(center, coords,offset):
              e.region = 'void1'
            #elif type == "robot" and array[x,y,z] > 0.9:
            elif valid_ring_position(center, coords,offset): # add material ring around holes 
              e.region = 'non-design'  
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
  print 'mesh has ' + str(number) + "design and non-design elements"
  print 'volume = ' +str(float(number)/float(number + void3_count))
  return mesh

