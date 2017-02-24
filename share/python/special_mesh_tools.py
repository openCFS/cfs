#!/usr/bin/env python
import mesh_tool 
from hdf5_tools import num_nodes_by_type

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

def create_mesh_for_lufo_bracket(meshfile, all_nodes = [], elements = [], offset = 0., forces = [], force2 = [],supports = []):
  mesh = mesh_tool.Mesh()
  mesh.nodes = all_nodes
  # insert elements
  for i in range(len(elements)):
    e = mesh_tool.Element()
    e.nodes = (elements[i][1:]) 
    for k in range(len(e.nodes)):
      e.nodes[k] -= offset
    e.density = 1.
    if len(e.nodes) == 4:
      e.type = mesh_tool.TET4
      num_nodes = 4
    elif len(e.nodes) == 6:
      e.type = mesh_tool.WEDGE6
      num_nodes = 6
    elif len(e.nodes) == 8:
      e.type = mesh_tool.HEXA8
      num_nodes = 8
    elif len(e.nodes) == 3:
      e.type = mesh_tool.TRIA3
      num_nodes = 3
      e.region = 'surface'
    shell = 0
    for k in range(num_nodes):
        coord = mesh.nodes[e.nodes[k]]
        if (coord[1] - 37.49281)**2 + (coord[2] + 30.1963)**2 < 29.**2 and abs(coord[0]) < 5.:
          shell += 1
    if shell > 3:
      e.region = 'non-design'
    else:
      e.region = 'design'
    mesh.elements.append(e)
  #if force1 == []:
  #  force1 = [69692]
  #if force2 == []:
  #  force2 = [69693]
  #if support == []:
  #  support = [69659,69669,69670,69671,69672,69673,69674,69675,69677,69689]
  for f in forces:
    for force in f:
      force -= offset
  for supps in supports:
    for support in supps:
      support -= offset

  mesh.bc = []
  for i in range(len(forces)):
    mesh.bc.append(('force'+str(i+1), forces[i]))
  #mesh.bc.append(('force2', force2))
  for i in range(len(supports)):
    mesh.bc.append(('support'+str(i+1), supports[i]))
  mesh = mesh_tool.convert_to_sparse_mesh(mesh)    
  return mesh

## helper function: creates special mesh for apod6 and includes forces and supports
#@ meshfile name of the meshfile
#@ all_nodes nodes of the mesh can be unstructured
#@ elements elements of the mesh can be unstructured
#@ force. force nodes
#@ support. support nodes
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
  mesh = mesh_tool.Mesh()()  
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
      e = mesh_tool.Element()()
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
      e = mesh_tool.Element()()
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
      e = mesh_tool.Element()()
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
