# this is a collection of legacy copy, most of it was copied from mesh_tool.py January 2022 on bigger refactoring effort
# Some code might be of interest, most will never be used.
# From Summer 2023 on this file shall be deleted without further notice! 

#!/usr/bin/env python
import mesh_tool 
try:
  from hdf5_tools import num_nodes_by_type
except:
  print("failed to import hdf5_tools in special_mesh_tools.py, hopefully we don't need it")    

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
    e.nodes = (elements[i][2:]) 
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
    #shell = 0
    #for k in range(num_nodes):
    #    coord = mesh.nodes[e.nodes[k]]
    #    if (coord[1] - 37.49281)**2 + (coord[2] + 30.1963)**2 < 29.**2 and abs(coord[0]) < 5.:
    #      shell += 1
    #if shell > 3:
    if elements[i][1] == 3:
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
      e.nodes = (elements[i][2:])
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
    write_ansys_mesh(mesh, meshfile+".mesh")     
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
  for i, n in enumerate(mesh.nodes):
    out.write('GRID%12d%8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
    #out.write('GRID    ' + '%-8d%-8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
  # Hexaeder elements
  for i, e in enumerate(mesh.elements):
    n = e.nodes
    if e.type == Ansys.HEXA8 and e.region == 'design':
      out.write('CHEXA%11d%8d%8d%8d%8d%8d%8d%8d+\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      out.write('+       %8d%8d\n'%(n[6]+1,n[7]+1))
      #out.write('CHEXA   ' + '%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+E%-6d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1,i+1))
      #out.write('+E%-6d%-8d%-8d\n'%(i+1,n[6]+1,n[7]+1))
    elif e.type == Ansys.HEXA8 and e.region == 'non-design':
      out.write('CHEXA%11d%8d%8d%8d%8d%8d%8d%8d+\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      out.write('+       %8d%8d\n'%(n[6]+1,n[7]+1))
      #out.write('CHEXA   ' + '%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d+E%-6d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1,i+1))
      #out.write('+E%-6d%-8d%-8d\n'%(i+1,n[6]+1,n[7]+1))
  # Wedge elements
  for i, e in enumerate(mesh.elements):
    n = e.nodes
    if e.type == Ansys.WEDGE6 and e.region == 'design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == Ansys.WEDGE6 and e.region == 'non-design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
  # write forces1
  for n in mesh.bc[0][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(1,n+1,0,5000./len(mesh.bc[0][1])))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(1,mesh.bc[0][1][i]+1,0,5000./len(mesh.bc[0][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

    # write forces2
  for n in mesh.bc[1][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,n+1,0,5000./len(mesh.bc[1][1])))

      # write forces3
  for n in mesh.bc[2][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,n+1,0,5000./len(mesh.bc[2][1])))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(2,mesh.bc[1][1][i]+1,0,5000./len(mesh.bc[1][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

  for n in mesh.bc[3][1]:
    out.write('SPC%13d%8d  13     \n'%(1,n+1))
  for n in mesh.bc[4][1]:
    out.write('SPC%13d%8d  2     \n'%(1,n+1))
  for n in mesh.bc[5][1]:
    out.write('SPC%13d%8d  2     \n'%(1,n+1))
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
  for n in mesh.nodes:
    out.write('GRID%12d        '% (i+1) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
    #out.write('GRID    ' + '%-8d%-8d'% (i+1,0) + str(n[0])[0:8] + str(n[1])[0:8] + str(n[2])[0:8] +'\n')
  # Hexaeder elements
  out.write('$$\n')
  out.write('$$  SPOINT Data\n')
  out.write('$$\n')
  out.write('$\n')
  out.write('$  CHEXA Elements: First Order\n')
  out.write('$\n')
  for e in mesh.elements:
    n = e.nodes
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
  for e in mesh.elements:
    n = e.nodes
    if e.type == WEDGE6 and e.region == 'design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,1,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == WEDGE6 and e.region == 'non-design':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
    elif e.type == WEDGE6 and e.region == 'non-design2':
      out.write('CPENTA%10d%8d%8d%8d%8d%8d%8d%8d\n'%(i+1,3,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
      #out.write('CPENTA  ' +'%-8d%-8d%-8d%-8d%-8d%-8d%-8d%-8d\n'%(i+1,2,n[0]+1,n[1]+1,n[2]+1,n[3]+1,n[4]+1,n[5]+1))
  # write forces1
  for bc in mesh.bc[0][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(1,bc[i]+1,0,5000./len(bc)))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(1,mesh.bc[0][1][i]+1,0,5000./len(mesh.bc[0][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

  # write forces2
  for bc in mesh.bc[1][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,bc[i]+1,0,5000./len(bc)))

  # write forces3
  for bc in mesh.bc[2][1]:
    out.write('FORCE%11d%8d%8d1.0     0.0     %-8f0.0\n'%(2,bc[i]+1,0,5000./len(bc)))
    #out.write('FORCE   ' + '%-8d%-8d%-8d%-8f'%(2,mesh.bc[1][1][i]+1,0,5000./len(mesh.bc[1][1])) + '%-8f%-8f%-8f'%(0.,1.,0.) + '\n')

  for bc in mesh.bc[3][1]:
    out.write('SPC%13d%8d  13     \n'%(1,bc[i]+1))
  for bc in mesh.bc[4][1]:
    out.write('SPC%13d%8d  2     \n'%(1,bc[i]+1))
  for bc in mesh.bc[5][1]:
    out.write('SPC%13d%8d  2     \n'%(1,bc[i]+1))
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

  for i in inp:
    if i[0:8].strip() == 'GRID':
      last_node_id = int(i[8:16].strip())

  nodes = [None] * (last_node_id+1)

  #rewind file
  file.seek(0)

  nodes_last_idx = 0 # count current last index of list 'nodes'

  for i in inp:
    #item = str.split(inp[i])
    #if i < len(inp)-1:
    #  item_n = str.split(inp[i+1])

    if len(i) == 0:
      continue
    # read and check header
    if i[0:8].strip() == 'GRID':
      if offset == -1: # set offset automatically
        offset = int(i[8:16].strip())
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
      elem.append([int(inp[i][8:16].strip()),int(inp[i][16:24].strip()), int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip())])
    elif inp[i][0:8].strip() == 'CPENTA':
      # read 3D wedge elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][16:24].strip()), int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip()),int(inp[i][56:64].strip()),int(inp[i][64:72].strip())])
    elif inp[i][0:8].strip() == 'CHEXA':
      # read 3D hexahedron elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][16:24].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip()),int(inp[i][48:56].strip()),int(inp[i][56:64].strip()),int(inp[i][64:72].strip()),int(inp[i][8:16].strip()),int(inp[i][16:24].strip())])
      i += 1
    elif inp[i][0:8].strip() == 'CTRIA3':
      # read 2D triangle elements
      elem.append([int(inp[i][8:16].strip()),int(inp[i][16:24].strip()),int(inp[i][24:32].strip()),int(inp[i][32:40].strip()),int(inp[i][40:48].strip())])
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
#   write_ansys_mesh(mesh, meshfile+".mesh",scale) # moved to create_mesh.py

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
    barycenter = mesh.calc_barycenter(e)
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
  # write_ansys_mesh(meshNew, filename[:-4]+"_voxelized_res_" + str(res) + ".mesh")

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



def add_bc_for_box_varel(mesh,bounds,pfem=None):
  print("bounds:",bounds)
  load = []
  support = []
  nodes = mesh.nodes
  eps = 1e-4
  if not pfem:
    for i, node in enumerate(nodes):
      # load on top panel y=1
      if numpy.isclose(node[1],1.0,1e-5):
      #if numpy.isclose(nodes[i][1],0.0):
      #if nodes[i][1] < 0.0001:
        load.append(i)
        continue
      # support edges
      #if numpy.isclose(nodes[i][1],0.0,1e-3):
      if node[1] < 1e-3:
      #if numpy.isclose(nodes[i][1],1.0):
      #if nodes[i][1] > 0.9999:
        if node[0] <= 1/20 + eps:
          support.append(i)
        elif node[2] <= 1/20 + eps:
          support.append(i)

    print("load:",len(load))
    print("support:",len(support))

    mesh.bc.append(("load", load))
    mesh.bc.append(("support", support))
  else:
    assert(pfem)
    for e in mesh.elements:
      bac = mesh.calc_barycenter(e)
      # define support surfaces
      if bac[1] < 1e-3:
        for b in [bac[0],bac[2]]:
          if b <= 1/20 + eps:
            surf = Element()
            surf.type = e.type
            surf.density = 1
            surf.nodes = e.nodes
            surf.region = "support_surf"
            mesh.elements.append(surf)

    print("added support surface elements for box varel")
    mesh = add_surf_elems_on_bb_faces(mesh, bounds)
    print("added support surface elements for bounding box faces")

  return mesh


def create_mesh_from_gmsh_special(meshfile,type):
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
  write_ansys_mesh(mesh, meshfile+".mesh")




def create_validation_mesh(coords,nondes_coords, s1, s2, s3, ip_nx, grad, dir, scale, d_f, valid_position, valid_ring_position, type = "apod6", thres = 0.0, csize = None, simp = None):
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
  if simp is False:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1, s2, s3, dx, dy, dz)
  else:
    ip_data, ip_near, out, ndim,scale_ = get_interpolation(coords, grad, s1, None, None, dx, dy, dz)
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
        if simp is False:
          s1, s2, s3 = ip_data[count][0:3]
        else:
          s1 = ip_data[count][0]
        l = [i,j,k]
        u = [i+dx_f,j+dy_f,k+dz_f]
        # if s1 < 0 point is out of the convex hull
        if s1 > 0.0 and simp is False:
          # 2scale optimization
          if s1 >= thres or s2 >= thres or s3 >= thres:
              create_cross_3D(array,l,u,s1,s2,s3,hole,res)
          else:
            create_cross_3D(array,l,u,void,void,void,hole,res)
        elif simp is False:
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
        for n in e.nodes:
          node = mesh.nodes[n]
          # test if element is above or below design region, bounds are given by non-design region
          if (node[1] >= -0.33403 - ((0.5*dy)/dy_f) and node[1] < -0.33303 + ((0.5*dy)/dy_f)) or (node[1] <= -0.35203 + ((0.5*dy)/dy_f) and node[1] > -0.35303 - ((0.5*dy)/dy_f)):
            count +=1
        # calculate center of element
        center = mesh.calc_barycenter(e)

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


def create_validation_mesh_for_pp_box(stlName,diffName,unionName):
  create_volume_mesh_with_gmsh(stlName)
  
  
# need this for pfem
# loop over all mesh elements and for each bb face, add respective surface elements
# for surface boundary conditions
def add_surf_elems_on_bb_faces(mesh,bounds):
  elems = mesh.elements()
  xmin,ymin,zmin,xmax,ymax,zmax = bounds

  for e in elems:
    baryc = mesh.calc_barycenter(e)
    surf = Element()
    reg = None
#     if numpy.isclose(baryc[0],xmin,1e-4):
#       reg = "left"
#     elif numpy.isclose(baryc[0],xmax,1e-4):
#       reg = "right"
#     elif numpy.isclose(baryc[1],ymin,1e-4):
#       reg = "front"
#     elif numpy.isclose(baryc[1],ymax,1e-4):
#       reg = "back"
#     elif numpy.isclose(baryc[2],zmin,1e-4):
#       reg = "bottom"
#     elif numpy.isclose(baryc[2],zmax,1e-4):
#       reg = "top"

    if numpy.isclose(baryc[1],ymax,1e-4):
      reg = "back"
    else:
      continue
#     # elem e not on bb face
#     if reg is None:
#       continue
#     else:
    surf.type = e.type
    surf.nodes = e.nodes
    mesh.elements.append(surf)


  return mesh

  
  

def create_gmsh_from_cfs_hdf5(hdf5_file, region, bcregions,output):
  # force names and support name has to be set manually, default force1, force2, force3, support, support2, support3
  mesh = create_mesh_from_hdf5(hdf5_file, region, bcregions)
  write_ansys_mesh(mesh, "test.mesh")
  out = open(output, "w")
  # gmsh header
  out.write('$MeshFormat \n')
  out.write('2.2 0 8\n')
  out.write('$EndMeshFormat \n')
  out.write('$Nodes \n')
  out.write(str(len(mesh.nodes))+' \n')
  dim = len(mesh.nodes[0])
  #write nodes
  for i, node in enumerate(mesh.nodes):  # write one based!
    out.write(str(i + 1) + "  " + str(node[0]) + "  " + str(node[1]))
    if dim == 3:
      out.write("  " + str(node[2]) + "\n")
    else:
      out.write("  0.0\n")
  #write elements
  out.write('$EndNodes \n')
  out.write('$Elements \n')
  out.write(str(len(mesh.elements)+len(mesh.bc[0][1]) + len(mesh.bc[1][1]) + len(mesh.bc[2][1]))+ '\n') #+ len(mesh.bc[3][1]) + len(mesh.bc[4][1]) + len(mesh.bc[5][1]) + len(mesh.bc[6][1]))+ '\n')
  # 1D boundary elements support, forces
  count = 0
  for bc in mesh.bc:
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
    for node in bc[1]:
      out.write(str(count+1) + ' ' +str(15) + ' 2 0 ' + str(id) + ' ' + str(node + 1)+' \n')
      count +=1
  # write 3D elements
  for e in mesh.elements:  # write one based!
    nodes = len(e.nodes)
    out.write(str(count + 1) + ' ' + str(5 if node_by_ansys_type(e.type) == 8 else 6) + ' 2 0 ' + str(2 if e.region == 'design' else 3))
    count +=1
    for n in range(node_by_ansys_type(e.type)):
      out.write(' '+ str(e.nodes[n] + 1))
    out.write('\n')

  out.write('$EndElements \n')
  #write forces, support

  out.write(' ')
  out.close()


# obviously this creates a mesh from density array
# @param input_array can be one of three cases: grayscale imgae (array of ints), color image (array of tuples (r,g,b(,a)), numpy.ndarray
def create_dense_mesh(input_array, nx, ny, mesh, threshold, scale, rhomin, multi_design=1, shearAngle=0, type=1, color_mode="random"):
  # convert angle to rad and check for feasibility
  angle = shearAngle / 180 * math.pi
  if (abs(angle) > math.pi / 2 - 1e-6):
    print('angle has to be between -pi/2 + 1e-6 and pi/2 - 1e-6')
    return 0
  dx = scale / nx / math.cos(angle)
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
        mesh.nodes.append((round(x * dx,14), round(y * dy,14)))
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
          if input_array[x, y][0] > 0 and input_array[x, y][1] == 0 and input_array[x, y][2] == 0:
            e.region = 'red'
            colorful_count += 1
          elif input_array[x, y][0] == 0 and input_array[x, y][1] > 0 and input_array[x, y][2] == 0:
            e.region = 'green'
            colorful_count += 1
          elif input_array[x, y][0] == 0 and input_array[x, y][1] == 0 and input_array[x, y][2] > 0:
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
  mesh.bc.append(("south", list(range(0, nx + 1))))
  mesh.bc.append(("north", list(range((nx + 1) * ny, (nx + 1) * (ny + 1)))))
  mesh.bc.append(("west", list(range(0, (nx + 1) * ny + 1, nx + 1))))
  mesh.bc.append(("east", list(range(nx, (nx + 1) * (ny + 1), nx + 1))))
  mesh.bc.append(("south_west", [0]))
  mesh.bc.append(("south_east", [nx]))
  mesh.bc.append(("north_west", [(nx + 1) * ny]))
  mesh.bc.append(("north_east", [(nx + 1) * (ny + 1) - 1]))



  # print mesh.bc
  msg = "dense resolution: " + str(nx) + " x " + str(ny) + " elements (" + str(scale) + "m x " + str(float(ny) / nx * scale) + "m)"
  msg += " -> " + str(mech_count) + " mech elements out of " + str(nx * ny) + " (" + str(float(mech_count) / (nx * ny) * 100.0) + " %)"
  msg += " with threshold " + str(threshold)
  if colorful_count > 0:
    msg += 'plus ' + str(colorful_count) + ' non gray elements'
  print(msg)  
  
  

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



  
