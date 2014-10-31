import sys
import h5py
import numpy
import operator
from mesh_tool import *


# open a hdf5 file as
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")

  
# checks if the region exists
# exits with an message if not
def validate_region(hdf5_file, region):
  regions = hdf5_file['/Mesh/Regions']
  if not any(k in regions.keys() for k in [region]):
    print "region '" + region + "' not within regions " + str(regions.keys())

def element_dimensions(elem_id, all_elements, all_nodes):
  node_coords = []
  for n in range(len(all_elements[elem_id])):
    node_coords.append(all_nodes[all_elements[elem_id][n] - 1])  # numbers are one-based
  ma = numpy.array([max(node_coords, key=operator.itemgetter(0))[0], max(node_coords, key=operator.itemgetter(1))[1], max(node_coords, key=operator.itemgetter(2))[2]])
  mi = numpy.array([min(node_coords, key=operator.itemgetter(0))[0], min(node_coords, key=operator.itemgetter(1))[1], min(node_coords, key=operator.itemgetter(2))[2]])
  elem_dim = ma - mi
  return elem_dim

def num_nodes_by_type(type_id):
  if type_id == 6:
    return 4
  if type_id == 16:
    return 6
  if type_id == 11:
    return 8
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
  assert(False)
    

# # give back elements with barycenters
# works 2D and 3D
# @return list barycenter tuple ordered by elements and min and max node coordinates and region element dimensions (first or all)
def centered_elements(hdf5_file, region, all_elem_dim=False, region_force=None, region_support=None):
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'].value  # for all regions
  reg_elements = hdf5_file['/Mesh/Regions/' + region + '/Elements'].value
  types = hdf5_file['/Mesh/Elements/Types'].value
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates'].value
  reg_nodes = hdf5_file['/Mesh/Regions/' + region + '/Nodes']
  if region_force <> None:
    reg_force_nodes = hdf5_file['/Mesh/Groups/' + region_force + '/Nodes']
  if region_support <> None:
    reg_support_nodes = hdf5_file['/Mesh/Groups/' + region_support + '/Nodes']
  
  # determine elem_dim from first region element dimensions or from all
  elem_dim = None
  if all_elem_dim:
    elem_dim = [0.0] * len(reg_elements)
    for i in range(len(reg_elements)):
      elem_dim[i] = element_dimensions(reg_elements[i] - 1, all_elements, all_nodes)
  else:
    elem_dim = element_dimensions(reg_elements[0] - 1, all_elements, all_nodes)        
    
  # determine region dimensions, we need to resort for the desired region! Due to 1 to zero based conversion we need to do it manually :(
  nodes = numpy.zeros((len(reg_nodes), 3))
  for e in range(len(reg_nodes)):
    nodes[e] = all_nodes[reg_nodes[e] - 1]
  # extract boundary force nodes from region_force if available
  if region_force <> None:
    nodes_force = numpy.zeros((len(reg_force_nodes), 3))
    for e in range(len(reg_force_nodes)):
      nodes_force[e] = all_nodes[reg_force_nodes[e] - 1]
  else:
    nodes_force = None
  
  # extract boundary support nodes from region_support if available
  if region_support <> None:
    # determine region dimensions, we need to resort for the desired region! Due to 1 to zero based conversion we need to do it manually :(
    nodes_support = numpy.zeros((len(reg_support_nodes), 3))
    for e in range(len(reg_support_nodes)):
      nodes_support[e] = all_nodes[reg_support_nodes[e] - 1]
  else:
    nodes_support = None
    
  min_dim = min(nodes[:, 0]), min(nodes[:, 1]), min(nodes[:, 2])  
  max_dim = max(nodes[:, 0]), max(nodes[:, 1]), max(nodes[:, 2])   
    
  result = []
  for e in range(len(reg_elements)):
    idx = reg_elements[e] - 1  # cfs writes one based
    nod = all_elements[idx]
    center = numpy.array([0.0, 0.0, 0.0])
    len_nod = num_nodes_by_type(types[idx])
    for n in range(len_nod):
      center += all_nodes[nod[n] - 1]  # numbers are one-based
      # print "el=" + str(e) + " n=" + str(n) + " node=" + str(nod[n]) + "->" + str(nodes[nod[n]-1]) + " center=" + str(center) 
    center *= 1.0 / len_nod
    result.append(center)
    # print "e=" + str(e) + " idx=" + str(idx) + " nod=" + str(nod) + " center=" + str(center)
  return result, min_dim, max_dim, elem_dim, nodes_force, nodes_support 
                
# # find minimal and maximal coordinate
# @param coordinates as from centered_elements
def find_corners(centers):
  min = [1e30, 1e30, 1e30]
  max = [-1e30, -1e30, -1e30]
  
  if len(centers) == 1:
    min = [0.0, 0.0, 0.0]
    max = [1.0, 1.0, 0.0]
    
  
  for e in range(len(centers)):
    test = centers[e]
    
    for c in range(3):
      if test[c] < min[c]:
        min[c] = test[c]
      if test[c] > max[c]:
        max[c] = test[c]
  return min, max      

def last_h5_step(hdf5_file):          
  ms = hdf5_file['/Results/Mesh/MultiStep_1']
  last = None
  for name in ms:
    if name.startswith('Step_'):
      last = name
      
  if last == None:
    raise Exception('no steps found in /Results/Mesh/MultiStep_1')     

  return int(last[5:])
          
def read_displacement(hdf5_file):
  u = hdf5_file['/Results/Mesh/MultiStep_1/Step_0/mechDisplacement/mech/Nodes/Real'].value
  return u

def read_density(hdf5_file):
  dens = hdf5_file['/Results/Mesh/MultiStep_1/Step_0/physicalPseudoDensity/mech/Elements/Real'].value
  return dens
    
# dumps meta data    
def dump_h5_meta(hdf5_file):   
  print 'Steps in "' + hdf5_file.filename + '":'
  ms = hdf5_file['/Results/Mesh/MultiStep_1']
  for name in ms:
    if name <> 'ResultDescription':
      print '  ' + name

  step = 'Step_' + str(last_h5_step(hdf5_file)) 
  ms = hdf5_file['/Results/Mesh/MultiStep_1/' + step]    
  print 'Results:'
  des = None
  for name in ms:
    print '  ' + name
    des = name
  if des == None:
    raise Exception("no design variables within last step " + name)
  
  ms = hdf5_file['/Results/Mesh/MultiStep_1/' + step + '/' + des]   
  print 'Regions (for ' + des + '):'
  for name in ms:
    size = len(hdf5_file['/Mesh/Regions/' + name + '/Elements'])
    print '  ' + name + ' with ' + str(size) + ' elements'
    
    des = None
          
# # Test for result
def has_element(hdf5_file, name, given_step=99999):
  try:   
    step = min((given_step, last_h5_step(hdf5_file)))
    ms = hdf5_file['/Results/Mesh/MultiStep_1/Step_' + str(step)]    
    for v in ms:
      if name == v:
        return True
  except Exception, e:
    print 'error probing for ' + name + ' in has_element: ', e

  return False
          
# returns a deep copied numpy array          
def get_element(hdf5_file, name, region, given_step=99999):
  step = min((given_step, last_h5_step(hdf5_file)))
  key = "/Results/Mesh/MultiStep_1/Step_" + str(step) + "/" + name + "/" + region + "/Elements/Real"
  try:
    return hdf5_file[key].value
  except:
    raise Exception("cannot access '" + key + "' in " + str(hdf5_file.filename))
# creates a mesh from hdf5 file  
def create_mesh_from_hdf5(hdf5_file, region, region_force=None, region_support=None):
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
  if region_support <> None:
    reg_support_nodes = hdf5_file['/Mesh/Groups/' + region_support + '/Nodes']
    mesh.bc.append((region_support, reg_support_nodes[:] - 1))

  
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
        e.region = region[0]
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
    # e.region = 'design'
    e.type = mesh_type_from_hdf5(types[i])
    mesh.elements.append(e) 
  return mesh
# f = h5py.File("/home/fwein/project/simp/hook.h5")
# s1 = get_element(f, "design_stiff1_plain", "mech")
# s2 = get_element(f, "design_stiff2_plain", "mech") 
# angle = get_element(f, "design_rotAngle_plain", "mech")
 
# t = to_mech_tensor(eval("[2.607069, 1.484705, 0.1626158, 0,0, 0.3030707]")) 
 
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")
# f = h5py.File("/home/fwein/tmp/l_sl-m_20-g_al-p_0.1-gamma_1e-07-tau_0.01.h5")
# r  = centered_elements(f)
# tensor = get_element(f, "mechTensor", "piezo")
# tensor = get_element(f, "mechTensor", "mech")
# rot = perform_rotations(tensor, 10)
# rot = perform_rotations(tensor, t21=True)
# orientational_stiffness(r, rot, 1500, 2.5).show()
