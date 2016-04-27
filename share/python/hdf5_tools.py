import sys
import h5py
import numpy
import operator


# open a hdf5 file as
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")

  
# checks if the region exists
# exits with an message if not
def validate_region(hdf5_file, region):
  regions = hdf5_file['/Mesh/Regions']
  if not any(k in regions.keys() for k in [region]):
    print "region " + region + " not within regions " + str(regions.keys())

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
  if type_id == 9:
    return 4
  assert(False)

# # give back elements with barycenters
# works 2D and 3D
# @return list barycenter tuple ordered by elements and min and max node coordinates and region element dimensions (first or all)
def centered_elements(hdf5_file, region, all_elem_dim=False, region_force=None, region_support=None,centered = True):
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
  if centered:
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
  else:
    # append nodes to result instead of element centers
    for i in range(len(nodes[:,0])):
      result.append([nodes[i,0],nodes[i,1],nodes[i,2]])
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
      if last != None:
        thisInt = int(name.split('_')[-1])
        if thisInt > lastInt:
          last = name
          lastInt = thisInt
      else:
        last = name
        lastInt = int(name.split('_')[-1])                
      
  if last == None:
    raise Exception('no steps found in /Results/Mesh/MultiStep_1')     

  return int(last[5:])

def read_displacement(hdf5_file,nr,region = None):
  if region is None:
    u = hdf5_file['/Results/Mesh/MultiStep_1/Step_'+str(nr)+'/mechDisplacement/mech/Nodes/Real'].value
  else:
    # special case for debuggin apod6. Generalize if necessary
    f = h5py.File(hdf5_file)
    u = f['/Results/Mesh/MultiStep_1/Step_'+str(nr)+'/mechDisplacement/non-design/Nodes/Real'].value
    non_des = f['/Mesh/Regions/non-design/Nodes'].value
      
    force_nodes = f['/Mesh/Groups/'+str(region)+'/Nodes'].value
    u_max = -100000.
    u_average = 0.
    for i in range(len(force_nodes)):
      j = force_nodes[i]
      for k in range(len(non_des)):
        if non_des[k] == j:
          break
      if u_max < u[k][1]:
        u_max = u[k][1]
      u_average += u[k][1]
    u_average /= len(force_nodes)
    print 'u_max = ' + str(u_max)
    print 'u_average = '+ str(u_average)
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

