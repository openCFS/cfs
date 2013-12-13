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
    print "region '" + region + "' not within regions " + str(regions.keys())

## give back elements with barycenters
# works 2D and 3D
# @return list barycenter tuple ordered by elements and min and max node coordinates and first region element dimensions
def centered_elements(hdf5_file, region):
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'] # for all regions
  reg_elements = hdf5_file['/Mesh/Regions/' + region + '/Elements']
  types = hdf5_file['/Mesh/Elements/Types']
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates']
  reg_nodes = hdf5_file['/Mesh/Regions/' + region + '/Nodes']
  
  # determine elem_dim from first region element dimensions
  elem_id = reg_elements[0] - 1
  node_coords = []
  for n in range(len(all_elements[elem_id])):
    node_coords.append(all_nodes[all_elements[elem_id][n]-1]) # numbers are one-based
  ma = numpy.array([max(node_coords,key=operator.itemgetter(0))[0], max(node_coords,key=operator.itemgetter(1))[1],  max(node_coords,key=operator.itemgetter(2))[2]])
  mi = numpy.array([min(node_coords,key=operator.itemgetter(0))[0], min(node_coords,key=operator.itemgetter(1))[1],  min(node_coords,key=operator.itemgetter(2))[2]])
  elem_dim = ma - mi
    
  # determine region dimensions, we need to resort for the desired region! Due to 1 to zero based conversion we need to do it manually :(
  nodes = numpy.zeros((len(reg_nodes), 3))
  for e in range(len(reg_nodes)):
    nodes[e] = all_nodes[reg_nodes[e] - 1]  
  min_dim = min(nodes[:,0]), min(nodes[:,1]), min(nodes[:,2])  
  max_dim = max(nodes[:,0]), max(nodes[:,1]), max(nodes[:,2])   
    
  result = []
  for e in range(len(reg_elements)):
    idx = reg_elements[e] - 1 # cfs writes one based
    nod = all_elements[idx]
    center = numpy.array([0.0, 0.0, 0.0])
    for n in range(len(nod)):
      center += all_nodes[nod[n]-1] # numbers are one-based
      # print "el=" + str(e) + " n=" + str(n) + " node=" + str(nod[n]) + "->" + str(nodes[nod[n]-1]) + " center=" + str(center) 
    center *= 1.0/len(nod)
    result.append(center)
    # print "e=" + str(e) + " idx=" + str(idx) + " nod=" + str(nod) + " center=" + str(center) 

  return result, min_dim, max_dim, elem_dim     
                
## find minimal and maximal coordinate
# @param coordinates as from centered_elements
def find_corners(centers):
  min = [1e30, 1e30, 1e30]
  max = [-1e30, -1e30, -1e30]
  
  if len(centers) == 1:
    min = [0.0,0.0,0.0]
    max = [1.0,1.0,0.0]
    
  
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
    print '  ' + name
    des = None
          
          
# returns a deep copied numpy array          
def get_element(hdf5_file, name, region, given_step=99999):
  ms = hdf5_file['/Results/Mesh/MultiStep_1']
  
  step = min((given_step, last_h5_step(hdf5_file)))
  key = "/Results/Mesh/MultiStep_1/Step_" + str(step) + "/" + name + "/" + region + "/Elements/Real"
  try:
    data = ms[key]
    return numpy.copy(data)
  except:
    raise Exception("cannot access '" + key + "' in " + str(hdf5_file.filename))
  

#f = h5py.File("/home/fwein/project/simp/hook.h5")
#s1 = get_element(f, "design_stiff1_plain", "mech")
#s2 = get_element(f, "design_stiff2_plain", "mech") 
#angle = get_element(f, "design_rotAngle_plain", "mech")
 
#t = to_mech_tensor(eval("[2.607069, 1.484705, 0.1626158, 0,0, 0.3030707]")) 
 
#f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")
# f = h5py.File("/home/fwein/tmp/l_sl-m_20-g_al-p_0.1-gamma_1e-07-tau_0.01.h5")
#r  = centered_elements(f)
#tensor = get_element(f, "mechTensor", "piezo")
# tensor = get_element(f, "mechTensor", "mech")
#rot = perform_rotations(tensor, 10)
# rot = perform_rotations(tensor, t21=True)
# orientational_stiffness(r, rot, 1500, 2.5).show()
