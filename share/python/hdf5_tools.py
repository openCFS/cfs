import sys
import h5py
import numpy
import operator


# open a hdf5 file as
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")

# u is an array of dim 2 with x as fast variable (left lower to right upper)
# hdf5_file = h5py.File(base + '.cfs', 'r')
# u = read_displacement(hdf5_file)
def cfs_displacement_to_vtk(u, nx, ny, name):
  from pyevtk.hl import gridToVTK
  
  assert len(u)==(nx+1)*(ny+1)
  assert u.shape[1] == 2
  x=np.linspace(0,1,(nx+1))      
  y=np.linspace(0,1,(ny+1))
  
  ux = np.zeros((nx+1,ny+1,1))
  uy = np.zeros((nx+1,ny+1,1))
   # fake third component of displacement vector
  uz = np.zeros((nx+1,ny+1,1))
  
  for j in range(0,ny+1):
    for i in range(nx, -1, -1):
      idx = j * (ny+1) + i
      print(i,j,idx,u[idx][0])
      ux[i,j,0] = u[idx][0]
      uy[i,j,0] = u[idx][1]

  print('x', x)      
  print('y', y)
  print('ux', ux)
  print('uy', uy)
  # gridToVTK expects 3D data, thus we fake the third dimension
  gridToVTK(name, x,y,np.zeros(1),pointData={"displacement":(ux,uy,uz)})
  print("# wrote '" + name + ".vtr'") 

  
# checks if the region exists
# exits with an message if not
def validate_region(hdf5_file, region):
  regions = hdf5_file['/Mesh/Regions']
  if not any(k in list(regions.keys()) for k in [region]):
    print("region '" + region + "' not within regions " + ", ".join(regions.keys()))

def element_dimensions(elem_id, all_elements, all_nodes):
  node_coords = []
  for n in range(len(all_elements[elem_id])):
    node_coords.append(all_nodes[all_elements[elem_id][n] - 1])  # numbers are one-based
#   print("all_nodes:",all_nodes)
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
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'][()]  # for all regions
  reg_elements = hdf5_file['/Mesh/Regions/' + region + '/Elements'][()]
  
  # check if reg_elements is list of list and flatten if necessary
  if any(isinstance(el, numpy.ndarray) for el in reg_elements):
    reg_elements = [el[0] for el in reg_elements]
  
  types = hdf5_file['/Mesh/Elements/Types'][()]
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates'][()]
  reg_nodes = hdf5_file['/Mesh/Regions/' + region + '/Nodes']
  if region_force != None:
    reg_force_nodes = hdf5_file['/Mesh/Groups/' + region_force + '/Nodes']
  if region_support != None:
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
  if region_force != None:
    nodes_force = numpy.zeros((len(reg_force_nodes), 3))
    for e in range(len(reg_force_nodes)):
      nodes_force[e] = all_nodes[reg_force_nodes[e] - 1]
  else:
    nodes_force = None
  
  # extract boundary support nodes from region_support if available
  if region_support != None:
    # determine region dimensions, we need to resort for the desired region! Due to 1 to zero based conversion we need to do it manually :(
    nodes_support = numpy.zeros((len(reg_support_nodes), 3))
    for e in range(len(reg_support_nodes)):
      nodes_support[e] = all_nodes[reg_support_nodes[e] - 1]
  else:
    nodes_support = None
    
  min_dim = [min(nodes[:, 0]), min(nodes[:, 1]), min(nodes[:, 2])]  
  max_dim = [max(nodes[:, 0]), max(nodes[:, 1]), max(nodes[:, 2])] 
  
  # element vertices described by node coords
  elements = []  
  # element vertices described by node ids
  elems_connectivity = []
  # coords all vertices in given region
  nodes_in_region = []
    
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
  
  # extract elements - each element is defined by its vertices
  for e in range(len(reg_elements)):
    elem = []
    idx = reg_elements[e] - 1
    nod = all_elements[idx]
    len_nod = num_nodes_by_type(types[idx])
    for n in range(len_nod):
      elem.append(all_nodes[nod[n] - 1])
      
    elements.append(elem)  
    elems_connectivity.append(nod)
    
  return result, min_dim, max_dim, elem_dim, nodes_force, nodes_support, elements, elems_connectivity, nodes, reg_nodes
                
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

def last_h5_step(hdf5_file,multistep=1):          
  ms = hdf5_file['/Results/Mesh/MultiStep_%i'%multistep]
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

# hdf5_file = h5py.File(name, 'r')
def read_displacement(hdf5_file,region = 'mech'):
  u = hdf5_file['/Results/Mesh/MultiStep_1/Step_' + str(last_h5_step(hdf5_file)) + '/mechDisplacement/'+region+'/Nodes/Real'][()]
  return u

# dumps meta data    
def dump_h5_meta(hdf5_file):   
  print('Steps in "' + hdf5_file.filename + '":')
  ms = hdf5_file['/Results/Mesh/MultiStep_1']
  for name in ms:
    if name != 'ResultDescription':
      print('  ' + name)

  step = 'Step_' + str(last_h5_step(hdf5_file)) 
  ms = hdf5_file['/Results/Mesh/MultiStep_1/' + step]    
  print('Results:')
  des = None
  for name in ms:
    print('  ' + name)
    des = name
  if des == None:
    raise Exception("no design variables within last step " + name)
  
  ms = hdf5_file['/Results/Mesh/MultiStep_1/' + step + '/' + des]   
  print('Regions (for ' + des + '):')
  for name in ms:
    size = len(hdf5_file['/Mesh/Regions/' + name + '/Elements'])
    print('  ' + name + ' with ' + str(size) + ' elements')
    des = None
    
# # Test for result
def has_element(hdf5_file, name, given_step=99999):
  try:   
    step = min((given_step, last_h5_step(hdf5_file)))
    ms = hdf5_file['/Results/Mesh/MultiStep_1/Step_' + str(step)]    
    for v in ms:
      if name == v:
        return True
  except Exception as e:
    print('error probing for ' + name + ' in has_element: ', e)

  return False
          
# returns a deep copied numpy array          
def get_element(hdf5_file, name, region, given_step=99999):
  step = min((given_step, last_h5_step(hdf5_file)))
  key = "/Results/Mesh/MultiStep_1/Step_" + str(step) + "/" + name + "/" + region + "/Elements/Real"
  try:
    return hdf5_file[key][()]
  except:
    raise Exception("cannot access '" + key + "' in " + str(hdf5_file.filename))

def get_step_values(hdf5_file) :
    """
    return the step values as a list of arrays for each multi-sequence step

    Parameters
    ----------
    h5: h5py.File object or str

    Returns
    -------
    out : list of arrays

    Example
    -------
    A single step at frequency 1
    >>> get_step_values('./TESTSUIT/Singlefield/Mechanics/uniaxStrainPwavePML3d/uniaxStrainPwavePML3d.h5ref')
    [array([1.])]

    A single eigenfrequecy step
    >>> get_step_values(dampedEV2D)
    [array([ 254.66060483,  636.80184511,  681.93162371, 1133.77662285,
            1144.0149628 , 1378.19737112, 1668.47028788, 1931.85306237,
            1987.25030061, 1992.65687422, 2194.52472254, 2305.77247832])]

    Check how many steps there are
    >>> step_val_list = get_step_values(dampedEV2D)
    >>> len(step_val_list)
    1
    """
    if type(hdf5_file) == str :
        from h5py import File
        with File(hdf5_file,'r') as h5:
            return get_step_values(h5)
    else:
        from numpy import allclose
        try:
          multisteps = hdf5_file['Results/Mesh']
        except:
          raise Exception("no Mesh results in %s"%hdf5_file.filename)
        step_values = []
        for msk in multisteps.keys():
          result_keys = list(multisteps[msk]['ResultDescription'].keys())
          rev_vals = multisteps[msk]['ResultDescription'][result_keys[0]]['StepValues'][:]
          # check step vals
          for rk in result_keys[1:]:
            if not allclose(rev_vals,multisteps[msk]['ResultDescription'][rk]['StepValues'][:]) :
              Exception("step values are different for '" + rk + "' and '"+result_keys[0]+"' in " + msk)
            #.values(
          step_values.append(rev_vals)
        return step_values

# returns nodal or elemental results as numpy array
def get_result(hdf5_file,result,region=None,step='last',multistep=1) :
    """
    read data from a hdf5-file
    
    Parameters
    ----------
    hdf5_file : h5py.File or str
        CFS++ hdf5 data file
    result: string 
        specifies the results to return: e.g. 'accuPressure','mechDisplacement',...
    region: string 
        region name
    step: integer, list, or string  
        defining the step as single integer, list of integersor or 'last' for laststep or
        'all' for all steps
      
    Returns
    -------
    out : ndarray

    Example
    -------
    Extract real data
    >>> U = get_result(Plate3D,'mechDisplacement')
    >>> U[-6:,:]
    array([[ 0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
       	   [-2.01311228e-05,   2.01311228e-05,   2.80932210e-05],
           [-4.28143434e-04,  -4.28143434e-04,   2.07611520e-02],
           [ 4.28143434e-04,  -4.28143434e-04,   2.07611520e-02],
           [-4.28143434e-04,   4.28143434e-04,   2.07611520e-02],
           [ 4.28143434e-04,   4.28143434e-04,   2.07611520e-02]])
    """
    if type(hdf5_file) == str :
        from h5py import File
        with File(hdf5_file,'r') as h5:
            return get_result(h5,result,region,step,multistep)
    else:
        from numpy import array, squeeze
        h5_ms = hdf5_file['Results/Mesh/MultiStep_%i'%multistep] # extract multistep
        if step=='last':
            steps=[last_h5_step(hdf5_file,multistep)]
        elif step=='all' :
            steps=h5_ms['ResultDescription/%s/StepNumbers'%(result)][:]
        elif type(step)==int :
            steps=[step]
        elif hasattr(step, '__iter__') :
            steps=step
        res = []
        for step in steps:
            h5_s = h5_ms['Step_%i'%step] # extract step
            h5_res = h5_s[result] # extract result
            if region==None :
                if len(h5_res.keys())>1 :
                    raise Exception("No region specified but more than one region present for result '"+result+"'in '"+hdf5_file.filename+"', MultiStep_%i, Step_%i"%(multistep,step)+" Available regions: "+", ".join(h5_res.keys()))
                else :
                    region = [k for k in h5_res.keys()][0]
            h5_res_reg = h5_res[region] # extraxt region
            res_type = list(h5_res_reg.keys())[0] # read result type (Nodes or Elements)
            if 'Imag' in h5_res_reg[res_type].keys() :
                res.append( h5_res_reg[res_type]['Real'][:] + 1j*h5_res_reg[res_type]['Imag'][:] )
            else :
                res.append( h5_res_reg[res_type]['Real'][:] )
        return squeeze(array(res))

def get_all_results(hdf5_file,region=None,step='last',multistep=1):
  ms = hdf5_file['/Results/Mesh/MultiStep_%i/ResultDescription'%multistep]
  res = dict()
  for name in ms:
    design = ms[name]
    if name == None:
      continue
    res[name] = get_result(hdf5_file,name,region,step,multistep)

  return res

def get_subregion_idx(hdf5_file,region,subregion,rtype='Nodes') :
    """
    returns indices of the elements in 'subregion' with respect to the indices in 'region'

    Parameters
    ----------
    hdf5_file : h5py.File or str
        CFS++ hdf5 data file
    region : string
        region name
    sugregion : string
        subregion name
    rtype : string (optional)
        Either 'Nodes'(default) or 'Elements', defines the region type

    Returns
    -------
    out : ndarray

    Example
    -------
    Extract data for subregion 'side'
    >>> U = get_result(Plate3D,'mechDisplacement','plate')
    >>> I = get_subregion_idx(Plate3D,'plate','side')

    # Of course, all displacements on the clamped side mst be zero
    >>> U[I,:] # doctest: +ELLIPSIS
    array([[0.,  0.,  0.],
           [0.,  0.,  0.],
      ...
           [0.,  0.,  0.]])
    """
    if type(hdf5_file) == str :
        from h5py import File
        with File(hdf5_file,'r') as h5:
            return get_subregion_idx(h5,region,subregion,rtype)
    else:
        from numpy import array, argwhere
        Is = hdf5_file['Mesh']['Regions'][subregion][rtype][:]-1
        Ir = hdf5_file['Mesh']['Regions'][region][rtype][:]-1
        Isr = array([argwhere(Ir==i).ravel() for i in Is if i in Ir]).ravel()
        return Isr

def get_coordinates(hdf5_file,region=None) :
  """
  return nodal coordinates, optional for nodes in a certain region

  Parameters
  ----------
  hdf5_file : h5py.File or str
      CFS++ hdf5 data file
  region : string, optional
      region name for a subset of coordinates

  Returns
  -------
  out : ndarray

  Examples
  --------
  >>> X = get_coordinates(Plate3D)

  # Show coordinates of nodes 2, 4, 6 (index shift) -> X[[1,3,5],:]
  >>> X[[1,3,5],:]
  array([[12. ,   7.2,   0. ],
         [11.5,   4.8,   0. ],
         [ 0.5,   7.2,   0. ]])
  """
  if type(hdf5_file) == str :
    from h5py import File
    with File(hdf5_file,'r') as h5:
      return get_coordinates(h5,region)
  else:
    if not region==None :
      I = hdf5_file['Mesh/Regions/%s'%region]['Nodes'][:] - 1
      return hdf5_file['Mesh/Nodes/Coordinates'][:][I,:]
    else :
      return hdf5_file['Mesh/Nodes/Coordinates'][:]

def get_centroids(hdf5_file,region=None) :
    """
    return the centroid coordinates for elemnts (in a certain region)

    The centroid is computed as the arithmetic mean of all nodal coodinates.

    Parameters
    ----------
    hdf5_file : h5py.File or str
        CFS++ hdf5 data file
    region : string, optional
        region name for a subset of elements

    Returns
    -------
    out : ndarray
        first dimension sorted like elements in hdf5_file

    Examples
    --------
    Get centroids for a single region result file
    >>> C = get_centroids(Plate3D)
    >>> C[0,:] # first element
    array([11.75 , 6. , 0.175])
    >>> C[-1,:] # last element
    array([9.35, 9.35, 0.35])

    Reading directly from a file path
    >>> get_centroids('TESTSUIT/Singlefield/Mechanics/uniaxStrainPwavePML3d/uniaxStrainPwavePML3d.h5ref') # doctest: +ELLIPSIS
    array([[0.5 , 0.  , 0.05],
           [0.5 , 0.  , 0.15],
           ...
           [1.  , 1.  , 2.85],
           [1.  , 1.  , 2.95]])
    """
    if type(hdf5_file) == str :
        from h5py import File
        with File(hdf5_file,'r') as h5:
            return get_centroids(h5,region)
    else:
        from numpy import sum, unique, argwhere, zeros, arange, mean
        # get connectivity and nodal coordinates
        conn = hdf5_file['Mesh/Elements/Connectivity'][:]
        coord = hdf5_file['Mesh/Nodes/Coordinates'][:]
        # determine indices of region
        if not region==None:
            validate_region(hdf5_file,region)
            I = hdf5_file['Mesh/Regions/%s/Elements'%region][:] - 1
        else :
            I = arange(conn.shape[0])
        # allocate result
        center = zeros([len(I),coord.shape[1]])
        # compute the centriods
        etype = hdf5_file['Mesh/Elements/Types'][:][I] # element types
        for et in unique(etype) : # sum operations for unique element types
            It = argwhere(etype == et)[:,0] # index of the type-elements in region
            Nnodes = sum(conn[I[It[0]]]>0) # determine how many nodes the element-type has
            nids = conn[I[It],:Nnodes] # node ids, only take used columns
            center[It,:] = mean(coord[nids-1,:],axis=1) # compute center as arithmetic mean
        return center

if __name__ == "__main__":

    # for testing in the testsuite, this file is run with the TESTSUIT_DIR
    # as working directory
    # $ python hdf5_tools.py -v
    # Thus, you must cd to the correct location for interactive testing.

    #%% Load some files (paths relative to TESTSUIT_DIR)
    from h5py import File
    # should load from TESTSUITE_DIR
    Plate3D = File('TESTSUIT/Singlefield/Mechanics/Plate3D/Plate3D.h5ref','r')
    dampedEV2D = File('TESTSUIT/Singlefield/Mechanics/dampedEV2D/dampedEV2D.h5ref','r')

    # finally run doctest
    import doctest
    result = doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # and exit with the number of failures as exit code
    sys.exit(result.failed)
