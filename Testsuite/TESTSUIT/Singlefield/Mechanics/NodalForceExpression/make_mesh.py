
from mesh_tool import *
import numpy as np

# this file is copied and modified from the PythonMesherVar test 


# distance of two coordinates given as 2D tutples
def dist(x1, x2):
  assert len(x1) == len(x2) == 2
  return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2) 

# cfs loads this file via embedded python and calls set_cfs_mesh()
# we expect as options x_res and inclusion_radius
def set_cfs_mesh(opt):
  x_res = int(opt['x_res'])
  rad = float(opt['inclusion_radius'])

  # generate a mesh of dimension 1x1 meters
  mesh = create_2d_mesh('bulk2d', x_res, x_res, 1.0)
  
  center = (.5, .5)

  # find inclusion elements
  for e in mesh.elements:
    coord = mesh.calc_barycenter(e)
    if dist(coord, center) <= rad:
      e.region = 'inclusion'
      # we could identify the bounary nodes also by checking the element's nodes coordinates
  
  # find nodes around inclusion
  # with a factor .5 instead of .5 we generate too much node points, which are within inclusion and not simulated 
  h = .3 * np.sqrt(2) * 1./x_res
  incl_nodes = []
  # also set outder boundary
  outer_nodes = []
  for i, node in enumerate(mesh.nodes):
    if dist(node, center) < rad + h and dist(node, center) > rad - h:
      #print(node)
      incl_nodes.append(i)
    x, y = node  
    if x == 0 or x == 1 or y == 0 or y == 1: # coorrdinates are rounded and can be compared
      outer_nodes.append(i) 
   
  mesh.bc.append(('inclusion_nodes',incl_nodes))
  mesh.bc.append(('outer_nodes',outer_nodes))

  return mesh

dic = {}
dic['x_res'] = 15
dic['inclusion_radius'] = .3
mesh = set_cfs_mesh(dic)

filename = 'inclusion_' + str(dic['x_res']) + '_' + str(dic['inclusion_radius']) + '.mesh'

write_ansys_mesh(mesh, filename)
print('wrote', filename)
