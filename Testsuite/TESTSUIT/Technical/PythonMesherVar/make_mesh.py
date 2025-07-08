# this script will be executed by embedded python from within cfs

from mesh_tool import *
import numpy as np

print('make_mesh.py from PythonMesherInclusion test')

# distance of two coordinates given as 2D tutples
def dist(x1, x2):
  assert len(x1) == len(x2) == 2
  return np.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2) 

# cfs loads this file via embedded python and calls set_cfs_mesh()
# we expect as options x_res and inclusion_radius
def set_cfs_mesh(opt):
  x_res = int(opt['x_res'])
  rad = float(opt['inclusion_radius'])
  import sys
  #print(sys.path)
  # generate a mesh of dimension 1x1 meters
  mesh = create_2d_mesh(x_res, x_res, 1.0)
  
  center = (.5, .5)

  # find inclusion elements
  for e in mesh.elements:
    coord = mesh.calc_barycenter(e)
    if dist(coord, center) <= rad:
      e.region = 'inclusion'
      # we could identify the bounary nodes also by checking the element's nodes coordinates
  
  # find nodes around inclusion
  h = .5 * np.sqrt(2) * 1./x_res
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

  # we can also write a ansys ascii mesh file
  #write_ansys_mesh(mesh, 'killme.mesh')
  
  # as we run within embedded python from cfs, call_cfs_input_reader in mesh_tool can import cfs and send the mesh to cfs
  # one could call theses functions directly without the mesh class. See call_cfs_input_reader() for the stuff to call
  call_cfs_input_reader(mesh)


   
