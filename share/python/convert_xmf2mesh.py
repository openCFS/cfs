#!/usr/bin/env python
import h5py
import argparse
import numpy as np
import xml.etree.ElementTree as ET

from mesh_tool import *

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a xmf file")
parser.add_argument("--output", help="a mesh file")
args = parser.parse_args()

if not args.output:
  args.output = args.input[:-4] + '.mesh'

# parse xmf file
tree = ET.parse(args.input)  
root = tree.getroot()
for elem in root.findall('.//Grid'):
  topo = elem.find('Topology')
  if topo:
    elems_data = topo.find('DataItem').text
    filename_elems, elems_tag = [a.strip(' \r\n') for a in elems_data.split(':')]
  geom = elem.find('Geometry')
  if geom:
    nodes_data = geom.find('DataItem').text
    filename_nodes, nodes_tag = [a.strip(' \r\n') for a in nodes_data.split(':')]
  atts = elem.find('Attribute')
  if atts:
    attributes_data = atts.find('DataItem').text
    attributes = [a.strip(' \r\n') for a in attributes_data.split(':')]
    filename_attributes  = attributes.pop(0)
    
# read nodes
f = h5py.File(filename_nodes, 'r')
nodes = np.array(np.squeeze(f[nodes_tag[1:]]))

# read elements
f = h5py.File(filename_elems, 'r')
elements = np.array(np.squeeze(f[elems_tag[1:]]))

# read data
f = h5py.File(filename_attributes, 'r')
for a in attributes:
  if a[1:]=='ID':
    regions = np.array(np.squeeze(f[a[1:]]))

nnodes, dim = np.shape(elements)
nelem, elem_type = np.shape(elements)

# get mesh resolution
res = [-1]*dim
if elem_type==4:
  for res_d, ii in enumerate(res):
    ne = len(np.unique(nodes[:,ii]))-1
    if ne > 0:
      res_d = ne

# generate mesh
mesh = Mesh(*res)

mesh.nodes = nodes

for ii in range(nelem):
  e = Element()
  e.type = TRIA3 if elem_type==3 else QUAD4
  e.nodes = list(elements[ii])
  e.region = str(int(regions[ii]))
  mesh.elements.append(e)

write_ansys_mesh(mesh, args.output)
print('Output written to ' + args.output)