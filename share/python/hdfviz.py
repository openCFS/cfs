#!/usr/bin/env python

# this utitility class helps to visualize cfs hdf5 result files

import h5py
import hdf5_tools as ht
import cfs_utils as cu
import numpy as np
import sys
import os
import matplotlib
import matplotlib.pyplot as plt
import argparse
import pathlib

# simple class holding an openCFS hdf5 result file and helping to visualize the data in a plane normal to x/y/z
# @param h5 h5py.File object or filename
# @param fieldname element name like 'vonMisesStress'
# @param regions list of regions like ['mech']
# @param plane eiter 'x', 'y' or 'z'
# @param res the resolution of the first axis of the plane
class HdfViz:
  
  def __init__(self, h5, fieldname, regions, plane, res = 1000, step = 0):

    if type(h5) == str:
      h5 = h5py.File(h5,'r')
    assert plane in ['x','y','z']
    assert type(regions) in [list,tuple]
      
    # the coordinates of the mesh nodes 0-based
    self.nodes = np.asarray(h5['/Mesh/Nodes/Coordinates'])

    # extremal coordinates of the nodes
    self.minimal = np.min(self.nodes,axis=0)
    self.maximal = np.max(self.nodes,axis=0)
    assert len(self.minimal) == 3 and len(self.maximal) == 3
    
    # the plane index 
    self.pd = 'xyz'.index(plane) # 'x' -> 0, 'z' -> 2
    self.d1 = 0 if self.pd in [1,2] else 1 # the coordinate index for the "x"-direction in the image
    self.d2 = 2 if self.pd in [0,1] else 1 # 'z'->1, 'y'->2, 'x'->2 for the image "y"-direction
    # print('plane',plane,'pd',self.pd,'d1',self.d1,'d2',self.d2) 
  
    # element connectivtiy
    self.elem_conn  = np.asarray(h5['/Mesh/Elements/Connectivity']) # first row: array([  1,   2,  13,  12, 122, 123, 134, 133], dtype=uint32)
    self.elem_types = np.asarray(h5['/Mesh/Elements/Types']) # [11, 11, 11, 11, 11, 11, 11, 11, ...]

    # read and keep elements and values by region idx
    self.elements = []
    self.values = []
    for reg in regions:
      self.elements.append(np.asarray(h5['/Mesh/Regions/' + reg + '/Elements']))
      self.values.append(np.asarray(h5['/Results/Mesh/MultiStep_1/Step_' + str(step) + '/' +  fieldname + '/' + reg + '/Elements/Real']))
    min_val = min([np.min(v) for v in self.values])
    max_val = max([np.max(v) for v in self.values])
    #print('values',fieldname, 'in',[min_val, max_val])
    self.val_norm = matplotlib.colors.Normalize(vmin=min_val, vmax=max_val,clip=True) # val_norm(x) in 0..1

    self.fieldname = fieldname
    self.step = step
    self.regions = regions
    self.h5 = h5 # cooler if we copy
    self.res = res
    
  # helper which gives for a 1-based elemnr the coordinates of the nodes
  def element(self, elemnr):
    coords = []
    for n in self.elem_conn[elemnr-1]: # 1-based nodenr
      coords.append(self.nodes[n-1])
    return coords  
  
  # helper which creates the matplotlib figure   
  def create_figure(self):
    d1, d2 = self.d1, self.d2
    max_x, max_y = self.maximal[d1], self.maximal[d2]
    min_x, min_y = self.minimal[d1], self.minimal[d2]
    dpi = self.res / ((max_x-min_x) * 100.0) 
    
    print('d1',d1,'d2',d2,'max_x',max_x,'max_y',max_y,'dpi',dpi)
      
    fig = plt.figure(figsize=(dpi*(max_x-min_x), dpi*(max_y-min_y)))
    ax = fig.add_subplot(111)
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    return fig, ax
  
    
  # we assume a very regular mesh, e.g. from create_mesh.py or mesh_tool.py
  def slice_regular_hexahedrons(self, offset):  
    pd = self.pd
    if offset  < self.minimal[pd] or offset > self.maximal[pd]:
      print('viz_plane: offset',offset,'out of range. Shall be within mesh bounds',[self.minimal[pd],self.maximal[pd]])
      sys.exit()
  
    fig, ax = self.create_figure()

    # cut the 1-based element with the offset plane and return the index of the matching nodes, might be empty
    def cut(coords):
      assert len(coords) == 8 # linear hexahedrons! Extend for higher order
      # sort for lower nodes and higer nodes
      lower = []
      higher = []
      for i, c in enumerate(coords): 
        if c[pd] <= offset:
          lower.append(i)
        else:
          higher.append(i)
      assert len(higher) in [0,4,8]
      #print('lower',lower,'higher',higher)
      return higher if len(higher) == 4 else []
    
    # faster version of cut() which simply returns if we are cut
    def within(coords):
      below = False
      above = False
      for c in coords:
        if below and above:
          return True
        if c[pd] <= offset:
          below = True
        else:
          above = True

      return False # we need half to be below

    # find the polygon index order for the four indices of an hexahedron.
    # shall be sufficient to do this for one hexahedron as it is expensive
    # the polygon needs to be in a proper order, which depends on the hexahedron numbering and the cutting plane
    def indices(higher, coords):
      assert len(higher) == 4
      d1, d2 = self.d1, self.d2
      # unsorted 2D coordinates
      u = np.array([[coords[i][d1],coords[i][d2]] for i in higher]) # [[0.1, 0.0], [0.0, 0.0], [0.1, 0.1], [0.0, 0.1]]
      #print(cu.nice(u))
      # barycenter
      bcx, bcy = [u[:,0].sum()/4, u[:,1].sum()/4]
      # sorted order for the higher indices
      s = []
      # we need the order (d1l, d2l), (d1h, d2l), (d1h, d2h), (d1l, d2h) for l = low and h = high
      s.append([i for i in range(4) if u[i,0] < bcx and u[i,1] < bcy][0])
      s.append([i for i in range(4) if u[i,0] > bcx and u[i,1] < bcy][0])
      s.append([i for i in range(4) if u[i,0] > bcx and u[i,1] > bcy][0])
      s.append([i for i in range(4) if u[i,0] < bcx and u[i,1] > bcy][0])
      return [higher[i] for i in s]
      
    # reuse indices with good hope for performance reasons
    idx = None

    for r in range(len(self.regions)):
      elem = self.elements[r]
      vals = self.values[r]
      for ei, e in enumerate(elem):
        # calulates a cut polygon if we cut
        coords = self.element(e)
        if within(coords):
          # assume that the element orientation is always the same, so reuse
          idx = indices(cut(coords), coords) if idx is None else idx
          #print(idx)
          x = [coords[i][self.d1] for i in idx]
          y = [coords[i][self.d2] for i in idx]
          #print('x',x,'y',y)
          ax.fill(x,y,color=matplotlib.cm.jet(self.val_norm(vals[ei])))

    return fig    

if __name__ == '__main__':       
  parser = argparse.ArgumentParser()
  parser.add_argument('input', help=".cfs input file")
  parser.add_argument('--field', help="fieldname to be visualized")
  parser.add_argument('--regions', nargs='+', help="list of regionnames")
  parser.add_argument('--plane', help="normal of cutting plane", choices=['x','y','z'], default='z')
  parser.add_argument('--offset', help='distance of cutting plane from origin', type=float, default=0.0)
  parser.add_argument('--res', help="first image dimension resolution", type=int, default=1000)
  parser.add_argument('--save', help="optionally save the output", action='store_true')
  parser.add_argument('--noshow', help="suppress popping up the image", action='store_true')
  args = parser.parse_args()
  if not os.path.exists(args.input):
    print("error: cannot open input file '" + args.input + "'")
    sys.exit()
  
  viz = HdfViz(args.input, args.field, args.regions, args.plane)
  fig = viz.slice_regular_hexahedrons(args.offset)
  #plt.show()
  if args.save:
    fn = pathlib.Path(args.input).stem 
    fn += '_' + args.field + '.png'
    fig.savefig(fn)
    print("create file '" + fn + "'") 
  if not args.noshow:
    plt.show(block=True)

