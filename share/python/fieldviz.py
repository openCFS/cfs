#!/usr/bin/env python
# this is a tool to visualize field data, similar to ParaView but for batch processing

import numpy as np
import argparse
import glob
import hdf5_tools as h5
import h5py
import scipy.interpolate
import os
import matplotlib.pyplot as plt


# give nodal data for a result over all regions
# list of coordinates and list of values
def get_nodal(cfs, result):
  c = []
  v = []
  
  full_coords = h5.get_coordinates(cfs)
  # all regions of our result
  regs = h5.get_result_regions(cfs, result)
  for r in regs:
    nodes = cfs['/Mesh/Regions/' + r + '/Nodes'][:] # 1-based , e.g. [211, 212, ...| 
    vals  = h5.get_node_result(cfs, args.field, r, args.step)
    assert len(nodes) == len(vals)
    for i in range(len(nodes)):
      #print(nodes[i],vals[i])
      t = full_coords[nodes[i]-1]
      c.append([t[0],t[1]])
      v.append(vals[i][0])

  return np.array(c),np.array(v)

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="openCFS .cfs files (hdf5")
parser.add_argument('--info', help="report meta data of file",action='store_true')
parser.add_argument('--field', help="nodal or element field to visualize")
parser.add_argument('--cmap', help="matplotlib colormap name", default='hot')
parser.add_argument('--step', help="data step, too large is last", type=int, default=99999)
parser.add_argument('--sequence', help="multi sequence step", type=int, default=1)
parser.add_argument('--save', help="optional filename to write image")
parser.add_argument('--saveall', help="don't show image but save all",action='store_true')
parser.add_argument('--noshow', help="suppress poping up an image",action='store_true')
args = parser.parse_args()

input = args.input if len(args.input) != 1 else glob.glob(args.input[0]) # for Windows potentially globalize 

if not args.info and not args.field:
  print('give either --info or --field for your', len(input),'files')
  os.sys.exit()

for file in input:
  problem = file[:-len('.cfs')]
  cfs = h5py.File(file,'r')
  
  try:
    res = h5.get_result_descriptions(cfs)

    if args.info:
      print('results and their regions for',file)
      max_name = max([len(r) for r in res])
      for r in res:
        print('  ',r.ljust(max_name),h5.get_result_regions(cfs, r))
        
    if args.field:
      if not args.field in res:
        print('error: given --field',args.field,'not in results',res,'for',file)
        os.sys.exit();
      
      # we assume a regualar mesh but cannot be sure, also we cannot be sure abotu the numerations
      # therfore we throw the data by their coordinates into an interpolation and sample from the interpolation
      # When the assumption holds, the data will stay original.
      
      coords, vals = get_nodal(cfs, args.field)
      
      X = np.unique(coords[:,0])
      Y = np.unique(coords[:,1])
      
      X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation
      interp = scipy.interpolate.LinearNDInterpolator(coords, vals) 
      V = interp(X, Y) # missing regions are filled by interpolation
     
      # use simple matplotlib to visualize 2d field data
      plt.pcolormesh(X, Y, V, shading='auto', cmap = args.cmap)
      plt.legend()
      plt.colorbar()
      plt.axis("equal")
  
      if args.saveall:
        plt.savefig(problem + '_' + args.field + '.png')
        plt.clf()
      elif not args.noshow:
        plt.show()
  except KeyError:
    print('h5 error reading',file)    
    
    
      
      
