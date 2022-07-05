#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser(description="generate a basic mesh. bulk2d/3d is unit square/cube by default. Further default geometries are cantilever 3x2(x1) and mbb2d (1x0.5) ")
parser.add_argument("--res", help="elements in x-direction", type=int, required = True )
parser.add_argument('--y_res', help="elements in y-direction", type=int, required = False )
parser.add_argument('--z_res', help="elements in x-direction", type=int, required = False )
parser.add_argument('--width', help="width in m, default 1 m", type=float, default = 1.0)
parser.add_argument('--height', help="optional height in m", type=float, required = False)
parser.add_argument('--depth', help="optional depth in m", type=float, required = False)
parser.add_argument('--type', help="mesh type", choices=['bulk2d', 'cantilever2d', 'mbb2d', 'triangles', 'bulk3d', 'cantilever3d'], required = True)
parser.add_argument('--file', help="optional give output file name. ")
parser.add_argument('--pfem', help="sets additional boundary elements for b.c.", action='store_true', default=False)
parser.add_argument('--numbering', help="numbering of nodes and elements (only 2D for now)", choices=['row_major', 'col_major'], default='row_major')

args = parser.parse_args()

mesh= None 
    
if args.type in ['bulk2d', 'cantilever2d', 'mbb2d']:
  if args.type != 'bulk2d':
    if args.y_res is not None or args.width != 1.0 or args.height is not None:
      print('Error: when using predefined geometry, set only --res or use bulk2d')
      sys.exit(-1)
  
    if args.type == 'cantilever2d': 
      args.width  = 3.0
      args.height = 2.0 
        
    if args.type == 'mbb2d': 
      args.width  = 1.0 
      args.height = 0.5 
  # create_2d_mesh(x_res, y_res = None, width = 1.0, height = None, pfem=False, row_major=True, triangles = False):
  mesh = create_2d_mesh(args.res, args.y_res, args.width, args.height, pfem=args.pfem, row_major = args.numbering == 'row_major', triangles = args.type == 'triangles')
elif args.type in ['bulk3d', 'cantilver3d']:
  if args.type == 'cantilever3d':
    if args.y_res is not None or args.z_res is not None or args.height is not None or args.depth is not None:
      print('Error: when using predefined geometry, set only --res or use bulk3d')
      sys.exit(-1)
    args.width  = 3.0
    args.height = 2.0
    args.depth  = 1.0 
  mesh = create_3d_mesh(args.res, args.y_res, args.z_res, args.width, args.height, args.depth, pfem = args.pfem) # no col_major implemented
  
res_name = '_' + str(args.res)
if args.type in ['bulk2d', 'bulk3d']:
  if args.y_res is not None:
    res_name += '_' + str(args.y_res)
  if args.z_res is not None:
    res_name += '_' + str(args.z_res)
  if args.width != 1.0:
    res_name += '-w_' + str(args.width).replace('.', '_')
  if args.height is not None:
    res_name += '-h_' + str(args.height).replace('.', '_')
  if args.depth is not None:
    res_name += '-d_' + str(args.depth).replace('.', '_')
if args.pfem:
  res_name += "_pfem" 

file = args.type + res_name + '.mesh' if args.file == None else args.file 

write_ansys_mesh(mesh, file)
print("created file '" + file + "' with " + str(len(mesh.elements)) + " elements")

