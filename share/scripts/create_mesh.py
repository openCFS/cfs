#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--y_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--z_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--type', help="predefined mesh type", choices=['bulk2d', 'bulk3d', 'cantilever2d', 'cantilever2d_reinforced','lbm2d', 'lbm3d'], required = True)
parser.add_argument('--lbm', help="subtype for 'lbm'", choices=['two_inlet_one_outlet', 'pipe_bend', 'extend_inlet','pipe','distributor'])
parser.add_argument('--file', help="optional give output file name. ")
parser.add_argument('--inclusion', help="inclusion for bulk2d and bulk3d", choices=["rect", "ball"])
parser.add_argument('--inclusion_size', help="mandatory size for inclusion as fraction of x-dimension (.9 is almost full)", type=float)

args = parser.parse_args()

mesh_name = args.type

if (args.inclusion or args.inclusion_size) and not (args.type == "bulk2d" or "cantilever2d"):
  print("inclusions currently only for --type buld2d")
  sys.exit()  
  
if (args.inclusion and not args.inclusion_size) or (args.inclusion_size and not args.inclusion):
  print("inclusions require both --inclusion and --inclusion_size")
  sys.exit()  

if args.type == 'bulk3d':
  mesh = create_3d_mesh(args.type, args.res, args.y_res, args.z_res, args.inclusion, args.inclusion_size)
elif args.type == 'bulk2d' or args.type.startswith('cantilever2d') or args.type.startswith('mbb'):
  mesh = create_2d_mesh(args.type, args.res, args.y_res)
elif args.type.startswith('lbm'):
  if args.lbm == None:
    print('error: --lbm subtype mandatory for --type lbm')
    sys.exit()
  if args.type == 'lbm2d':
     mesh = create_lbm2d(args.res, args.lbm)
  elif args.type == 'lbm3d':
     mesh = create_lbm3d(args.res, args.y_res, args.z_res, args.lbm)
  else:
     mesh = create_3d_mesh(args.res, args.y_res, args.z_res)
  mesh_name = args.type +"_" + args.lbm
else:
  assert(False)  
  
res_name = '_' + str(args.res)
if (args.type == 'bulk2d' or args.type == 'bulk3d') and args.y_res <> None:
  res_name  += '_' + str(args.y_res)
if args.type == 'bulk3d' and args.z_res:
  res_name  += '_' + str(args.z_res)
if args.inclusion:
  res_name += '_' + args.inclusion + '_' + str(args.inclusion_size).replace('.', '_') 
  
file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"
