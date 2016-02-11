#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--y_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--z_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--width', help="width in m", type=float, default = 1.0)
parser.add_argument('--height', help="optional height in m", type=float, required = False)
parser.add_argument('--type', help="predefined mesh type", choices=['bulk2d', 'bulk3d', 'cantilever2d', 'cantilever2d_reinforced','lbm2d', 'lbm3d','msfem_two_load'], required = True)
parser.add_argument('--lbm', help="subtype for 'lbm'", choices=['two_inlet_one_outlet', 'pipe_bend','pipe','distributor','backstep','diffuser','two_inlet_two_outlet'])
parser.add_argument('--patch', help="define many regions", choices=['3x3', '4x4'])
parser.add_argument('--inclusion', help="inclusion for bulk2d and bulk3d", choices=["rect", "ball"])
parser.add_argument('--inclusion_size', help="mandatory size for inclusion as fraction of x-dimension (.9 is almost full)", type=float)
parser.add_argument('--file', help="optional give output file name. ")

args = parser.parse_args()

mesh_name = args.type

# sanity checks
if args.lbm and args.type <> "lbm":
  print "error: --lbm only for --type lbm"
  sys.exit()
  
if (args.inclusion or args.inclusion_size) and not (args.type == "bulk2d" or "cantilever2d" or "lbm3d"):
  print("inclusions currently only for --type {bulk2d,cantilever2d,lbm3d,lbm2d}")
  sys.exit()  
  
if (args.inclusion == 'ball' and (args.type == "lbm3d" or "lbm2d")):
  print("inclusion ball not implemented yet for lbm meshes")
  sys.exit()
  
if (args.inclusion and not args.inclusion_size) or (args.inclusion_size and not args.inclusion):
  print("inclusions require both --inclusion and --inclusion_size")
  sys.exit()  

if args.inclusion and args.patch:
  print "--inclusion and --patch don't go concurrently"
  sys.exit() 
  
if args.type == 'bulk3d':
  mesh = create_3d_mesh(args.type, args.res, args.y_res, args.z_res, args.inclusion, args.inclusion_size)  
elif args.type == 'bulk2d' or args.type.startswith('cantilever2d') or args.type.startswith('mbb') or args.type == 'msfem_test' or args.type == 'ghost' or args.type == 'triangle_msfem' or args.type == 'pressure2' or args.type == 'msfem_two_load':
  mesh = create_2d_mesh(args.type, args.res, args.y_res, args.width, args.height, args.inclusion, args.inclusion_size, args.patch)
elif args.type.startswith('lbm'):
  if args.lbm == None:
    print('error: --lbm subtype mandatory for --type lbm')
    sys.exit()
  if args.type == 'lbm2d':
     mesh = create_lbm2d(args.res, args.lbm,args.inclusion,args.inclusion_size)
  elif args.type == 'lbm3d' and args.lbm == 'backstep':
    mesh = create_backstep(args.res, args.y_res, args.z_res)
  elif args.type == 'lbm3d':
     mesh = create_lbm3d(args.res, args.y_res, args.z_res, args.lbm,args.inclusion,args.inclusion_size)
  else:
     mesh = create_3d_mesh(args.res, args.y_res, args.z_res)
  mesh_name = args.type +"_" + args.lbm
elif args.type == '3D':
  mesh = create_regular3d_mesh(args.type, args.res)
else:
  assert(False)  
  
res_name = '_' + str(args.res)
if (args.type == 'bulk2d' or args.type == 'bulk3d') and args.y_res <> None:
  res_name += '_' + str(args.y_res)
if args.type == 'bulk3d' and args.z_res:
  res_name += '_' + str(args.z_res)
if args.width <> 1.0:
  res_name += '-w_' + str(args.width).replace('.', '_')
if args.height is not None:
  res_name += '-h_' + str(args.height).replace('.', '_')
if args.inclusion:
  res_name += '_' + args.inclusion  + '_' + str(args.inclusion_size).replace('.', '_')
if args.patch:
  res_name += '_' + args.patch   

file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"
