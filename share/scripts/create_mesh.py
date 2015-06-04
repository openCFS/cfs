#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--res", help="x-discretization of length 1m", type=int, required=True)
parser.add_argument('--y_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required=False)
parser.add_argument('--z_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required=False)
parser.add_argument('--type', help="predefined mesh type", choices=['bulk2d', 'bulk3d', 'cantilever2d', 'cantilever2d_reinforced', 'lbm', 'msfem_test', 'ghost','triangle_msfem','pressure2','msfem_two_load'], required=True)

parser.add_argument('--lbm', help="subtype for 'lbm'", choices=['two_inlet_one_outlet', 'pipe_bend'])
parser.add_argument('--file', help="optional give output file name. ")

args = parser.parse_args()

mesh_name = args.type

if args.type == 'bulk3d':
  mesh = create_3d_mesh(args.res, args.y_res, args.z_res)
elif args.type == 'bulk2d' or args.type.startswith('cantilever2d') or args.type.startswith('mbb') or args.type == 'msfem_test' or args.type == 'ghost' or args.type == 'triangle_msfem' or args.type == 'pressure2' or args.type == 'msfem_two_load':
  mesh = create_2d_mesh(args.type, args.res, args.y_res)
elif args.type == 'lbm':
  if args.lbm == None:
    print('error: --lbm subtype mandatory for --type lbm')
    sys.exit()
  mesh = create_lbm(args.res, args.lbm)
  mesh_name = args.lbm
elif args.type == '3D':
  mesh = create_regular3d_mesh(args.type, args.res)
else:
  assert(False)  
  
res_name = '_' + str(args.res)
if (args.type == 'bulk2d' or args.type == 'bulk3d') and args.y_res <> None:
  res_name += '_' + str(args.y_res)
if args.type == 'bulk3d' and args.z_res:
  res_name += '_' + str(args.z_res)
  
file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"
