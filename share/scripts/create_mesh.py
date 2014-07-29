#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="long side resolution of mesh if action is 'mesh", type=int, required = True )
parser.add_argument('--type', help="predefined mesh type", choices=['cantilever2d', 'cantilever2d_reinforced', 'mbb', 'mbb_reinforced', 'lbm'], required = True)
parser.add_argument('--lbm', help="subtype for 'lbm'", choices=['two_inlet_one_outlet', 'pipe_bend'])
parser.add_argument('--file', help="optional give output file name")

args = parser.parse_args()

mesh_name = args.type

if args.type == 'cantilver2d' or args.type == 'cantilver2d_reinforced':
  mesh = create_cantilever2d_mesh(args.type, args.res)
elif args.type == 'mbb' or args.type == 'mbb_reinforced':  
  mesh = create_mbb_mesh(args.type, args.res)
elif args.type == 'lbm':
  if args.lbm == None:
    print('error: --lbm subtype mandatory for --type lbm')
    sys.exit()
  mesh = create_lbm(args.res, args.lbm)
  mesh_name = args.lbm
else:
  assert(False)  
  
file = mesh_name + '_' + str(args.res) + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"
