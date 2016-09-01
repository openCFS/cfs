#!/usr/bin/env python
import argparse
from draw_profile_functions import *
from mesh_tool import *

# def create_mesh_with_profiles(x1, x2, y1, y2, z1, z2, xres, yres, zres,ipo):
def create_mesh_with_profiles(args):
  array = create_profiles_array(args)
  if args.z1 == 0.0 and args.z2 == 0.0:
    mesh = create_2d_mesh_from_array(array)
  else:
    mesh = create_3d_mesh_from_array(array)
  
#   visualize_structure(array,nx,ny,nz)
  
  return mesh

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--type', help="predefined mesh type", choices=['profiles2d', 'profiles3d'], required = True)
parser.add_argument('--file', help="optional give output file name. ")
# for profile functions
parser.add_argument('--x1', help="first radius for profile of bar in x-direction; 0 <= x1 <= 1", type=float, default=0.5)
parser.add_argument('--x2', help="second radius for profile of bar in x-direction; 0 <= x2 <= 1", type=float, default=0.5)
parser.add_argument('--y1', help="first radius for profile of bar in y-direction; 0 <= y1 <= 1", type=float, default=0.5)
parser.add_argument('--y2', help="second radius for profile of bar in y-direction; 0 <= y2 <= 1", type=float, default=0.5)
parser.add_argument('--z1', help="first radius for profile of bar in z-direction; 0 <= z1 <= 1", type=float, default=0.5)
parser.add_argument('--z2', help="second radius for profile of bar in z-direction; 0 <= z2 <= 1", type=float, default=0.5)
parser.add_argument('--profile', help="type of profile functions", choices=["linear","circular"])
parser.add_argument('--skip_x', help="don't show bar in x direction", action='store_true')
parser.add_argument('--skip_y', help="don't show bar in y direction", action='store_true')
parser.add_argument('--skip_z', help="don't show bar in z direction", action='store_true')

args = parser.parse_args()

mesh_name = args.type

# sanity checks
if args.type == "profiles2d" and not (args.x1 and args.x2 and args.y1 and args.y2) and (args.z1 or args.z2):
  print("error: profiles2d needs values for x1,x2 and y1,y2")
  sys.exit()

if args.type == "profiles3d" and not (args.x1 and args.x2 and args.y1 and args.y2 and args.z1 and args.z2):
  print("error: profiles3d needs values for x1,x2 and y1,y2 and z1,z2")
  sys.exit()
  
if args.type == 'profiles2d':
  mesh = create_mesh_with_profiles(args)
elif args.type == 'profiles3d':  
  mesh = create_mesh_with_profiles(args)
else:
  print("available profile mesh type '" + args.type + "' not available!")
  sys.exit()
  mesh_name = args.type + "_" + args.profile
  
res_name = '_' + str(args.res)

file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"

