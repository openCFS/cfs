#!/usr/bin/env python
import argparse
from draw_profile_functions import *
from mesh_tool import *

def create_mesh_with_profiles(x1, x2, y1, y2, z1, z2, xres, yres, zres,ipo):
  nx = xres
  ny = yres if yres <> None else xres 
  nz = zres if zres <> None else xres
  
  if ipo == "linear":    
    array = set_lin_profiles_array(nx, ny, nz, x1, x2, y1, y2, z1, z2)
  if ipo == "cubic":
    array = set_cubic_profiles_array(nx, ny, nz, x1, x2, y1, y2, z1, z2)
  if z1 == 0.0 and z2 == 0.0:
    mesh = create_2d_mesh_from_array(array)
  else:
    mesh = create_3d_mesh_from_array(array)
  
#   visualize_structure(array,nx,ny,nz)
  
  return mesh

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--y_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--z_res', help="y-discretization of bulk2s and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--type', help="predefined mesh type", choices=['profiles2d', 'profiles3d'], required = True)
parser.add_argument('--file', help="optional give output file name. ")
# for profile functions
parser.add_argument('--x1', help="first radius for profile of bar in x-direction; 0 <= x1 <= 1", type=float, required = False)
parser.add_argument('--x2', help="second radius for profile of bar in x-direction; 0 <= x2 <= 1", type=float, required = False)
parser.add_argument('--y1', help="first radius for profile of bar in y-direction; 0 <= y1 <= 1", type=float, required = False)
parser.add_argument('--y2', help="second radius for profile of bar in y-direction; 0 <= y2 <= 1", type=float, required = False)
parser.add_argument('--z1', help="first radius for profile of bar in z-direction; 0 <= z1 <= 1", type=float, required = False)
parser.add_argument('--z2', help="second radius for profile of bar in z-direction; 0 <= z2 <= 1", type=float, required = False)
parser.add_argument('--interpolation', help="interpolation type of profile functions (linear or cubic)", choices=["linear","cubic"], required = False)

args = parser.parse_args()

mesh_name = args.type

# sanity checks
if args.type.startswith("profiles") and not args.interpolation:
  print("interpolation type for profile functions required!")
  sys.exit()
  
if args.type == "profiles2d" and not (args.x1 and args.x2 and args.y1 and args.y2) and (args.z1 or args.z2):
  print("error: profiles2d needs values for x1,x2 and y1,y2")
  sys.exit()

if args.type == "profiles3d" and not (args.x1 and args.x2 and args.y1 and args.y2 and args.z1 and args.z2):
  print("error: profiles3d needs values for x1,x2 and y1,y2 and z1,z2")
  sys.exit()
  
if args.type == 'profiles2d':
  mesh = create_mesh_with_profiles(args.x1, args.x2, args.y1, args.y2, 0, 0, args.res, args.y_res, 1,args.interpolation)
elif args.type == 'profiles3d':  
  mesh = create_mesh_with_profiles(args.x1, args.x2, args.y1, args.y2, args.z1, args.z2, args.res, args.y_res, args.z_res,args.interpolation)
else:
  print("available profile mesh type '" + args.type + "' not available!")
  sys.exit()
  mesh_name = args.type + "_" + args.interpolation
  
res_name = '_' + str(args.res)

file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"

