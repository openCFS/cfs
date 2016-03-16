#!/usr/bin/env python
from mesh_tool import *
import argparse

# for inclusion_opverlapp finds inclusion_size for create_2d_mesh via bisection
def find_inclusion_overlap(args):
  if args.inclusion_overlap < 0.0 or args.inclusion_overlap > 1.0:
    print "inclusion_overlap shall be within 0 ... 1"  
    sys.exit()
    
  mesh = None  
  lower = 0.9
  upper = numpy.sqrt(2.0) # at least for square
  err = 10 # we normalize to pixtes
  ny = args.res if args.y_res is None else args.y_res
  while numpy.abs(err) > 1.0:
    radius = lower + (upper - lower)/2
    mesh = create_2d_mesh(args.type, args.res, args.y_res, args.width, args.height, args.inclusion, radius, args.patch)
    
    overlap = [item for item in mesh.bc if item[0] == 'left_iface']
    overlap = 0 if len(overlap) == 0 else len(overlap[0][1])
    #print overlap
    #print args.inclusion_overlap
    err = overlap - args.inclusion_overlap * ny 
    print "overlap for inclusion_size " + str(radius) + " is " + str(overlap) + " -> error " + str(err) + "\n"
    
    if err > 0:
      upper = radius
    else:
      lower = radius
  
  return mesh    


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
parser.add_argument('--inclusion_size', help="possible mandatoryy size for inclusion as fraction of x-dimension (.9 is almost full)", type=float)
parser.add_argument('--inclusion_overlap', help="alternative to inclusion_size for ball. Give fraction of overlapping to boundary", type=float)
parser.add_argument('--file', help="optional give output file name. ")

args = parser.parse_args()

mesh_name = args.type

# sanity checks
if args.lbm and not (args.type == "lbm2d" or "lbm3d"):
  print "error: --lbm only for --type lbm2d or lbm3d"
  sys.exit()
  
if (args.inclusion or args.inclusion_size or args.inclusion_overlap) and (args.type == "bulk3d" or args.type ==  "cantilever2d_reinforced"): 
  print "inclusions currently not for your type implemented" 
  
if args.inclusion_overlap and args.type <> "bulk2d":
  print "--inclusion_overlap only for bulk2d implemented"
  sys.exit()  
  
if args.inclusion == 'ball' and (args.type == "lbm3d" or args.type == "lbm2d"):
  print("inclusion ball not implemented yet for lbm meshes")
  sys.exit()
  
if (args.inclusion and not (args.inclusion_size or args.inclusion_overlap)) or ((args.inclusion_size or args.inclusion_overlap) and not args.inclusion):
  print "inclusions require both --inclusion and --inclusion_size or --inclusion_overlap"  
  sys.exit()  

if args.inclusion and args.patch:
  print "--inclusion and --patch don't go concurrently"
  sys.exit() 
  
if args.type == 'bulk3d':
  mesh = create_3d_mesh(args.type, args.res, args.y_res, args.z_res, args.inclusion, args.inclusion_size)  
elif args.type == 'bulk2d' or args.type.startswith('cantilever2d') or args.type == 'triangle_msfem' or args.type == 'pressure2' or args.type == 'msfem_two_load':
  if not args.inclusion_overlap:  
    mesh = create_2d_mesh(args.type, args.res, args.y_res, args.width, args.height, args.inclusion, args.inclusion_size, args.patch)
  else:
    mesh = find_inclusion_overlap(args) 
  
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
if args.inclusion_size:
  res_name += '_' + args.inclusion  + '_' + str(args.inclusion_size).replace('.', '_')
if args.inclusion_overlap:
  res_name += '_' + args.inclusion + '_ol_' + str(args.inclusion_overlap).replace('.', '_')
if args.patch:
  res_name += '_' + args.patch   

file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"

