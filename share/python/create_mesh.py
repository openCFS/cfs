#!/usr/bin/env python
# this is to prevent annoying  
# /usr/lib64/python3.6/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
#  from ._conv import register_converters as _register_converters
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from mesh_tool import *
import argparse

# for inclusion_opverlapp finds inclusion_size for create_2d_mesh via bisection
def find_inclusion_overlap(args):
  if args.inclusion_overlap < 0.0 or args.inclusion_overlap > 1.0:
    print("inclusion_overlap shall be within 0 ... 1")  
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
    print("overlap for inclusion_size " + str(radius) + " is " + str(overlap) + " -> error " + str(err) + "\n")
    
    if err > 0:
      upper = radius
    else:
      lower = radius
  
  return mesh    

# print sys.argv

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--y_res', help="y-discretization of bulk2d and bulk3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--z_res', help="z-discretization of bulk3d and matlab3d for quadratic/ cubic elements", type=int, required = False )
parser.add_argument('--width', help="width in m", type=float, default = 1.0)
parser.add_argument('--height', help="optional height in m", type=float, required = False)
parser.add_argument('--depth', help="optional depth in m", type=float, required = False)
parser.add_argument('--type', help="predefined mesh type", choices=['bulk2d', 'bulk3d', 'matlab3d', 'cantilever2d', 'cantilever2d_reinforced', 'cantilever3d', 'lbm2d', 'lbm3d', 'msfem_two_load', 'two_load', 'validation_test', 'force_inverter', 'force_inverter_half', 'gripper', 'gripper_half', 'voxels_from_optistruct', 'convert_optistruct', 'traegerblz', 'box_lufo', 'mbb'], required = True)
parser.add_argument('--lbm', help="subtype for 'lbm'", choices=['two_inlet_one_outlet', 'pipe_bend', 'pipe', 'distributor', 'backstep', 'diffuser', 'two_inlet_two_outlet', 'low_in_high_out'])
parser.add_argument('--patch', help="define many regions", choices=['3x3', '4x4'])
parser.add_argument('--inclusion', help="inclusion for bulk2d and bulk3d", choices=["rect", "ball", "top_panel"])
parser.add_argument('--inclusion_size', help="possible mandatoryy size for inclusion as fraction of x-dimension (.9 is almost full)", type=float)
parser.add_argument('--inclusion_overlap', help="alternative to inclusion_size for ball. Give fraction of overlapping to boundary", type=float)
parser.add_argument('--file', help="optional give output file name. ")
parser.add_argument('--optistruct', help="optistruct file name")
parser.add_argument('--optistruct_type', help="optistruct mesh type", choices=['cell_opt', 'apod6', 'lufo_bracket'], default='cell_opt')
parser.add_argument('--optistruct_scaling', help="optistruct scaling factor for unit conversion.", type=float, default=1.)
parser.add_argument('--pfem', help="sets additional boundary elements for b.c.", action='store_true', default=False)
parser.add_argument('--numbering', help="numbering of nodes and elements (only 2D for now)", choices=['row_major', 'col_major'], default='row_major')


args = parser.parse_args()

mesh_name = args.type

# sanity checks
if args.lbm and not (args.type == "lbm2d" or "lbm3d"):
  print("error: --lbm only for --type lbm2d or lbm3d")
  sys.exit()
  
if (args.inclusion or args.inclusion_size or args.inclusion_overlap) and (args.type == "bulk3d" or args.type ==  "cantilever2d_reinforced"): 
  print("inclusions currently not for your type implemented") 
  
if args.inclusion_overlap and args.type != "bulk2d":
  print("--inclusion_overlap only for bulk2d implemented")
  sys.exit()  
  
if args.inclusion == 'ball' and (args.type == "lbm3d" or args.type == "lbm2d"):
  print("inclusion ball not implemented yet for lbm meshes")
  sys.exit()
  
if (args.inclusion and not (args.inclusion_size or args.inclusion_overlap)) or ((args.inclusion_size or args.inclusion_overlap) and not args.inclusion):
  print("inclusions require both --inclusion and --inclusion_size or --inclusion_overlap")  
  sys.exit()  

if args.inclusion and args.patch:
  print("--inclusion and --patch don't go concurrently")
  sys.exit() 
  
if args.type == "voxels_from_optistruct" or args.type == "convert_optistruct":
  if args.optistruct == None:
    print("Need optistruct file for conversion")
    sys.exit()
    
mesh= None 
    
if args.type == 'bulk3d' or args.type == 'validation_test' or args.type == 'cantilever3d' or args.type == 'traegerblz' or args.type == 'box_lufo':
  mesh = create_3d_mesh(args.type, args.res, args.y_res, args.z_res, args.width, args.height, args.depth, args.inclusion, args.inclusion_size,scale=args.width,pfem=args.pfem) 
elif args.type == 'matlab3d':
  mesh = create_3d_matlab_mesh(args.type, args.res, args.y_res, args.z_res, args.width, args.height, args.depth)  
elif args.type.startswith('lbm'):
  if args.lbm == None:
    print('error: --lbm subtype mandatory for --type lbm')
    sys.exit()
  if args.type == 'lbm2d':
     mesh = create_lbm2d(args.res, args.lbm, args.inclusion, args.inclusion_size)
  elif args.type == 'lbm3d' and args.lbm == 'backstep':
    mesh = create_backstep(args.res, args.y_res, args.z_res)
  elif args.type == 'lbm3d':
     mesh = create_lbm3d(args.res, args.y_res, args.z_res, args.lbm, args.inclusion, args.inclusion_size)
  else:
     mesh = create_3d_mesh(args.res, args.y_res, args.z_res)
  mesh_name = args.type +"_" + args.lbm
elif args.type == '3D':
  mesh = create_regular3d_mesh(args.type, args.res)
elif args.type == 'voxels_from_optistruct':
  mesh = voxelize_mesh_from_optistruct(args.optistruct, args.res)
  mesh_name = args.optistruct[:-4] + "_voxelized"
elif args.type == 'convert_optistruct':
  mesh = create_mesh_from_optistruct(args.optistruct, args.optistruct_scaling, args.optistruct_type)
  mesh_name = args.optistruct
else: # default case 2d_mesh
  if not args.inclusion_overlap:  
    mesh = create_2d_mesh(args.type, args.res, args.y_res, args.width, args.height, args.inclusion, args.inclusion_size, args.patch, args.pfem, args.numbering)
  else:
    mesh = find_inclusion_overlap(args) 
  
res_name = '_' + str(args.res)
if (args.type == 'bulk2d' or args.type == 'bulk3d' or args.type =='matlab3d') and args.y_res != None:
  res_name += '_' + str(args.y_res)
if (args.type == 'bulk3d' or args.type == 'matlab3d') and args.z_res:
  res_name += '_' + str(args.z_res)
if args.width != 1.0:
  res_name += '-w_' + str(args.width).replace('.', '_')
if args.height is not None:
  res_name += '-h_' + str(args.height).replace('.', '_')
if args.depth is not None:
  res_name += '-d_' + str(args.depth).replace('.', '_')
if args.inclusion_size:
  res_name += '_' + args.inclusion  + '_' + str(args.inclusion_size).replace('.', '_')
if args.inclusion_overlap:
  res_name += '_' + args.inclusion + '_ol_' + str(args.inclusion_overlap).replace('.', '_')
if args.patch:
  res_name += '_' + args.patch
if args.type == 'convert_optistruct':
  res_name = ""  
if args.pfem:
  res_name += "_pfem" 

file = mesh_name + res_name + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file, args.optistruct_scaling)
print("created file '" + file + "' with " + str(len(mesh.elements)) + " elements")

