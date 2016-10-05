#!/usr/bin/env python
import argparse
from draw_profile_functions import *
from mesh_tool import *
import numpy as np
import matviz_rot
from matviz_vtk import *
import vtk

def calc_volume(array):
  res, res, res = array.shape
  
  elems = np.where(array <> -1,1,0).sum() # np.where() delivers array with info on if condition <> -1 is fulfilled
  
  print "volume:", float(elems)/float(res**3)

def visualize_structure(array,singRegion,show,save):
  print "starting visualization..."
  # create vtk cells and points
  cells = vtk.vtkCellArray()
  points = vtk.vtkPoints()
  
#   count = 0
  
  res, res, res = array.shape
  
  h = 1.0 / res
  
  centers = []
  
  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):
        x = i * h + h / 2.0
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        
        if array[i,j,k] > 0:
          if i > 0 and j > 0 and k > 0 and i < res-1 and j < res-1 and k < res - 1:   
            if array[i-1,j,k] < 0 or array[i+1,j,k] < 0 or array[i,j-1,k] < 0 or array[i,j+1,k] < 0 or array[i,j,k-1] < 0 or array[i,j,k+1] < 0:
              centers.append([x,y,z])
          else:
            centers.append([x,y,z]) 
  
  create_centered_bars(cells,points,centers,[h,h,h])
  
  polydata = vtk.vtkPolyData()
  polydata.SetPoints(points)
  polydata.SetPolys(cells)

#   print count," surface elements"
  
  if save:
    show_write_vtk(polydata,1000,save)
  if show: 
    show_vtk(polydata, 1000, [], True)
  
def give_radiusFunction():
  r = np.linspace(0.5, 0.5*np.sqrt(2),100)
  # area F = circle area - 4*circle segment (outside of bounding box)
  alpha = 2.0 * np.arccos(0.5/r) # 0,5 is half of side length of unit cube
  areaSegment = 0.5 * r**2 * (alpha - np.sin(alpha)) # area of circle segment
  areaCircle = np.pi * r**2  # area of circle
  diff = areaCircle - 4.0 * areaSegment # effective area
  
#   plt.plot(r,diff)
#   plt.show()
  
  f = interpolate.interp1d(diff,r) # approximate inverse function by linear interpolation
  
  return f

# for given relative stiffness (from 0 to 1), calculate respective radius 
def calc_radius(stiff):
  val = 0
  if stiff <= 0.5**2*np.pi:
    val = 2*np.sqrt(stiff/np.pi)
  else:
    f = give_radiusFunction()
    val = 2*f(stiff)
  
#   print val/2.0  
  return val 

# def create_mesh_with_profiles(x1, x2, y1, y2, z1, z2, xres, yres, zres,ipo):
def create_mesh_with_profiles(args):
  print "stiffnesses: ",args.x1,args.x2,args.y1,args.y2,args.z1,args.z2
    # calculating radii in relation to given stiffnesses x1,x2,y1,...
  args.x1 = calc_radius(args.x1)
  args.x2 = calc_radius(args.x2)
  args.y1 = calc_radius(args.y1)
  args.y2 = calc_radius(args.y2)
  args.z1 = calc_radius(args.z1)
  args.z2 = calc_radius(args.z2)
  
  print "radii: ",args.x1/2.0,args.x2/2.0,args.y1/2.0,args.y2/2.0,args.z1/2.0,args.z2/2.0    
  
  array = create_profiles_array(args)
  
  if args.target.startswith("volume"):
    calc_volume(array)
  
  if args.z1 == 0.0 and args.z2 == 0.0:
    mesh = create_2d_mesh_from_array(array)
  else:
    mesh = create_3d_mesh_from_array(array,args.single_region)
    
  mesh = convert_to_sparse_mesh(mesh)
  
  validate_periodicity(mesh)
  
  if args.show and args.target.startswith("volume"):
    save = "volume.vtp" if not args.save else args.save
    assert(save.endswith('.vtp'))
    visualize_structure(array,args.single_region,args.show,save)  
  
  return mesh

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True )
parser.add_argument('--type', help="predefined mesh type", choices=['profiles2d', 'profiles3d'], required = True)
# for profile functions
parser.add_argument('--stiffness', help="stiffness for profile of bar in all directions (x1,x2,y1,...); in [0,1]", type=float)
parser.add_argument('--x1', help="first stiffness for profile of bar in x-direction; 0 <= x1 <= 1", type=float)
parser.add_argument('--x2', help="second stiffness for profile of bar in x-direction; 0 <= x2 <= 1", type=float)
parser.add_argument('--y1', help="first stiffness for profile of bar in y-direction; 0 <= y1 <= 1", type=float)
parser.add_argument('--y2', help="second stiffness for profile of bar in y-direction; 0 <= y2 <= 1", type=float)
parser.add_argument('--z1', help="first stiffness for profile of bar in z-direction; 0 <= z1 <= 1", type=float)
parser.add_argument('--z2', help="second stiffness for profile of bar in z-direction; 0 <= z2 <= 1", type=float)
parser.add_argument('--profile', help="type of profile functions", choices=["linear","circular","spline"])
parser.add_argument('--bend', help="bending factor for spline (0-1)", type=float, default=0.5)
parser.add_argument('--skip_x', help="don't show bar in x direction", action='store_true')
parser.add_argument('--skip_y', help="don't show bar in y direction", action='store_true')
parser.add_argument('--skip_z', help="don't show bar in z direction", action='store_true')
parser.add_argument('--show', help="show final structure in new window", action='store_true')
parser.add_argument('--single_region', help="create mesh with only one region", action='store_true')
parser.add_argument('--verbose', help="show spline plots",choices=["off","all_profiles","bisec","profile_map"], default='off')
parser.add_argument('--target', help="what to generate",choices=["volume_vtk","volume_mesh","3dlines","None","surface_mesh"], required=True)
parser.add_argument('--save', help="overwrite default target name")


args = parser.parse_args()

if args.stiffness <> None:
  args.x1 = args.stiffness
  args.x2 = args.stiffness
  args.y1 = args.stiffness
  args.y2 = args.stiffness
  args.z1 = args.stiffness
  args.z2 = args.stiffness
else:
  if args.x1 == None or args.y1 == None or args.z1 == None:
    print "Error:stiffness or x1, y1, z1 necessary!"
    sys.exit(1)

  if args.x2 == None:
    args.x2 = args.x1
    
  if args.y2 == None:
    args.y2 = args.y1
    
  if args.z2 == None:
    args.z2 = args.z1
    
val = args.x1
if args.x2 == val and args.y1 == val and args.y2 == val and args.z1 == val and args.z2 == val:
  mesh_name = args.type + "_" + args.profile + "_stiff_" + str(args.x1) + "_bend_" + str(args.bend) + "_" + str(args.res)
else: 
  mesh_name = args.type + "_" + args.profile + "_stiff_" + str(args.x1) + "_" + str(args.x2) + "_" + str(args.y1) + "_" + str(args.y2) + "_" + str(args.z1) + "_" + str(args.z2) + "_bend_" + str(args.bend) + "_" + str(args.res)

# sanity checks
if args.type == "profiles2d" and not (args.x1 and args.x2 and args.y1 and args.y2) and (args.z1 or args.z2):
  print("error: profiles2d needs values for x1,x2 and y1,y2")
  sys.exit()

if args.type == "profiles3d" and not (args.x1 and args.x2 and args.y1 and args.y2 and args.z1 and args.z2):
  print("error: profiles3d needs values for x1,x2 and y1,y2 and z1,z2")
  sys.exit()

if args.type == 'profiles2d':
  assert(False)
  #mesh = create_mesh_with_profiles(args)
elif args.type == 'profiles3d':
  mesh = create_mesh_with_profiles(args)

if args.target == 'volume_mesh':   
  file = mesh_name + '.mesh' if args.save == None else args.save
  assert(file.endswith('.mesh')) 

  write_gid_mesh(mesh, file)

  print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"



