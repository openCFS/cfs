#!/usr/bin/env python
# import mesh_tool first because of h5py
import mesh_tool
import argparse
from draw_profile_functions import *
import numpy as np
import matviz_rot
from matviz_vtk import *

def calc_volume(array,infoXml):
  res, res, res = array.shape
  
  elems = np.where(array != -1,1,0).sum() # np.where() delivers array with info on if condition <> -1 is fulfilled
  
  vol = float(elems)/float(res**3)
  
  if infoXml != None:
    infoXml.write('  <volume value="' + str(vol) + '"/>\n')
  
  print("volume:" + str(vol))

def visualize_structure(array,singRegion,show,save):
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
        
        if array[i,j,k] >= 0:
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
    print("starting visualization...")
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
def create_mesh_with_profiles(args,infoXml,log):
  print("stiffnesses: "  +str(args.x1) + "," + str(args.x2) + "," + str(args.y1) + "," + str(args.y2) + "," + str(args.z1) + "," + str(args.z2))
  
  mesh = None
  
  # calculating radii in relation to given stiffnesses x1,x2,y1,...
  args.x1 = calc_radius(args.x1)
  infoStr = '  <radii rx1="' + str(args.x1) + '" '
  args.x2 = calc_radius(args.x2)
  infoStr += ' rx2="' + str(args.x2) + '" '
  args.y1 = calc_radius(args.y1)
  infoStr += ' ry1="' + str(args.y1) + '" '
  args.y2 = calc_radius(args.y2)
  infoStr += ' ry2="' + str(args.y2) + '" '
  args.z1 = calc_radius(args.z1)
  infoStr += ' rz1="' + str(args.z1) + '" '
  args.z2 = calc_radius(args.z2)
  infoStr += ' rz2="' + str(args.z2) + '"'
  
  if infoXml is not None:
    assert(infoStr)
    infoXml.write(infoStr + "/>\n\n")
  
  print("radii: " + str(args.x1/2.0) + "," + str(args.x2/2.0) + "," + str(args.y1/2.0) + "," + str(args.y2/2.0) + "," + str(args.z1/2.0) + "," + str(args.z2/2.0))    
  
  array = generate_basecell(args, infoXml,log)
  
  if args.target.startswith("volume"):
    calc_volume(array,infoXml)
  
    if args.z1 == 0.0 and args.z2 == 0.0:
      mesh = create_2d_mesh_from_array(array)
    else:
      mesh = mesh_tool.create_3d_mesh_from_array(array,args.single_region)
      
    mesh_tool.validate_periodicity(mesh)
  
  if (args.show or args.target.startswith("volume")) and not args.target.startswith("surface") and not args.target.startswith("3dlines"):
    if args.save_vtp:
      save = "volume.vtp" if not args.save else args.save
      if not save.endswith('.vtp'):
        save += ".vtp"
      visualize_structure(array,args.single_region,args.show,save)
    elif args.show:
      visualize_structure(array,args.single_region,args.show,False)
  if infoXml != None:
    infoXml.write('</basecell>')  
  
  return mesh

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="x-discretization of length 1m", type=int, required = True)
parser.add_argument("--res_surf_lines", help="resolution for surface lines, must be <= 360", type=int)
# for profile functions
parser.add_argument('--stiffness', help="stiffness for profile of bar in all directions (x1,x2,y1,...); in [0,1]", type=float)
parser.add_argument('--x1', help="first stiffness for profile of bar in x-direction; 0 < x1 < 1", type=float)
parser.add_argument('--x2', help="second stiffness for profile of bar in x-direction; 0 < x2 < 1", type=float)
parser.add_argument('--y1', help="first stiffness for profile of bar in y-direction; 0 < y1 < 1", type=float)
parser.add_argument('--y2', help="second stiffness for profile of bar in y-direction; 0 < y2 < 1", type=float)
parser.add_argument('--z1', help="first stiffness for profile of bar in z-direction; 0 < z1 < 1", type=float)
parser.add_argument('--z2', help="second stiffness for profile of bar in z-direction; 0 < z2 < 1", type=float)
parser.add_argument('--bend', help="bending factor for spline (0-1)", type=float, default=0.5)
parser.add_argument('--skip_x', help="don't show bar in x direction", action='store_true')
parser.add_argument('--skip_y', help="don't show bar in y direction", action='store_true')
parser.add_argument('--skip_z', help="don't show bar in z direction", action='store_true')
parser.add_argument('--show', help="show final structure in new window", action='store_true')
parser.add_argument('--single_region', help="create mesh with only one region", action='store_true', default=True)
parser.add_argument('--verbose', help="show spline plots",choices=["off","all_profiles","bisec","profile_map","polar_plot","interpolation","all_splines"], default="off")
parser.add_argument('--target', help="what to generate",choices=["volume_vtk","volume_mesh","3dlines","None","surface_mesh"], required=True)
parser.add_argument('--save', help="overwrite default target name")
parser.add_argument('--save_vtp', help="write volume mesh data to .vtp file", action='store_true',default=False)
parser.add_argument('--to_info_xml', help="writes information on profile funcs to .info.xml", action='store_true', default=False)
parser.add_argument('--export', help="export different stuff", choices=['radius_maps','surface_points'], required=False)
parser.add_argument('--force_bisec', help="take given bisec curve", choices=['bicubic','bspline','linear','heaviside'], required=False)
parser.add_argument('--beta', help="steepness of heaviside function", type=float, default=10)
parser.add_argument('--eta', help="midpoint heaviside function", type=float, default=0.5)
parser.add_argument('--logging',help="print logging while fixing surface gaps to log_fix_surface_gaps.txt", action='store_true',default=False,required=False)
parser.add_argument('--interpolation', help="interpolation type between splines and bisecs", choices=['linear','heaviside'], default="linear")


args = parser.parse_args()

if args.res_surf_lines is None:
  args.res_surf_lines = args.res

if args.stiffness is not None:
  args.x1 = args.stiffness
  args.x2 = args.stiffness
  args.y1 = args.stiffness
  args.y2 = args.stiffness
  args.z1 = args.stiffness
  args.z2 = args.stiffness
else:
  if args.x1 == None or args.y1 == None or args.z1 == None:
    print("Error:stiffness or x1, y1, z1 necessary!")
    sys.exit(1)

  if args.x2 == None:
    args.x2 = args.x1
    
  if args.y2 == None:
    args.y2 = args.y1
    
  if args.z2 == None:
    args.z2 = args.z1
    
val = args.x1

meshName = None
if args.save is None: # set default gid mesh name
  if args.force_bisec:
    meshName = "basecell_interp_" + args.interpolation + "_force_" + args.force_bisec
  else:
    meshName = "basecell_interp_" + args.interpolation
    
  assert(meshName is not None)
  meshName += "_stiff_" + str(args.x1)
  
  if not (args.x2 == val and args.y1 == val and args.y2 == val and args.z1 == val and args.z2 == val):
    meshName += "_" + str(args.x2) + "_" + str(args.y1) + "_" + str(args.y2) + "_" + str(args.z1) + "_" + str(args.z2)
  
  meshName += "_bend_" + str(args.bend) + "_" + str(args.res)  
else:
  meshName = args.save  

infoXml = None

if args.to_info_xml:
  # command line
  cmd = sys.argv[0].split('/')[-1]
  for i in range(1, len(sys.argv)):
    cmd += ' ' + sys.argv[i] 
  
  infoXmlName = meshName + ".info.xml"
  infoXml = open(infoXmlName,"w") 
  infoXml.write('<?xml version="1.0"?>\n\n')
  infoXml.write('<basecell nx="' + str(args.res) + '" ny="' + str(args.res) + '" nz="' + str(args.res) +'">\n')
  infoXml.write('  <cmd value="' + cmd + '"/>\n')
  infoXml.write('  <input x1="' + str(args.x1) + '" x2="' + str(args.x2) + '" y1="' + str(args.y1) + '" y2="' + str(args.y2) + '" z1="' + str(args.z1) + '" z2="' + str(args.z2) + '"/>\n')

log = None 

if args.logging:
  log = open(meshName+".log","w")
  
# sanity checks
if not (args.x1 and args.x2 and args.y1 and args.y2 and args.z1 and args.z2):
  print("error: need values for x1,x2 and y1,y2 and z1,z2!")
  sys.exit(1)

mesh = create_mesh_with_profiles(args,infoXml,log)

if args.target == 'volume_mesh':   
  file = meshName + '.mesh'
  assert(file.endswith('.mesh'))
  
  mesh_tool.write_gid_mesh(mesh, file)

  print("created file '" + file + "' with " + str(len(mesh.elements)) + " elements")

############## info xml scheme #####################
# <basecell>
# <input x1="" x2="" y1="" y2="" z1="" z2=""/>
# <radii rx1="" rx2="" ry1="" ry2="" rz1="" rz2=""/>
# <profile type="{bspline,bisectionSpline,linear}" dir"{x,y,z}">
#  <bSpline rad1="" rad2="" bend="">
#    <controlPolgygon>
#      <P1 x="" y=""/>
#      <P2 x="" y=""/>
#      <P3 x="" y=""/>
#      <P4 x="" y=""/>
#    </controlPolgygon>
#  </bSpline>
#  <bisectionSpline type="{biquadratic,bspline,linear}" angle="" >
#
#  </bisectionSpline>
# </profile>
# /<basecell>