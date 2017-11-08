#!/usr/bin/env python
# import mesh_tool first because of h5py

import matplotlib
#matplotlib.use('tkagg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import mesh_tool
import cfs_utils
import argparse
import draw_profile_functions
import numpy as np
from scipy import interpolate
import matviz_vtk
import vtk
import sys

class Surface_Mesh():
  def __init__(self,points,triangles):
    self.points = np.asarray(points)
    self.cells = triangles     
    
  def setConnectivity(self):
    # for each point ( id = array index), store tuple with all neighboring id nodes
    self.connectivity = [set() for index in range(len(self.points))]  
    for t in self.cells:
      
      self.connectivity[t[0]].add(t[1])
      self.connectivity[t[0]].add(t[2])
      
      self.connectivity[t[1]].add(t[0])
      self.connectivity[t[1]].add(t[2])
      
      self.connectivity[t[2]].add(t[0])
      self.connectivity[t[2]].add(t[1])
      
#     print("point 0  with coords ", self.points[0]," has neighbors ",self.connectivity[0])  
#     
#     sys.exit()   
  def smooth(self,niter):
    self.setConnectivity()
    lamb = 0.6
    for i in range(niter):
      self.laplacian_smoothing(lamb)   
      self.laplacian_smoothing(-lamb-0.03)
      
    print("Taubin smoothing with ", niter, " iterations")  
  
  # laplacian smoothing: p_i = p_i + \lambda * L(p_i)
  # using weighted average: L(p_i) = (w_ij*p_j + w_ik*p_k) / (w_ik+w_ik) - p_i, assuming neighbors are p_j,p_k   
  def laplacian_smoothing(self,lamb):
    new_points = [None] * len(self.points) 
    for i,p in enumerate(self.points):
      # calculate gradient   
      if np.isclose(p[0], 0) or np.isclose(p[1], 0) or np.isclose(p[2], 0) or np.isclose(p[0], 1.0) or np.isclose(p[1], 1.0) or np.isclose(p[2], 1.0):
        new_points[i] = p
      else:
        # calculate L(p_i)
        # w_ij*p_j + w_ik*p_k
        num = np.asarray([0,0,0])
        denom = 0
        L = 0
        # nid is id of a neighbor node
        neighborhood = self.connectivity[i] 
        for nid in neighborhood:
          # n are coords of neighbor with id nid
          n = self.points[nid]
          # distance between neighbor and this node
          w = 1.0 / np.linalg.norm(len(neighborhood))
          L = L + w * (n-p)
          #w = numpy.linalg.norm(n - p)
          #num = num + w * n
          #denom  = denom + w
            
        
        #L(p_i) = (w_ij*p_j + w_ik*p_k) / (w_ik+w_ik) - p_i
#         assert(not np.isclose(denom,0,1e-6))
#         L = num / denom - p   
        new_points[i] = p + lamb * L   
    
    self.points = new_points

def calc_volume(array,infoXml=None):
  res, res, res = array.shape
  
  elems = np.where(array != -1,1,0).sum() # np.where() delivers array with info on if condition <> -1 is fulfilled
  
  vol = float(elems)/float(res**3)
  
  if infoXml != None:
    infoXml.write('  <volume value="' + str(vol) + '"/>\n')
  
  print("volume:" + str(vol))
  
  return vol

def visualize_structure(points,cells,show,save):
  polydata = matviz_vtk.fill_vtk_polydata(points,cells)
  
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

if __name__ == "__main__":
#   import doctest, draw_profile_functions
#   doctest.testmod(draw_profile_functions)
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
  parser.add_argument('--input', help="thickness bascell parameters for lazy people, e.g. x1,x2,y1,y2,z1,z2", default="")
  parser.add_argument('--bend', help="bending factor for spline (0-1)", type=float, default=0.5)
  parser.add_argument('--skip_x', help="don't show bar in x direction", action='store_true')
  parser.add_argument('--skip_y', help="don't show bar in y direction", action='store_true')
  parser.add_argument('--skip_z', help="don't show bar in z direction", action='store_true')
  parser.add_argument('--show', help="show final structure in new window", action='store_true')
  parser.add_argument('--multiple_regions', help="create mesh with only one region", action='store_true', default=False)
  parser.add_argument('--verbose', help="show spline plots",choices=["off","all_bisecs","profile_map","polar_plot","interpolation","all_splines"], default="off")
  parser.add_argument('--plot_bisec', help="plot a bisec function {x,y,z}{0...8}, e.g. x7")
  parser.add_argument('--target', help="what to generate",choices=["volume_mesh","3dlines","None","surface_mesh"], required=True)
  parser.add_argument('--save', help="overwrite default target name")
  parser.add_argument('--save_vtp', help="write volume mesh data to .vtp file", action='store_true',default=False)
  parser.add_argument('--to_info_xml', help="writes information on profile funcs to .info.xml", action='store_true', default=False)
  parser.add_argument('--export', help="export different stuff", choices=['radius_maps','surface_points'], required=False)
  parser.add_argument('--force_bisec', help="take given bisec curve", choices=['bicubic','bspline','linear','heaviside'], required=False)
  parser.add_argument('--beta', help="steepness of heaviside function", type=float, default=10)
  parser.add_argument('--eta', help="midpoint heaviside function", type=float, default=0.5)
  parser.add_argument('--interpolation', help="interpolation type between splines and bisecs", choices=['linear','heaviside'], default="linear")
  parser.add_argument('--stiffness_as_diameter',help="interprete values for x1, x2, y1, ... directly as radii", action='store_true',default=False,required=False)
  parser.add_argument('--tets', help="tetrahedralize surface mesh", action='store_true',default=False)
  parser.add_argument('--smooth_iter', help="number of steps for Taubin's surface smoothing",type=int, default=50)
  
  args = parser.parse_args()
  
  parser.add_argument('--stiffness_as_radius',help="interprete values for x1, x2, y1, ... directly as radii", action='store_true',default=False,required=False)

  if args.res_surf_lines is None:
    args.res_surf_lines = args.res
  
  if args.stiffness is not None:
    args.x1 = args.stiffness
    args.x2 = args.stiffness
    args.y1 = args.stiffness
    args.y2 = args.stiffness
    args.z1 = args.stiffness
    args.z2 = args.stiffness
  elif args.input:
    n = args.input.split(',')
    if len(n) != 6:
      print("Error: Need 6 basecell parameters x1,x2,...!")
      sys.exit(1)
    args.x1 = float(n[0])
    args.x2 = float(n[1])
    args.y1 = float(n[2])
    args.y2 = float(n[3])
    args.z1 = float(n[4])
    args.z2 = float(n[5])
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
  
  if args.target == "volume_mesh":
    args.save_vtp = True

  print("stiffnesses: "  +str(args.x1) + "," + str(args.x2) + "," + str(args.y1) + "," + str(args.y2) + "," + str(args.z1) + "," + str(args.z2))  
    
  # sanity checks
  if not (args.x1 and args.x2 and args.y1 and args.y2 and args.z1 and args.z2):
    raise Exception("error: need values for x1,x2 and y1,y2 and z1,z2!")  
    
  print("radii: " + str(args.x1/2.0) + "," + str(args.x2/2.0) + "," + str(args.y1/2.0) + "," + str(args.y2/2.0) + "," + str(args.z1/2.0) + "," + str(args.z2/2.0))
  
  ###################### file name stuff ############
  fileNameBase = None
  if args.save is None: # set default gid mesh name
    if args.force_bisec:
      fileNameBase = "basecell_interp_" + args.interpolation + "_force_" + args.force_bisec
    else:
      fileNameBase = "basecell_interp_" + args.interpolation
      if args.interpolation == 'heaviside':
        fileNameBase += "_beta_" + str(args.beta) + "_eta_" + str(args.eta)
      
    assert(fileNameBase is not None)
    if args.stiffness_as_diameter:
      fileNameBase += "_diam_" + str(args.x1)
    else:
      fileNameBase += "_stiff_" + str(args.x1)
    
    if not (args.x2 == val and args.y1 == val and args.y2 == val and args.z1 == val and args.z2 == val):
      fileNameBase += "_" + str(args.x2) + "_" + str(args.y1) + "_" + str(args.y2) + "_" + str(args.z1) + "_" + str(args.z2)
    
    fileNameBase += "_bend_" + str(args.bend) + "_res_" + str(args.res)
    
    fileNameBase += "_skip_x" if args.skip_x else ""
    fileNameBase += "_skip_y" if args.skip_y else ""
    fileNameBase += "_skip_z" if args.skip_z else ""
    fileNameBase += "_" + args.target
    args.save = fileNameBase  
  else:
    if args.save.endswith(".stl") or args.save.endswith(".vtp"):
      fileNameBase = args.save[:-4]
    else:
      fileNameBase = args.save  
  
  ########### stuff for logging with .info.xml #########################
  infoStr = None
  infoXml = None
  
  if args.to_info_xml:
    # command line
    cmd = sys.argv[0].split('/')[-1]
    for i in range(1, len(sys.argv)):
      cmd += ' ' + sys.argv[i] 
    
    infoXmlName = fileNameBase + ".info.xml"
    infoXml = open(infoXmlName,"w") 
    infoXml.write('<?xml version="1.0"?>\n\n')
    infoXml.write('<basecell nx="' + str(args.res) + '" ny="' + str(args.res) + '" nz="' + str(args.res) +'">\n')
    infoXml.write('  <cmd value="' + cmd + '"/>\n')
    infoXml.write('  <input x1="' + str(args.x1) + '" x2="' + str(args.x2) + '" y1="' + str(args.y1) + '" y2="' + str(args.y2) + '" z1="' + str(args.z1) + '" z2="' + str(args.z2) + '" interpolation="' + args.interpolation + '"/>\n')
      
  if not args.stiffness_as_diameter:
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
  else:
    infoStr = '  <radii rx1="' + str(args.x1/2.0) + '" '
    infoStr += ' rx2="' + str(args.x2/2.0) + '" '
    infoStr += ' ry1="' + str(args.y1/2.0) + '" '
    infoStr += ' ry2="' + str(args.y2/2.0) + '" '
    infoStr += ' rz1="' + str(args.z1/2.0) + '" '
    infoStr += ' rz2="' + str(args.z2/2.0) + '"'  
  
  if infoXml is not None:
    assert(infoStr)
    infoXml.write(infoStr + "/>\n\n")  
  
  ################### actual work starts here ##############################
  # we need voxel array for gid mesh writing 
  array, points, cells, volume = draw_profile_functions.generate_basecell(args,infoXml)
  
  if args.target.startswith("surface"):
    smesh = Surface_Mesh(points, cells)
    smesh.smooth(args.smooth_iter)
    points = smesh.points
    cells = smesh.cells
  
  ############### writing files ############################################
  #mesh = create_mesh_with_profiles(args,infoXml,log)
  mesh = None
  polydata = matviz_vtk.fill_vtk_polydata(points,cells)
  if args.show: # show it only
    print("starting visualization...")
    show_vtk(polydata, 1000, [], True)
  ################### take care of gid mesh ##############
  if args.target.startswith("volume"):
    if args.z1 == 0.0 and args.z2 == 0.0:
      mesh = create_2d_mesh_from_array(array)
    else:
      mesh = mesh_tool.create_3d_mesh_from_array(array,args.multiple_regions)
    
    mesh_tool.validate_periodicity(mesh)
  elif args.target.startswith("surface") or args.target.startswith("surface"):
    stlName = fileNameBase + ".stl"
    if args.tets: # create tetrahedralized volume mesh from surface description
      mesh = mesh_tool.create_volume_mesh_with_gmsh(stlName)    
  
  if args.target == "volume_mesh" or args.target == "surface_mesh" and args.tets:   
    file = fileNameBase + '.mesh'
    assert(file.endswith('.mesh'))
    
    assert(mesh is not None)
    mesh_tool.write_gid_mesh(mesh, file)
  
  ################ take care of stl and vtp files ##############  
  if args.target.startswith("volume"):
    vtpName = fileNameBase + ".vtp"
    matviz_vtk.show_write_vtk(polydata,1000,vtpName)
  elif args.target.startswith("surface"):
    stlName = fileNameBase + ".stl"
    # make sure normals are oriented consistently
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(polydata)
    normals.SetConsistency(1)
    normals.SetAutoOrientNormals(1)
    normals.Update()
    matviz_vtk.write_stl(normals.GetOutput(),stlName)
    if args.save_vtp:
      matviz_vtk.show_write_vtk(normals.GetOutput(),1000,args.save+".vtp")
  else:
    print("Ohohohoh....")
    sys.exit()
              
  if infoXml != None:
    infoXml.write('</basecell>')    
    
class Basecell_Data():
  x1 = x2 = y1 = y2 = z1 = z2 = bend = beta = eta = res = None
  def __init__(self,res,bend,x1,x2,y1,y2,z1,z2,interpolation,beta=None,eta=None,offset=0,target="surface",res_surf_lines=None,tets=False):
    self.res = res
    self.x1 = x1
    self.x2 = x2
    self.y1 = y1
    self.y2 = y2
    self.z1 = z1
    self.z2 = z2
    self.bend = bend
    self.tets = tets
    assert(interpolation == "linear" or interpolation == "heaviside")
    if interpolation == "heaviside":
      assert(beta is not None and eta is not None)
    
    eps = 1e-12  
    assert(x1 > eps and x2 > eps and y1 > eps and y2 > eps and z1 > eps and z2 > eps)  
      
    self.beta = beta
    self.eta = eta  
    self.interpolation = interpolation
    self.target = target
    
    assert(self.target == "volume_mesh" or "volume_mesh")
    
    # set debugging stuff 
    self.verbose = "off"
    self.skip_x = False
    self.skip_y = False
    self.skip_z = False
    self.plot_bisec = False
    self.force_bisec = False
    self.export = None
    self.save = None
    self.save_vtp = False
    self.to_info_xml = False
    self.stiffness_as_diameter = False
    
    if not res_surf_lines:
      self.res_surf_lines = res
    
class Basecell():
  # data is an object of type Basecell_Data()
  def __init__(self,data):
    assert(type(data) is Basecell_Data)
    self.data = data
    self.points, self.cells, self.volume = draw_profile_functions.generate_basecell(data,None,None)
    
    # xmin, ymin, zmin, xmax, ymax, zmax
    self.bounds = np.zeros(6)
    if data.target != "volume_mesh":
      self.update()
      
    self.center = np.asarray([0.5,0.5,0.5])
    
  def scale(self,scalex,scaley,scalez):
    for i in range(len(self.points)):
      self.points[i] = np.asanyarray(self.points[i]) * np.asarray([scalex,scaley,scalez])
    
    self.center = self.center * np.asarray([scalex,scaley,scalez])
  def translate(self,x,y,z):
    for i in range(len(self.points)):
      self.points[i] = np.asanyarray(self.points[i]) + np.asanyarray([x,y,z])
    
    self.center = self.center + np.asanyarray([x,y,z])
    
  # recalculate bounds after rescaling and translating
  def update(self):
    self.bounds[0] = np.min(np.asarray(self.points)[:,0])
    self.bounds[1] = np.min(np.asarray(self.points)[:,1])
    self.bounds[2] = np.min(np.asarray(self.points)[:,2])
    self.bounds[3] = np.max(np.asarray(self.points)[:,0])
    self.bounds[4] = np.max(np.asarray(self.points)[:,1])
    self.bounds[5] = np.max(np.asarray(self.points)[:,2])
    
  
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
        