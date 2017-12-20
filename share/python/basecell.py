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
from draw_profile_functions import Vertex
from draw_profile_functions import calc_distance
import numpy as np
from scipy import interpolate
import matviz_vtk
import vtk
import sys

# return list with length len(points)
# for each point, this list gives the ids of its neighbors
def getConnectivity(points,cells):
  # for each point ( id = array index), store tuple with all neighboring id nodes
  connectivity = [set() for index in range(len(points))]
  for i,t in enumerate(cells):
    assert(t[0] < len(points))
    assert(t[1] < len(points))
    assert(t[2] < len(points))
    connectivity[t[0]].add(t[1])
    connectivity[t[0]].add(t[2])
    
    connectivity[t[1]].add(t[0])
    connectivity[t[1]].add(t[2])
    
    connectivity[t[2]].add(t[0])
    connectivity[t[2]].add(t[1])
  
  return connectivity

def taubin_smoothing(points,connectivity,niter,bounds=None):
  assert(niter > 0)
  # smoothing parameter: p_i = p_i + lambda*L(p_{i,j})
  lamb = 0.8
  new_points = points
  for i in range(niter):
    new_points = laplacian_smoothing(laplacian_smoothing(new_points,connectivity,lamb),connectivity,-lamb-0.04,bounds)
    
#   print("Taubin smoothing with ", niter, " iterations")
    
  return new_points

# laplacian smoothing: p_i = p_i + \lambda * L(p_i)
# using weighted average: L(p_i) = (w_ij*p_j + w_ik*p_k) / (w_ik+w_ik) - p_i, assuming neighbors are p_j,p_k
# bounds: tuple/list with 6 entries - points on these boundaries are not smoothed   
# bounds order: xmin,ymin,zmin,xmax,ymax,zmax
def laplacian_smoothing(points,connectivity,lamb,bounds=None):
  new_points = [None] * len(points) 
  if bounds is None:
    bounds = [0.0,0.0,0.0,1.0,1.0,1.0] # default unit cube
  for i,vertex in enumerate(points):
    p = vertex.coords
    # don't smooth at the boundary
    if np.isclose(p[0], bounds[0]) or np.isclose(p[1], bounds[1]) or np.isclose(p[2], bounds[2]) or np.isclose(p[0], bounds[3]) or np.isclose(p[1], bounds[4]) or np.isclose(p[2], bounds[5]):
      assert(i == vertex.idx)
      new_points[i] = Vertex(p,vertex.idx)
    else:
      # calculate L(p_i)
      # w_ij*p_j + w_ik*p_k
      num = np.asarray([0,0,0])
      denom = 0
      L = 0
      # nid is id of a neighbor node
      neighborhood = connectivity[i] 
      for nid in neighborhood:
        # n are coords of neighbor with id nid
        n = points[nid].coords
        # distance between neighbor and this node
        w = 1.0 / np.linalg.norm(len(neighborhood))
        L = L + w * (n-p)
      new_points[i] = Vertex(p + lamb * L,vertex.idx)   
  
  for p in new_points:
    assert(type(p) is Vertex)
  return new_points  

def calc_volume(array,infoXml=None):
  res, res, res = array.shape
  
  elems = np.where(array != -1,1,0).sum() # np.where() delivers array with info on if condition <> -1 is fulfilled
  
  vol = float(elems)/float(res**3)
  
  if infoXml != None:
    infoXml.write('  <volume value="' + str(vol) + '"/>\n')
  
  print("volume:" + str(vol))
  
  return vol

def visualize_structure(points,cells,show,save):
  polydata, _ = matviz_vtk.fill_vtk_polydata(points,cells)
  
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
  parser.add_argument('--target', help="what to generate",choices=["volume_mesh","3dlines","marching_cubes","surface_mesh"], required=True)
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
  parser.add_argument('--smooth_iter', help="number of steps for Taubin's surface smoothing",type=int, default=30)
  
  
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
  
  
  ############# debugging ###########################
  #read non-design coords from csv file
#   import csv
#   nondes = []
#   with open("nondes.csv","r") as f:
#     reader = csv.reader(f,delimiter=',')
#     for i,row in enumerate(reader):
#       # skip header
#       if i == 0:
#         continue
#       nondes.append((float(row[0]),float(row[1]),float(row[2])))
#   
#   left_front_corner = [0, -163, 110.5]
#   right_back_corner = [290,163,206.5]
      
  ################### actual work starts here ##############################
  # we need voxel array for gid mesh writing 
  args.bc_flags = None
  array, points, cells, _ = draw_profile_functions.generate_basecell(args,infoXml,None)
  volume = calc_volume(array, infoXml)
  
  if args.target == "surface_mesh":
    connectivity = getConnectivity(points,cells)
    points = taubin_smoothing(points,connectivity,20) 
  
  ############### writing files ############################################
  #mesh = create_mesh_with_profiles(args,infoXml,log)
  mesh = None
  # --volume_mesh gives list of points with only coordinate entries
  # --surface_mesh gives list of Vertex() objects
  coords = points
  if args.target == "surface_mesh" or args.target == "marching_cubes":
    coords = [p.coords for p in points]
  polydata = matviz_vtk.fill_vtk_polydata(coords,cells)
  if args.show: # show it only
    print("starting visualization...")
    show_vtk(polydata, 1000, [], True)
  ################### take care of gid mesh ##############
  if args.target == "volume_mesh":
    if args.z1 == 0.0 and args.z2 == 0.0:
      mesh = create_2d_mesh_from_array(array)
    else:
      mesh = mesh_tool.create_3d_mesh_from_array(array,args.multiple_regions)
    
    mesh_tool.validate_periodicity(mesh)
    
    mesh_tool.write_gid_mesh(mesh, fileNameBase+".mesh") 
  
  ################ take care of stl and vtp files ##############  
  if args.target == "volume_mesh":
    vtpName = fileNameBase + ".vtp"
    matviz_vtk.show_write_vtk(polydata,1000,vtpName)
  elif args.target == "surface_mesh" or args.target == "marching_cubes":
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
    if args.tets: # create tetrahedralized volume mesh from surface description
      print("here")
      mesh_tool.create_volume_mesh_with_gmsh(stlName)
  else:    
    print("Error: Missing if branch for writing output!")
    sys.exit()
              
  if infoXml != None:
    infoXml.write('</basecell>')    
    
class Basecell_Data():
  x1 = x2 = y1 = y2 = z1 = z2 = bend = beta = eta = res = None
  def __init__(self,res,bend,x1,x2,y1,y2,z1,z2,interpolation,beta=None,eta=None,offset=0,target="surface_mesh",res_surf_lines=None,tets=False,bc_flags=None):
    self.res = res
    self.x1 = x1
    self.x2 = x2
    self.y1 = y1
    self.y2 = y2
    self.z1 = z1
    self.z2 = z2
    self.bend = bend
    self.tets = tets
    # flags = [0,nx-1,0,ny-1,0,nz-1]
    self.bc_flags = bc_flags
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
  # idx = (i,j,k)
  # nondes: list of non-design elements' barycenters
  # need global_grid_coords for nondes: base cell is embedded in global uniform grid - store left lower and right upper corner (i,j,k) for this cell
  # e.g. global_grid_coords[0] = (0,0,0) and global_grid_coords[1] = (40,40,40) if base cell voxel resolution is 40
  def __init__(self,data,idx=None,nondes=None,global_grid_coords=None):
    assert(type(data) is Basecell_Data)
    if global_grid_coords:
      assert(len(global_grid_coords) == 2)
      assert(len(global_grid_coords[0]) == 3)
    self.data = data
    self.idx = idx
    self.voxels, self.points, self.cells, self.boundary_points = draw_profile_functions.generate_basecell(data,None,nondes,global_grid_coords)
    assert(len(self.boundary_points) > 0)
    assert(self.points is not None)
    assert(self.cells is not None)
    for p in self.points:
      assert(type(p) is Vertex)
    for list in self.boundary_points:
      for p in list:
        assert(type(p) is Vertex)  
#     print("self.points:",self.points,"\n")

    # xmin, ymin, zmin, xmax, ymax, zmax
    self.bounds = np.zeros(6)
    if data.target != "volume_mesh":
      self.update()
    if data.target == "surface_mesh":
      connectivity = getConnectivity(self.points,self.cells)
      self.points = taubin_smoothing(self.points,connectivity,20)  
      
    self.center = np.asarray([0.5,0.5,0.5])
    
  def scale(self,scalex,scaley,scalez):
    scale = np.asarray([scalex,scaley,scalez])
    self.points = [Vertex(p.coords*scale,p.idx) for p in self.points]
    new_bp = []
    for list in self.boundary_points:
      new_list = [Vertex(p.coords*scale,p.idx) for p in list]
      new_bp.append(new_list)
      
    self.boundary_points = new_bp 
    assert(len(self.boundary_points) > 0)
    self.center = self.center * np.asarray([scalex,scaley,scalez])
  def translate(self,x,y,z):
    shift = np.asarray([x,y,z])
    self.points = [Vertex(p.coords+shift,p.idx) for p in self.points]
    new_bp = []
    for list in self.boundary_points:
      new_list = [Vertex(p.coords+shift,p.idx) for p in list]
      new_bp.append(new_list)
      
    self.boundary_points = new_bp 
    assert(len(self.boundary_points) > 0)
    self.center += shift
    
  # recalculate bounds after rescaling and translating
  def update(self):
    coords = np.asarray([p.coords for p in self.points])
    self.bounds[0] = np.min(coords[:,0])
    self.bounds[1] = np.min(coords[:,1])
    self.bounds[2] = np.min(coords[:,2])
    self.bounds[3] = np.max(coords[:,0])
    self.bounds[4] = np.max(coords[:,1])
    self.bounds[5] = np.max(coords[:,2])
  
  # take list with 'new' coords and replace mylist at entry list_idx  
  def replace_boundary_points(self,new,list_idx):
    print("\n start replacing")
    for i in range(len(self.points)):
      if self.points[i].idx != i:
        print("self.points[i].idx != i  idx:",self.points[i].idx,"  i:",i)
      assert(self.points[i].idx == i)
    print("\nlist_idx:",list_idx)
    for p in self.boundary_points[list_idx]:
      print(p)
    
#     if len(new) != len(self.boundary_points[list_idx]):  
#       out = open("cells.txt","w")
#       for c in self.cells:
#         out.write(str(c) + "\n")
#       out.close()    

    # if number of nodes on both sides differ, we have hanging nodes
    while len(new) != len(self.boundary_points[list_idx]):
      assert(np.abs(len(new)-len(self.boundary_points[list_idx])) < 10)
      # detect gap by finding first and last point (indices) of difference
      long_first_id = None
      long_last_id = None
      long_list = new if len(new) > len(self.boundary_points[list_idx]) else self.boundary_points[list_idx]
      short_list = new if len(new) < len(self.boundary_points[list_idx]) else self.boundary_points[list_idx]
      assert(long_list is not short_list)
      
      short_first_id = None
      short_last_id = None
      for i in range(len(long_list)):
    #     print("\nlong_list[i]:",long_list[i])
    #     print("short_list[i]:",short_list[i])
        if i < len(short_list) and calc_distance(long_list[i].coords,short_list[i].coords) > 1e-6:
          # set only once
          assert(long_first_id is None)
          long_first_id = i-1
          # store idx of last short_list element with match in long_list
          short_first_id = i-1
          print("long_first_id:",long_first_id," coords:",long_list[i-1].coords)
          print("short_first_id:",short_first_id," coords:",short_list[i-1].coords)
          break
      
      assert(short_first_id is not None)
      for i in range(short_first_id+1,len(short_list)):
        for j,p in enumerate(long_list):
    #       print("short_list[i]:",short_list[i])
    #       print("p:",p)
          if np.isclose(calc_distance(short_list[i].coords,p.coords),0):
            # found points where both lists match again
            short_last_id = i
            long_last_id = j
            print("long_last_id:",long_last_id," coords:",long_list[i].coords)
            print("short_last_id:",short_last_id," coords:",short_list[j].coords)
            break       
        if short_last_id is not None:
          break
      
      short_point_ids = []
      long_point_ids = []
      print("\n short_point_ids:")
      for p in range(short_first_id,short_last_id+1):
        short_point_ids.append(p)
        print(short_point_ids[-1],short_list[short_point_ids[-1]].idx)
      
      print("\n long_point_ids:")  
      for p in range(long_first_id,long_last_id+1):
        long_point_ids.append(p)
        print(long_point_ids[-1],long_list[long_point_ids[-1]].idx)
      
      assert(len(long_point_ids) > len(short_point_ids))
      # assert first and last points match
      assert(np.isclose(calc_distance(short_list[short_point_ids[0]].coords,long_list[long_point_ids[0]].coords),0))
      assert(np.isclose(calc_distance(short_list[short_point_ids[-1]].coords,long_list[long_point_ids[-1]].coords),0))
      
      merged = False
      for i in range(1,len(long_point_ids)-1):
        for j in range(1,len(short_point_ids)-1):
          # merge point if too close
          dist = calc_distance(long_list[long_point_ids[i]].coords,short_list[short_point_ids[j]].coords)
          if dist < 1e-3:
            print("point ",short_list[short_first_id+j], " close to ",long_list[long_first_id+i]," dist: ",dist)
            # don't replace ids to ensure connectivity of cells
            short_list[short_first_id+j].coords = long_list[long_first_id+i].coords
            merged = True
      
      # real hanging nodes
      if not merged:    
        assert(len(long_point_ids)-1 == len(short_point_ids))
#         print("\nlong list:")
#         for p in long_list:
#           print(p)
#         print("\nshort list:")
#         for p in short_list:
#           print(p)  
        
        print("insert ", long_list[long_first_id+1], " in short list at pos ", short_first_id+1)
        new_p = long_list[long_first_id+1]
        new_p.idx = len(self.points)   
        self.points.append(new_p)
        short_list.insert(short_first_id+1,new_p)
        # added additional point to shorter list, now create new cells
        left_id = short_list[short_first_id].idx
        right_id = short_list[short_first_id+2].idx
        new_id = new_p.idx
        cell = None
        print("left_id:",left_id," coords:",short_list[short_first_id].coords)
        print("right_id:",right_id, " coords:", short_list[short_first_id+2].coords)
        print("new_id:",new_id, " coords:",new_p.coords)
        # get cell with edge left-right
        for i,c in enumerate(self.cells):
          if left_id in c and right_id in c:
            # get cell and delete it directly
            cell = self.cells.pop(i)
            break
        
        assert(cell is not None)
        print("cell:",cell,self.points[cell[0]],self.points[cell[1]],self.points[cell[2]])
        for i,p in enumerate(self.points):
          if p.idx == cell[0]:
            print("cell[0] points list id:",id," point:",p)
            break
        for i,p in enumerate(self.points):
          if p.idx == cell[1]:
            print("cell[1] points list id:",id," point:",p)
            break
        for i,p in enumerate(self.points):
          if p.idx == cell[2]:
            print("cell[2] points list id:",id," point:",p)
            break    
        mid_id = [ vert for vert in cell if vert not in [left_id,right_id]][0]
        print("mid_id:",mid_id)
        assert(mid_id != left_id)
        assert(mid_id != right_id)
        # split cell into two cells
        self.cells.append((left_id,new_id,mid_id))
        self.cells.append((new_id,right_id,mid_id))
        print("appending cell ",(left_id,new_id,mid_id))
        print("appending cell ",(new_id,right_id,mid_id))
        
    bp = self.boundary_points[list_idx]
    # remember point ids that we replace
    old_ids = [p.idx for p in bp]
    # replace coords, keep ids
    for i,p in enumerate(new):
#       print("replace ",self.boundary_points[list_idx][i].coords, " with ", p.coords)
      self.boundary_points[list_idx][i].coords = p.coords
      
    new_points = [None] * len(self.points)
    # changing has also to be done in global points list
    for i,p in enumerate(self.points):
      # identify boundary point by id
      if p.idx in old_ids:
        # give new coords
        # find correct entry (index) in boundary list
        bp_id = old_ids.index(p.idx)
        new_points[p.idx] = new[bp_id]
        new_points[p.idx].idx = p.idx 
      else:
        # keep point if not replacing
        new_points[i] = p

    assert(len(self.points) == len(new_points))
     
    for i in range(len(new_points)):
      if np.linalg.norm(np.asarray(self.points[i].coords) - np.asarray(new_points[i].coords) > 1e-6):
        print("replace point ", self.points[i], " with point ", new_points[i].coords)
                 
    self.points = new_points
    
#     print("replaced ", len(old_ids), " points with ", len(new))      
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
        
