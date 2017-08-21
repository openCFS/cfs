#1/usr/bin/env python

import basecell
import draw_profile_functions
import numpy as np
import matviz_vtk
import pymp
import vtk
import math

try:
  import meshpy.triangle as triangle
#   import meshpy.tet as tet
except:
  print("Failed to load meshpy - need it for tetrahedralizing basecell mesh")

# similar to create_3d_cross_ip; # without rotation and shearing
# returns

def create_3d_interpretation_ortho(args,coords,s1,s2,s3,scale,samples,grad,thresh):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  # csize: size of one cell, e.g. [8,8,8]
  
  # point coordinates from h5 file
  centers, minimum, maximum = coords[0:3]
  
  print("min:",minimum," max:",maximum)
  
  # appendind cells
  appends = vtk.vtkAppendPolyData()
  
  if scale <= 0:
    scale = 1.0
  
  # assume we always start at (0,0,0)
  # order: min_x,min_y,min_z,max_x,max_y,max_z
  bounds = np.ones(6) * (-1)
  bounds[0] = bounds[1] = bounds[2] = 0
  bounds[3] = maximum[0] + minimum[0]
  bounds[4] = maximum[1] + minimum[1]
  bounds[5] = maximum[2] + minimum[2]
  
  # set size dx/dy/dz of one cell
  dx = bounds[3] / samples[0]
  dy = bounds[4] / samples[1]
  dz = bounds[5] / samples[2]
  
  min_thresh = 0.1
  max_thresh = 0.82
    
  delta = (abs(maximum[0] + minimum[0]), abs(maximum[1] + minimum[1]), abs(maximum[2] + minimum[2]))
  # where we want nodes
  nx = int(delta[0] / dx)
  ny = int(delta[1] / dy)
  nz = int(delta[2] / dz)    
  
  data_grid, data_grid_near, sample_coords, ndim, scale_ = matviz_vtk.get_interpolation_row_major(coords, grad, s1, s2, s3, dx, dy, dz)
  count = 0
  basecells = pymp.shared.list()
  with pymp.Parallel(args.bc_num_threads) as p:
    for id in p.range(nx*ny*nz):
      count = count + 1
      i, j, k = get_3d_grid_coords(id,nx,ny,nz)
      
      this = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k))
      east = get_interp_3darray_elem(data_grid,data_grid_near,(i+1,j,k))
      top = get_interp_3darray_elem(data_grid,data_grid_near,(i,j+1,k))
      front = get_interp_3darray_elem(data_grid,data_grid_near,(i,j,k+1))
      
      assert(this is not None and east is not None and top is not None and front is not None)

      # if one of the values is < min_thresh, set it to min_thresh        
      # if one of the values is > max_thresh, set it to max_thresh
      x1 = min(max(this[0],min_thresh),max_thresh)
      x2 = min(max(east[0],min_thresh),max_thresh)
      y1 = min(max(this[1],min_thresh),max_thresh)
      y2 = min(max(top[1],min_thresh),max_thresh)
      z1 = min(max(this[2],min_thresh),max_thresh)
      z2 = min(max(front[2],min_thresh),max_thresh)
      
      # translate cell to correct position
      left_front_corner  = np.asarray(sample_coords[i,j,k])
      
      bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta)
      bc_input.eta = 0.7
      bc_input.stiffness_as_diameter = True
      cell_obj = basecell.Basecell(bc_input)
      cell_obj.scale(dx, dy, dz)
      cell_obj.translate(left_front_corner[0],left_front_corner[1],left_front_corner[2])
      cell_obj.update()
      if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
        continue
      
      # flags for meshing circles on the boundary
      flags = [None] * 6
      grid_bounds = [0,0,0,nx-1,ny-1,nz-1]
      grid_coords = [i,j,k]
#       flags[0] = True if i == 0 else False
#       flags[1] = True if j == 0 else False
#       flags[2] = True if k == 0 else False
#       flags[3] = True if i == nx-1 else False
#       flags[4] = True if j == ny-1 else False
#       flags[5] = True if k == nz-1 else False
      for c in range(len(flags)):      
        flags[c] = True if grid_coords[c%3] == grid_bounds[c] else False 

      cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
      cell_obj.center = cell_center
      boundary_list = []
      # at least one boundary circle needs to be triangulated
      if any(flags):
        boundary_list = mesh_basecell_boundaries(flags,cell_obj,bounds,cell_center)
      with p.lock:
        basecells.append((cell_obj,boundary_list,flags))
        print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)
  
  flags = []
  # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
  for i,obj in enumerate(basecells):
    cell = obj[0]
    flags.append(obj[2])
    vtk_points = vtk.vtkPoints()
    pd = matviz_vtk.fill_vtk_polydata(cell.points, cell.cells)
    
    appends.AddInputData(pd)
    # meshed boundary circles
    circles = obj[1]
#     assert(len(circles) > 0)
    # list with meshes on the boundaries
    for l in circles:
      points = np.asarray(l[0])
      pd = matviz_vtk.fill_vtk_polydata(l[0], l[1])
        
      appends.AddInputData(pd)
      appends.Update()
      
  appends.Update() # not sure if we have to do this in each loop iteration
  
  # merge duplicated points etc.
  cleanFilter = vtk.vtkCleanPolyData()
  cleanFilter.SetInputConnection(appends.GetOutputPort())
  cleanFilter.Update()
  
  boundaryPts,loops = detect_boundary_edges(cleanFilter.GetOutput())
  
  appends.AddInputData(fill_boundary_loops(boundaryPts,loops))
  
  cleanFilter.Update()
  
  return cleanFilter.GetOutput()
    
# @param idx: tuple of three ints storing array indices(i,j,k)
# @param array: return element of array at position idx if exist
# if out of range, return None 
# @param fallback for array, if array elem at idx has value -1      
def get_interp_3darray_elem(array,fallback,idx):
  # must be a ndarray so that indexing with tuples works
  assert(type(array) == np.ndarray)
  assert(type(fallback) == np.ndarray)
#   assert(array.ndim == 3)
  nx, ny, nz, _ = array.shape
  
  if idx[0] >= nx or idx[0] < 0  or idx[1] >= ny or idx[1] < 0  or idx[2] >= nz  or idx[2] < 0 :
    return None
  else:
    if array[idx][0] == -1 or array[idx][1] == -1 or array[idx][2] == -1:
      return fallback[idx]
    else:
      return array[idx]

# @param id: idx in 1D array (row major)
# @return i,j,k: indices of equivalent 3D array
def get_3d_grid_coords(index,nx,ny,nz):
  i = index % nx
  j = (index - i)/nx % ny
  k = ((index - i)/nx-j)/ny
  
  return int(i),int(j),int(k)

# returns list of coordinates for given boundary 'bound'
# using order {0:"x_left", 1:"y_left", 2:"z_left", 3:"x_right", 4:"y_right", 5:"z_right"}
# @param bc_bounds: minimum/maximum coordinate value of interpretation domain
# order: min_x,min_y,min_z,max_x,max_y,max_z
# @param cc_3d: basecell center; has also to be projected onto right plane
def extract_2d_bc_boundary_coords(coords,bound,bc_bounds,cc_3d):
  result = []
  # depending on boundary flag, determine which coordinate component to compare
  # 0 -> 0; 1 -> 1; 2 -> 2; 3 -> 0; 4 -> 2; 5 -> 3
  major_dir = bound%3
  minor_1, minor_2 = draw_profile_functions.give_normal_plane_axes(major_dir)
  cc_2d = (cc_3d[minor_1],cc_3d[minor_2])
  for p in coords:
#     print("bc_bounds[bound]:",bc_bounds,bound,bc_bounds[bound])
    if np.isclose(p[major_dir],bc_bounds[bound],1e-6):
      result.append((p[minor_1],p[minor_2]))
  
  return result, cc_2d      

# mesh boundary circles of given basecell
# flages is a list with 6 entries, each indicating whether a basecell face should be meshed or not
# bounds order: min_x,min_y,min_z,max_x,max_y,max_z
# @param: basecell center, need this for sorting points in circle order
def mesh_basecell_boundaries(flags,basecell,bounds,cell_center):
  list = []
  # order: {0:"x_left", 1:"y_left", 2:"z_left", 3:"x_right", 4:"y_right", 5:"z_right"}
  # save which boundaries we are going to mesh
  mesh_bound = np.where(flags)[0]
  for b in mesh_bound:
    bound_coords, cc_2d = extract_2d_bc_boundary_coords(basecell.points,b,bounds,cell_center)
    points, cells = mesh_basecell_boundary(bound_coords,b,bounds,cc_2d)
    list.append((points,cells))
    
  return list

# here we only deal with 3 dimensions
def mesh_basecell_boundary(coords_2d,bound,bc_bounds,cell_center):
  assert(len(cell_center) == 2)
  # need this for mapping between planar and space coordinate
  major_dir = bound%3
  minor_dir_1, minor_dir_2 = draw_profile_functions.give_normal_plane_axes(major_dir)
  # sort points in circle order
  coords_2d.sort(key=lambda c:math.atan2(c[0]-cell_center[0], c[1]-cell_center[1]))
  
  info = triangle.MeshInfo()
#   import matplotlib
#   matplotlib.use('tkagg')
#   from matplotlib import pyplot as plt
#   coords_2d = np.asarray(coords_2d)
#   plt.plot(coords_2d[:,0],coords_2d[:,1],'o')
#   plt.show()
  test = [ [elem[0],elem[1]] for elem in coords_2d]
  info.set_points(test)
  info.set_facets(draw_profile_functions.round_trip_connect(0,len(coords_2d)-1))
  
  mesh = triangle.build(info,generate_faces=True)
  mesh_points = np.array(mesh.points)
  mesh_tris = np.array(mesh.elements)
  # check whether we're on the left or right side
  comp = bc_bounds[bound]
  points = []
  cells = []
  # map from 2d point to 3d point
  for p in mesh_points:
    new_p = np.zeros(3)
    new_p[major_dir] = comp
    new_p[minor_dir_1] = p[0]
    new_p[minor_dir_2] = p[1]
    points.append(new_p)
  for tri in mesh_tris:
    cells.append((tri[0], tri[1], tri[2])) 
  
#   import matplotlib
#   matplotlib.use('tkagg')
#   from matplotlib import pyplot as plt
#   plt.gcf()
#   labels = []
#   for i in range(len(points)):
#     labels.append(i)
#   points = np.asanyarray(points)
#   plt.plot(points[:,minor_dir_1],points[:,minor_dir_2],'o')
#   for i, label in enumerate(labels):
#     plt.text(points[i,minor_dir_1],points[i,minor_dir_2],labels[i])
#   plt.triplot(points[:,minor_dir_1],points[:,minor_dir_2],cells)
#   plt.show()
    
  return points,cells

def detect_boundary_edges(polydata):
  featureEdges = vtk.vtkFeatureEdges()
  featureEdges.SetInputData(polydata)
  featureEdges.BoundaryEdgesOn()
  featureEdges.FeatureEdgesOff()
  featureEdges.NonManifoldEdgesOff()
  featureEdges.ManifoldEdgesOff()
  featureEdges.Update()
  
  from vtk.util.numpy_support import vtk_to_numpy
  boundaryPts = vtk_to_numpy(featureEdges.GetOutput().GetPoints().GetData())
  lines = featureEdges.GetOutput().GetLines()
  boundaryEdges = []
  for i in range(lines.GetNumberOfCells()):
    ids = vtk.vtkIdList()
    c = lines.GetNextCell(ids)
    boundaryEdges.append((ids.GetId(0),ids.GetId(1)))

  print("bounds:",boundaryEdges)
  print("num boundary edges:",len(boundaryEdges))
  from itertools import cycle # circular list iterator    
  edgeLoops = []
  loop = []  
  end = -1
  numEdges = len(boundaryEdges)
  handled = [False] * numEdges
  while not all(handled):
    for i,this in cycle(enumerate(boundaryEdges)):
      if handled[i]:
        continue 
      
      if len(loop) == 0:
        loop.append(this)
        end = this[0]
#         print("\nstarting with:",this)
        handled[i] = True
      else:
        next = None
        tail = loop[-1][1]
        # we have at least one edge in loop list, find next neighbor edge
        for j,other in cycle(enumerate(boundaryEdges)):
          if handled[j]:
            continue
#             print("already handled:",boundaryEdges[j])

          if other[0] == tail:
            next = other
            handled[j] = True
#             print("next of ",loop[-1], " is ", next)
            break
            
        assert(next is not None)
        loop.append(next)    
      
#         print(loop)
        if len(loop) > 0 and loop[-1][1] == end:
          edgeLoops.append(loop)
#           print("end: ",loop[-1][1],"=",end)
#           print("handled:",handled)
          loop = []
          end = -1
          break
      
#   print(edgeLoops)
#   for l in edgeLoops:
#     print("\nhole:")
#     for t in l:
#       print(boundaryPts[t[0]])
#     print("len:",len(l))
  return boundaryPts,edgeLoops 

def fill_boundary_loops(points,loops):
  appendPd = vtk.vtkAppendPolyData()
  # each loop consists of a chain of polygons, e.g. (0,1),(1,3),(3,0)
  for i,l in enumerate(loops):
    major = -1 # coordinate component that we can omit
    this = points[l[0][0]]
    next = points[l[0][1]]
    nnext = points[l[1][1]]
    comp = None
    # find out which component we can remove for 2d triangulation
    # assume check first three points is sufficient
    if np.isclose(this[0],next[0]) and np.isclose(this[0],nnext[0]):
      major = 0
      comp = this[0]
    elif np.isclose(this[1],next[1]) and np.isclose(this[1],nnext[1]):
      major = 1
      comp = this[1]
    elif np.isclose(this[2],next[2]) and np.isclose(this[2],nnext[2]):
      major = 2
      comp = this[2]
        
    assert(major > -1)
    minor_1, minor_2 = draw_profile_functions.give_normal_plane_axes(major)
    
#     print("major:",major)
    coords_2d = []
    coords_3d = []
    for t in l:
      p = points[t[0]]
      assert(len(p) == 3)
      coords_2d.append([p[minor_1],p[minor_2]])
      coords_3d.append(p)
    
    assert(len(coords_2d) == len(l))  
    info = triangle.MeshInfo()
    
    test = [ [elem[0],elem[1]] for elem in coords_2d]
    info.set_points(test)
    info.set_facets((draw_profile_functions.round_trip_connect(0,len(test)-1)))
    mesh = triangle.build(info,generate_faces=True)
    mesh_tris = np.array(mesh.elements)
    
    if len(mesh_tris) == 0:
      print("\ncould not mesh hole:")
      print(test)
      print(coords_3d)
      
    appendPd.AddInputData(matviz_vtk.fill_vtk_polydata(coords_3d, mesh_tris))
    appendPd.Update()
  
  print("closed ",len(loops), " holes")  
  return  appendPd.GetOutput()