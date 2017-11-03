#1/usr/bin/env python

import basecell
import draw_profile_functions
import numpy as np
import matviz_vtk
import pymp
import vtk
import math
import sys

try:
  import meshpy.triangle as triangle
#   import meshpy.tet as tet
except:
  print("Failed to load meshpy - need it for tetrahedralizing basecell mesh")

try:
  from mpi4py import MPI
except:
  print("WARNING: Could not load mpi4py!")


# similar to create_3d_cross_ip; # without rotation and shearing
# returns

def create_3d_interpretation_ortho(args,coords,min_bb,max_bb,s1,s2,s3,scale,samples,grad,thresh):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # thres: threshold value for design variables s1/s2/s3. The cell is not visualized if s1,s2,s3 <= thres
  
  # MPI_Init() or MPI_Init_thread() is actually called when you import the MPI
  # use the standard communicator
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  commSize = comm.Get_size()
  
  # point coordinates from h5 file
  centers, _, _ = coords[0:3]
  
  # appending cells
  appends = vtk.vtkAppendPolyData()
  
  if scale <= 0:
    scale = 1.0

  # assume we always start at (0,0,0)
  # order: min_x,min_y,min_z,max_x,max_y,max_z
  bounds = np.ones(6) * (-1)
  bounds[0:3] = min_bb[0:3]
  bounds[3:6] = max_bb[0:3]
  
  # 0:xmin,1:ymin,2:zmin,3:xmax,4:ymax,5:zmax  
  delta = (abs(bounds[3] - bounds[0]), abs(bounds[4] - bounds[1]), abs(bounds[5] - bounds[2]))
  
  # set size dx/dy/dz of one cell
  dx = delta[0] / samples[0]
  dy = delta[1] / samples[1]
  dz = delta[2] / samples[2]
  
  min_thresh = 0.1
  max_thresh = 0.82
    
  # where we want nodes
  nx = int(delta[0] / dx)
  ny = int(delta[1] / dy)
  nz = int(delta[2] / dz)

#   print("samples:",samples)  
#   print("delta:",delta)
  print("bounds:",bounds)  
  print("dx,dy,dz:",dx,dy,dz)
#   print("nx,ny,nz:",nx,ny,nz) 
  
  data_grid, data_grid_near, sample_coords= matviz_vtk.get_interpolation_row_major(coords, bounds, grad, s1, s2, s3, nx, ny, nz, dx, dy, dz)
  
  basecells = list()
  boundary_flags = list()
  
  nProblem = nx*ny*nz

  # distribute loop evenly over number of processes and take care
  # of remainder -> difference between work chunks can be 1
  # e.g. 10 loop runs, 4 processes, number of runs per process: [3,3,2,2] 
  # number of loop runs per process:
  nRuns = [int(nProblem/commSize) + (1 if p < nProblem%commSize else 0) for p in range(0,commSize)]
  start = int(np.sum(nRuns[0:rank]))
  end = int(start + nRuns[rank])
  
  for id in range(start,end):
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
#     flags[0] = True if i == 0 else False
#     flags[1] = True if j == 0 else False
#     flags[2] = True if k == 0 else False
#     flags[3] = True if i == nx-1 else False
#     flags[4] = True if j == ny-1 else False
#     flags[5] = True if k == nz-1 else False
    for c in range(len(flags)):      
      flags[c] = True if grid_coords[c%3] == grid_bounds[c] else False 
    cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
    cell_obj.center = cell_center
    boundary_list = []
    # at least one boundary circle needs to be triangulated
    if any(flags):
      boundary_list = mesh_basecell_boundaries(flags,cell_obj,bounds,cell_center)

    basecells.append((cell_obj,boundary_list,flags))
    print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)

  # broadcast all data to master and exit script
  if rank != 0:
    sys.stdout.flush()
    result = (basecells,boundary_flags)
    comm.send(result,dest=0,tag=1)
    status = MPI.Status()
    tmp = comm.recv(source=0,status=status)
    # make sure communication wors
    if status.Get_tag() == 999:
      print("\n             rank ",rank," exiting now")
      sys.stdout.flush()
      sys.exit()
  else:
    finished_workers = 0
    # a status object containing source, tag and size of a message
    status = MPI.Status()
    while finished_workers != commSize-1:
      data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
      # we store the result
      tmp_bc, tmp_bf = data
      basecells.extend(tmp_bc)
      boundary_flags.extend(tmp_bf)
      #get the message tag
      tag = status.Get_tag()
      source = status.Get_source()
      comm.send(None,dest=source,tag=999)
      finished_workers += 1
    
    print("finished_workers ",finished_workers," commsize: ",commSize)
    sys.stdout.flush()
    
    import time   
    flags = []
    start = time.time()
    # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
    for i,obj in enumerate(basecells):
      if (i%100==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
        
      cell = obj[0]
      flags.append(obj[2])
      vtk_points = vtk.vtkPoints()
      pd = matviz_vtk.fill_vtk_polydata(cell.points, cell.cells,cell.center)
      
      appends.AddInputData(pd)
      # meshed boundary circles
      circles = obj[1]
  #     assert(len(circles) > 0)
      # list with meshes on the boundaries
      for l in circles:
        points = np.asarray(l[0])
        pd = matviz_vtk.fill_vtk_polydata(l[0], l[1])
          
        appends.AddInputData(pd)
  
    appends.Update() # not sure if we have to do this in each loop iteration
    
    # merge duplicated points etc.
    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appends.GetOutputPort())
    cleanFilter.Update()
    
    boundaryPts,loops = detect_boundary_edges(cleanFilter.GetOutput())
    appends.AddInputData(fill_boundary_loops(boundaryPts,loops))
     
    cleanFilter.Update()
    pd = cleanFilter.GetOutput()
    
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(pd)
    normals.SetConsistency(1)
    normals.SetAutoOrientNormals(1)
    normals.Update()
  
    return normals.GetOutput()
#   return cleanFilter.GetOutput()
    
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
  import matplotlib
  from matplotlib import pyplot as plt
  coords_2d = np.asarray(coords_2d)
#   plt.plot(coords_2d[:,0],coords_2d[:,1],'o')
#   plt.show()
  test = [ [np.float64(elem[0]),np.float64(elem[1])] for elem in coords_2d]
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

#   print("boundaryEdges:",boundaryEdges)
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
    nnnext = points[l[2][1]]
    comp = None
    eps = 1e-12
    flag_x = True
    flag_y = True
    flag_z = True
    # find out which component we can remove for 2d triangulation
    # assume check first three points is sufficient
    for id in range(len(l)):
      next = points[l[id][1]]
      if not np.isclose(this[0],next[0],eps):
        flag_x = False
      if not np.isclose(this[1],next[1],eps):
        flag_y = False 
      if not np.isclose(this[2],next[2],eps):
        flag_z = False
        
    if flag_x:
      major = 0
      comp = this[0]
    elif flag_y:
      major = 1
      comp = this[1]
    else:
      assert(flag_z)
      major = 2
      comp = this[2]
        
    assert(major > -1)
    minor_1, minor_2 = draw_profile_functions.give_normal_plane_axes(major)
    
#     print("major:",major)
    coords_2d = []
    for t in l:
      p = points[t[0]]
      assert(len(p) == 3)
      coords_2d.append([p[minor_1],p[minor_2]])
    
    assert(len(coords_2d) == len(l))  
    info = triangle.MeshInfo()
    
    test = [ [np.float64(elem[0]),np.float64(elem[1])] for elem in coords_2d]
    info.set_points(test)
    info.set_facets((draw_profile_functions.round_trip_connect(0,len(test)-1)))
    mesh = triangle.build(info,generate_faces=True)
    mesh_points = np.array(mesh.points)
    mesh_tris = np.array(mesh.elements)
    
#     import matplotlib
#     from matplotlib import pyplot as plt
#     matplotlib.use('tagg')
#     plt.plot(mesh_points[:,0],mesh_points[:,1],'o')
#     plt.triplot(mesh_points[:,0],mesh_points[:,1],mesh_tris)
#     plt.show()
    
    if len(mesh_tris) == 0:
      print("\ncould not mesh hole:")
      print(test)
      sys.exit()
    
    new_points = []
    # map from 2d point to 3d point
    for p in mesh_points:
      new_p = np.zeros(3)
      new_p[major] = comp
      new_p[minor_1] = p[0]
      new_p[minor_2] = p[1]
      new_points.append(new_p)
      
    appendPd.AddInputData(matviz_vtk.fill_vtk_polydata(new_points, mesh_tris))
    appendPd.Update()
  
  print("closed ",len(loops), " holes")  
  return  appendPd.GetOutput()
    