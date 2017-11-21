#!/usr/bin/env python


import numpy as np
import matviz_vtk
import vtk
import math
import sys

try:
  import meshpy.triangle as triangle
except:
  print("Warning: Failed to load meshpy - need it for basecell surface mesh")

try:
  from mpi4py import MPI
except:
  print("Warning: Could not load mpi4py!")

import basecell
  
try:
  import basecell
  import draw_profile_functions
except:
  print("Warning: Could not load basecell and draw_profile_functions!")
  
    
# similar to create_3d_cross_ip; # without rotation and shearing
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
  
  # thread local data, each entry is a tuple with Basecell object and its boundary circle meshes
  tl_data = list()
  basecells = []
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
    
    # flags for meshing circles on the boundary
    flags = [None] * 6
#     grid_bounds = [0,0,0,nx-1,ny-1,nz-1]
    grid_bounds = [0,nx-1,0,ny-1,0,nz-1]
    grid_coords = [i,j,k]
#     flags[0] = True if i == 0 else False
#     flags[1] = True if j == 0 else False
#     flags[2] = True if k == 0 else False
#     flags[3] = True if i == nx-1 else False
#     flags[4] = True if j == ny-1 else False
#     flags[5] = True if k == nz-1 else False
    for c in range(len(flags)):      
      flags[c] = True if grid_coords[c%3] == grid_bounds[c] else False 
    
    bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="surface_mesh",bc_flags=flags)
    bc_input.eta = 0.7
    bc_input.stiffness_as_diameter = True
    cell_obj = basecell.Basecell(bc_input,id,(i,j,k))
    cell_obj.scale(dx, dy, dz)
    cell_obj.translate(left_front_corner[0],left_front_corner[1],left_front_corner[2])
    cell_obj.update()
    if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
      continue
    cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
    cell_obj.center = cell_center
    
    basecells.append(cell_obj)
    print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)

  # broadcast all data to master and exit script
  if rank != 0:
    sys.stdout.flush()
    result = local_basecells
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
#     basecells = [None] * nx*ny*nz
    while finished_workers != commSize-1:
      data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
      print("data:",data)
      # we store the result
      basecells.extend(data)
      #get the message tag
      tag = status.Get_tag()
      source = status.Get_source()
      comm.send(None,dest=source,tag=999)
      finished_workers += 1
    
    print("finished_workers ",finished_workers," commsize: ",commSize)
    sys.stdout.flush()
    
    import time   
    start = time.time()
    # length of shortest edge
    short = 99999 
    short_e = 99999   
    print(basecells)
    # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
    for i,cell in enumerate(basecells):
      if (i%100==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
        
      vtk_points = vtk.vtkPoints()
      coords = [p.coords for p in cell.points]
      pd, short_e = matviz_vtk.fill_vtk_polydata(coords, cell.cells)
      if short_e < short:
        short = short_e
      appends.AddInputData(pd)
      # meshed boundary circles
  
    appends.Update() # not sure if we have to do this in each loop iteration
    
    # merge duplicated points etc.
#     cleanFilter = vtk.vtkCleanPolyData()
#     cleanFilter.SetInputConnection(appends.GetOutputPort())
#     cleanFilter.Update()
    
    assert(short > -1 and short < 99999)
#     points, cells = matviz_vtk.vtk_polydata_to_numpy(cleanFilter.GetOutput())
    points, cells = matviz_vtk.vtk_polydata_to_numpy(appends.GetOutput())
    points, cells = merge_duplicated_points(points,cells,short)
    
#     points, cells = resolve_hanging_nodes(points,cells,(dx,dy,dz),(nx,ny,nz))
    
#     boundaryPts,loops = detect_boundary_edges(cleanFilter.GetOutput())
#     appends.AddInputData(fill_boundary_loops(boundaryPts,loops))

#     cleanFilter.Update()
#     pd = cleanFilter.GetOutput()
    pd, _ = matviz_vtk.fill_vtk_polydata(points,cells)
    
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
    
    appendPd.AddInputData(matviz_vtk.fill_vtk_polydata(new_points, mesh_tris)[0])
    appendPd.Update()
    
  print("closed ",len(loops), " holes")  
  return  appendPd.GetOutput()

# merge duplicated points (with given tolerance)
# at interface of two base cells
def merge_duplicated_points(points,cells,short):
  from scipy import spatial
  tol = 0.8 * short
  tree = spatial.KDTree(points)
  #find all pairs of points within a distance
  # gives indices with tuples of points indices
  idx = tree.query_pairs(tol)
  
  # stores tuples of point ids
  # first tuple entry gives point id that should be replaced with second tuple entry
  replace = []
#   print("\nnumber of duplicated points:",len(idx))
  new_points = []
  for i in idx:
#     print("pair:",points[i[0]],points[i[1]])
    replace.append((i[0],i[1]))
  
  
  return rearrange_points_cells(points,cells,replace)

# @param replace: list with 2-element-tuples storing info which id should be replaced with which one
# returns points and cells list with 0-based and consecutive ids 
def rearrange_points_cells(points,cells,replace):   
  # create new points list with ordered and unique point ids (no gaps)
  # store old valid ids                                    
  old_ids = set()
  
  print("found ", len(replace), " duplicates")
  # for each cell
  for i in range(len(cells)):
    cell = cells[i]
    assert(len(cell) == 3)
    # for each replacement
    for r in replace:
      assert(len(r) == 2)
      # for each cell vertex
      for i in range(len(cell)):
        if cell[i] == r[0]:
#           print("replaced id ", cell[i], " with ", r[1]) 
          cell[i] = r[1]
           
    old_ids.add(cell[0])    
    old_ids.add(cell[1])
    old_ids.add(cell[2])
  
  # set does not support indexing
  old_ids = list(old_ids)
  # create map to connect old valids ids with new ones (0-based)
  map = len(points) * [-1]
  for i in range(len(old_ids)):
    map[old_ids[i]] = i
    
  # now correct node ids
  newpoints = len(old_ids) * [None]
  for oi in old_ids:
    newpoints[map[oi]] = points[oi]
  
  for n in newpoints:
    assert(n is not None)
  
  # correct element vertices
  newcells = []
  for i in range(len(cells)):
    newcells.append((map[cells[i][0]],map[cells[i][1]],map[cells[i][2]]))
    
  return newpoints, newcells

def get_interface_points(points,cells,dx,nx):
  from scipy import spatial
  # kd tree for efficient neighbor search
  tree = spatial.KDTree(points)

  r = np.max([dx[0],dx[1],dx[2]])
  tol = 3 * [0.0]
  tol[0] = dx[0] /20.0
  tol[1] = dx[1] /20.0
  tol[2] = dx[2] /20.0
  
#   print("tol x,y,z", tol)  
  
  lists = []
  # store a point list for each cell interface
  # each list entry is a tuple with (id,coords)
  check_lists = [] 
  for i in range(1,nx[0]+1):
    for j in range(1,nx[1]+1):
      for k in range(1,nx[2]+1):

        print("i,j,k:",i,j,k)
        
        # midpoint of an interface
        midpoint = None
        new_points = []
        for count in range(0,3):
          if count == 0: #x
            major = i
            midpoint = [i*dx[0],j*dx[1]/2,k*dx[2]/2]
          elif count == 1: #y
            major = j
            midpoint = [i*dx[0]/2,j*dx[1],k*dx[2]/2]
          else: #z
            major = k
            midpoint = [i*dx[0]/2,j*dx[1]/2,k*dx[2]]
            
          # neglect boundary circles  
          if major >= nx[count]:
            continue
          
          assert(midpoint is not None)
          ids = tree.query_ball_point(midpoint,r)
          assert(len(ids) > 0)
          cands = [points[id] for id in ids]
          
          for idx,cand in enumerate(cands):
            if midpoint[count]-tol[count] <= cand[count] <= midpoint[count]+tol[count]:
              # store id and coords of vertex
              new_points.append(draw_profile_functions.Vertex(cand,ids[idx]))   
          
          assert(new_points is not None)
          assert(len(new_points) > 0)    
          lists.append(new_points)    

  print("number of interfaces:",len(lists))
#   for l in lists:
#     print("\nnew_points:\n",l)
#   sys.exit()
    
  return lists      
# detect hanging node and resolve them by adding additional triangles;
# look only at points and triangles at interfaces between two base cells:
# a hanging node is a vertex with distance ~0 to another triangle edge
# distance from a vertex P to edge AB:
# d = 2 * triangle area ABC / length of vector AB
# @param dx: three lattice spacings dx,dy,dz in a tuple/list
# @param nx: number of cells in 3 directions nx,ny,nz in a tuple/list
def resolve_hanging_nodes(points,cells,dx,nx):
  assert(dx[0] > 0)
  assert(dx[1] > 0)
  assert(dx[2] > 0)
  
  # each entry is a list
  # each list has entries of type Vertes()
  lists = get_interface_points(points,cells,dx,nx)
  # store all neighbors for all points
  connectivity = basecell.getConnectivity(points, cells)
#   print("connectivity of node 4871=",points[4871],": ")
#   for id in connectivity[4871]:
#     print(id,points[id])
#   print("connectivity of node 2420=",points[2420],": ")
#   for id in connectivity[2420]:
#     print(id,points[id])
#   sys.exit()
  assert(len(connectivity) == len(points))
  
  count = 0
  for l,list in enumerate(lists):
    print("interface ",l)
    for vertex in list:
      p1 = np.asarray(vertex.coords)
      # get all connected neighbors
      neighbors = connectivity[vertex.id]
      # for each neighbor, form an edge and compare
      # if any vertex in the surroundings has distance 0 (hanging node)
      for neighbor in neighbors:
#         print("\nneighbor of ",p1," is ",neighbor,points[neighbor])
        # edge is defined by p1 and p2
        p2 = np.asarray(points[neighbor])
        for other in list:
          # prevent self-comparison
          if other.id == neighbor or other.id == vertex.id:
            continue
          
          p = np.asarray(other.coords)
          if np.linalg.norm(p2-p) > np.max(dx):
            continue
          # vertices p1 and p2, point to test p
          # distance between P and edge v1v2
          # d = 2 * triangle area ABC / length of vector AB
          d = np.linalg.norm(np.cross(p-p1, p-p2))
          d /= np.linalg.norm(p2-p1)
          
          if np.isclose(d,0,1e-6):
            print("\ndetect hanging node: edge=",p1,p2," point:",p)
            print("ids:",vertex.id,neighbor,other.id)
            count += 1
  
  print("found ",count," hanging nodes")
#   # detect hanging nodes
#   for list in lists:
#     # each list contains triangle vertices for one cell interface
            
  return points, cells
