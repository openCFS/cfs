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
  from draw_profile_functions import Face_Name
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
#   print("bounds:",bounds)  
#   print("dx,dy,dz:",dx,dy,dz)
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
    cell_obj = basecell.Basecell(bc_input,id)
    cell_obj.scale(dx, dy, dz)
    cell_obj.translate(left_front_corner[0],left_front_corner[1],left_front_corner[2])
    cell_obj.update()
    if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
      continue
    cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
    cell_obj.center = cell_center
    
    basecells.append((cell_obj,(i,j,k)))
    print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)

  # broadcast all data to master and exit script
  if rank != 0:
    sys.stdout.flush()
    result = basecells
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
    pos = [c[1] for c in basecells ]
    basecells = [c[0] for c in basecells]
#     basecells = [None] * nx*ny*nz
    while finished_workers != commSize-1:
      data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
      cells = [d[0] for d in data ]
      p = [d[1] for d in data ]
      # we store the result
      basecells.extend(cells)
      pos.extend(p)
      #get the message tag
      tag = status.Get_tag()
      source = status.Get_source()
      comm.send(None,dest=source,tag=999)
      finished_workers += 1
    
    print("finished_workers ",finished_workers," commsize: ",commSize)
    sys.stdout.flush()

    # store base cells in a 3d array to have info about connectivity
    basecell_grid = np.empty((nx,ny,nz), dtype=object)
     
    for i,p in enumerate(pos):
      basecell_grid[p[0],p[1],p[2]] = basecells[i]
      
    for i in range(nx):
      for j in range(ny):
        for k in range(nz):
          this = basecell_grid[i,j,k]
          assert(len(this.boundary_points) > 0)
          assert(this is not None)
          if i < nx-1:
            this_bp = [len(l) for l in this.boundary_points]
            right = basecell_grid[i+1,j,k]
            right_bp = [len(l) for l in right.boundary_points]
            assert(len(this.boundary_points[Face_Name.XMAX.value]) == len(right.boundary_points[Face_Name.XMIN.value]))
            right.replace_boundary_points(this.boundary_points[Face_Name.XMAX.value],Face_Name.XMIN.value)
#             right.boundary_points = this.boundary_points
          if j < ny-1:  
            top = basecell_grid[i,j+1,k]
            assert(len(this.boundary_points[Face_Name.YMAX.value]) == len(top.boundary_points[Face_Name.YMIN.value]))
#             top.boundary_points = this.boundary_points
            top.replace_boundary_points(this.boundary_points[Face_Name.YMAX.value],Face_Name.YMIN.value)
          if k < nz-1:
            front = basecell_grid[i,j,k+1]
            assert(len(this.boundary_points[Face_Name.ZMAX.value]) == len(front.boundary_points[Face_Name.ZMIN.value]))
            front.replace_boundary_points(this.boundary_points[Face_Name.ZMAX.value],Face_Name.ZMIN.value)
            front.boundary_points = this.boundary_points 
      
    import time   
    start = time.time()
    # length of shortest edge
    short = 99999 
    short_e = 99999
    # for appending base cells
    offset = 0 
    vertices = []   
    faces = []
    # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
    for i,cell in enumerate(basecells):
      if (i%100==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
      
      verts, cells, short_e = add_offset(cell.points, cell.cells, offset)
      
      vertices.extend(verts)
      faces.extend(cells)
      offset = len(vertices)
      
      for f in faces:
        assert(type(f[0]) is int)
        assert(type(f[1]) is int)
        assert(type(f[2]) is int)
      if short_e < short:
        short = short_e

    assert(short > -1 and short < 99999)
#     points, cells = matviz_vtk.vtk_polydata_to_numpy(cleanFilter.GetOutput())
    coords = [ v.coords for v in vertices]
#     points, cells = matviz_vtk.vtk_polydata_to_numpy(coords,faces)
    print("shortest edge:",short)
    for i, v in enumerate(vertices):
      assert(v.idx == i)
    
#     vertices, faces = merge_duplicated_points(vertices,faces,short)
    
#     points, cells = resolve_hanging_nodes(points,cells,(dx,dy,dz),(nx,ny,nz))
    
#     boundaryPts,loops = detect_boundary_edges(cleanFilter.GetOutput())
#     appends.AddInputData(fill_boundary_loops(boundaryPts,loops))

#     cleanFilter.Update()
#     pd = cleanFilter.GetOutput()
    pd, _ = matviz_vtk.fill_vtk_polydata(coords,faces)
    
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
def merge_duplicated_points(verts,faces,short):
  from scipy import spatial
  tol = 0.5 * short
  coords = [v.coords for v in verts]
  tree = spatial.KDTree(coords)
  #find all pairs of points within a distance
  # gives indices with tuples of points indices
  idx = tree.query_pairs(tol)
#   print("duplicated ids:",idx)
  # stores new point ids
  # list idx gives point id that should be replaced with value of element idx
  # max id of first tuple entry
  replace = [-1] * len(verts)
#   print("\nnumber of duplicated points:",len(idx))
  for i in idx:
#     print("pair:",points[i[0]],points[i[1]])
    replace[i[0]]=i[1]
#     replace.append((i[0],i[1]))
  
  
  return rearrange_points_cells(verts,faces,replace)

# @param replace: list with 2-element-tuples storing info which id should be replaced with which one
# returns points and cells list with 0-based and consecutive ids 
def rearrange_points_cells(verts,faces,replace):   
  # create new points list with ordered and unique point ids (no gaps)
  # store old valid ids                                    
  old_ids = set()
  
  # for each cell
  for i in range(len(faces)):
    new_face = [None] * 3
    for j in range(3):
      if replace[faces[i][j]] > -1:
        # if there is something to correct
        new_face[j] = replace[faces[i][j]]
#         print("replace ",faces[i][j]," with ",replace[faces[i][j]])
#         print("replace ",verts[faces[i][j]].coords," with ",verts[replace[faces[i][j]]].coords)
      else:
        new_face[j] = faces[i][j]   
      old_ids.add(new_face[j])
      
    faces[i] = new_face  
  
  # set does not support indexing
#   old_ids = list(old_ids)
#   # create map to connect old valids ids with new ones (0-based)
#   map = len(verts) * [-1]
#   for i in range(len(old_ids)):
#     map[old_ids[i]] = i
#      
#   # now correct node ids
#   new_verts = len(old_ids) * [None]
#   for oi in old_ids:
#     new_verts[map[oi]] = draw_profile_functions.Vertex(verts[oi].coords,map[oi])
# #     print("verts[oi].idx:",verts[oi].idx," oi:",oi)
#     assert(verts[oi].idx == oi)
#    
#   for n in new_verts:
#     assert(n is not None)
#   
#   # correct element vertices
#   new_faces = []
#   for f in faces:
#     new_faces.append((map[f[0]],map[f[1]],map[f[2]]))
    
  return verts,faces  
#   return new_verts, new_faces

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
  
  # store all neighbors for all points
  coords = [p.coords for p in points]
  connectivity = basecell.getConnectivity(coords, cells)
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
      p = np.asarray(vertex.coords)
      # get all connected neighbors
      neighbors = connectivity[vertex.idx]
      # for each neighbor, form an edge and compare
      # if any vertex in the surroundings has distance 0 (hanging node)
      for neighbor in neighbors:
#         print("\nneighbor of ",p1," is ",neighbor,points[neighbor])
        # edge is defined by p1 and p2
        p2 = np.asarray(points[neighbor])
        for other in list:
          # prevent self-comparison
          if other.idx == neighbor or other.idx == vertex.idx:
            continue
          
          p = np.asarray(other.coords)
          if np.linalg.norm(p2-p) > np.max(dx):
            continue
          
          # check if projection of point A lies between 2 points P and S
          # defining the edge
          
          
          d = np.linalg.norm(np.cross(p-p1, p-p2))
          d /= np.linalg.norm(p2-p1)
          
          if np.isclose(d,0,1e-6):
            print("\ndetect hanging node: edge=",p1,p2," point:",p)
            print("ids:",vertex.idx,neighbor,other.idx)
            count += 1
  
  print("found ",count," hanging nodes")
#   # detect hanging nodes
#   for list in lists:
#     # each list contains triangle vertices for one cell interface
            
  return points, cells

# @param verts: list of objects of type Vertex
# @param cells: list of triangles
# @param offset: value to shift vertex ids
def add_offset(verts,cells,offset):
  assert(offset >=0)
  coords = [None] * len(verts)
  # length of shortest edge 
  short = 999999
  # add offset to all triangle vertices
  for i in range(len(verts)):
    assert(type(verts[i]) is draw_profile_functions.Vertex)
    verts[i].idx += offset
  
  for i in range(len(cells)):
    cell = cells[i]
    e01 = np.linalg.norm(np.asarray(verts[cell[1]].coords) - np.asarray(verts[cell[0]].coords))
    e12 = np.linalg.norm(np.asarray(verts[cell[2]].coords) - np.asarray(verts[cell[1]].coords))
    e02 = np.linalg.norm(np.asarray(verts[cell[2]].coords) - np.asarray(verts[cell[0]].coords))
    if np.min([e01,e12,e02]) < short:
      short = np.min([e01,e12,e02])
      
    cells[i]= (int(cells[i][0])+offset,int(cells[i][1]+offset),int(cells[i][2]+offset))  
  
  return verts, cells, short
