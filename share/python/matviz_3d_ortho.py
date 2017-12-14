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
def create_3d_interpretation_ortho(args,coords,min_bb,max_bb,s1,s2,s3,scale,samples,grad,nondes=None):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # ip_nx: number of uniform cells in x-direction, can be replaced by csize (size of cell in each direction)
  # grad: type of interpolation ('linear', 'nearest')
  # scale: parameter for scaling the cell size if necessary
  # nondes: (centers, min_bb, max_bb, elem_dim)
  # nondes[centers]: list of elements (corner vertices) that define non-design regions
  
  # writes array with nondes to vtk file with extension *.vtr
#   write_nondes_to_vtr_file(args,min_bb,max_bb,nondes)
  
  # MPI_Init() or MPI_Init_thread() is actually called when you import the MPI
  # use the standard communicator
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  commSize = comm.Get_size()
  
  # point coordinates from h5 file
  centers, _, _ = coords[0:3]
  
  nondes_coords = nondes[0]
  nondes_min = nondes[1]
  nondes_max = nondes[2]
  assert(len(nondes_min) == 3)
  assert(len(nondes_max) == 3) 
  
  # appending cells
  appends = vtk.vtkAppendPolyData()
  
  if scale <= 0:
    scale = 1.0

  # order: min_x,min_y,min_z,max_x,max_y,max_z
  bounds = np.ones(6) * (-1)
  bounds[0:3] = min_bb[0:3]
  bounds[3:6] = max_bb[0:3]
  
  design_bounds = bounds
  
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
  
  assert(len(min_bb) == 3)
  assert(len(max_bb) == 3)

#   print("min_bb:",min_bb)
#   print("nondes_min:",nondes_min)
#   print("max_bb:",max_bb)
#   print("nondes_max:",nondes_max)
  # if any non-design bound exceeds design domain, we need additional ghost layers for non-design  
  if np.any(np.less(np.asarray(nondes_min),np.asarray(min_bb))) or np.any(np.greater(np.asarray(nondes_max),np.asarray(max_bb))):
    diff_min = np.abs(np.asarray(nondes_min) - np.asarray(min_bb))
    diff_max = np.abs(np.asarray(nondes_max) - np.asarray(max_bb))
    
#     print("diff_min:",diff_min)
#     print("diff_max:",diff_max)
    
    # how many ghost cells to add for non-design shell?
    add_left = diff_min / np.asarray(dx,dy,dz)
    add_right = diff_max / np.asarray(dx,dy,dz)
    add_sum = add_left + add_right
    
    add_left = [int(v) for v in add_left]
    add_right = [int(v) for v in add_right]
    
    print("add_left:",add_left)
    print("add_right:",add_right)

  # np.minimum gives elementwise min value
  bounds[0:3] = np.minimum(np.asarray(min_bb),np.asarray(nondes_min))
  bounds[3:6] = np.minimum(np.asarray(max_bb),np.asarray(nondes_max))
  
  print("min_dim:",bounds[0:3])
  print("max_dim:",bounds[4:6])
  
  resolution = (int(samples[0]*args.bc_res+add_sum[0]),int(samples[1]*args.bc_res+add_sum[1]),int(samples[2]*args.bc_res+add_sum[2]))
  print("resolution:",resolution)
  voxel_grid  = np.full(resolution,-1,dtype=int)
  
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  
  hx = width[0] / float(resolution[0])
  hy = width[1] / float(resolution[1])
  hz = width[2] / float(resolution[2])
  
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
    right_back_corner = left_front_corner + np.asarray([dx,dy,dz])
    #print("left_front_corner:",left_front_corner)
    #print("right_back_corner:",right_back_corner)
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
    
    bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="surface_mesh",bc_flags=flags)
#     bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
    bc_input.eta = 0.7
    bc_input.stiffness_as_diameter = True
    cell_obj = basecell.Basecell(bc_input,id,None,(left_front_corner,right_back_corner))
#     cell_obj = basecell.Basecell(bc_input,id,nondes,(left_front_corner,right_back_corner))
#     cell_obj.scale(dx, dy, dz)
#     cell_obj.translate(left_front_corner[0],left_front_corner[1],left_front_corner[2])
#     cell_obj.update()
#     if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
#       continue
#     cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
#     cell_obj.center = cell_center
    voxel_grid[add_left[0]+i*args.bc_res:(i+1)*args.bc_res,add_left[1]+j*args.bc_res:(j+1)*args.bc_res,add_left[2]+k*args.bc_res:(k+1)*args.bc_res] = cell_obj.voxels
    
    basecells.append((cell_obj,(i,j,k)))
    print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)
  
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  
  x = np.arange(bounds[0],bounds[3]+hx,hx)
  y = np.arange(bounds[1],bounds[4]+hy,hy)
  z = np.arange(bounds[2],bounds[5]+hz,hz)
  
  for p in nondes_coords:
    i,j,k = draw_profile_functions.cartesian_to_voxel_coords(p,bounds[0],bounds[1],bounds[2],hx,hy,hz)
    voxel_grid[i,j,k] = -2
  
  from pyevtk.hl import gridToVTK
  gridToVTK("voxels",x,y,z,cellData={"voxels":voxel_grid})
  
  # binary helper array
  helper = np.zeros(voxel_grid.shape[0:3],dtype=int)
  # use voxel info for Marching cubes algorithm
  # set voxels on boundary wit value 0
  # set voxels inside structure with value 1
  # voxels outside structure have value -1
  for i in range(0,resolution[0]):
    for j in range(0,resolution[1]):
      for k in range(0,resolution[2]):
        if voxel_grid[i,j,k] != -1: # valid voxel
          helper[i,j,k] = 1
  
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, normals, values = measure.marching_cubes(helper,spacing=(np.float32(hx),np.float32(hy),np.float32(hz)),allow_degenerate=False)
  verts = np.asarray(verts)
  verts += (hx/2.0,hy/2.0,hz/2.0)
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(bounds[0],bounds[1],bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
  matviz_vtk.show_write_vtk(transformFilter.GetOutput(), 10, "marching.vtp")
  
  connectivity = basecell.getConnectivity(verts,faces)
  points = []
  for i,v in enumerate(verts):
    points.append(draw_profile_functions.Vertex(v,i))
  points = basecell.taubin_smoothing(points,connectivity,30,design_bounds)
  verts = [p.coords for p in points]  
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(bounds[0],bounds[1],bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
  
  matviz_vtk.show_write_vtk(transformFilter.GetOutput(), 10, "smooothed_marching.vtp")
    
  sys.exit()
  
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
 
#     basecells = align_cell_interfaces(nx,ny,nz,basecells,pos)
#     resolve_hanging_nodes(nx,ny,nz,basecells,pos)
      
    import time   
    start = time.time()
    # length of shortest edge
    short = 99999 
    short_e = 99999
    # for appending base cells
    offset = 0 
    vertices = []   
    faces = []
    # appending cells
    appends = vtk.vtkAppendPolyData()
    # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
    for i,cell in enumerate(basecells):
      if (i%100==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
      
#       verts, cells, short_e = add_offset(cell.points, cell.cells, offset)
      
      coords = [ v.coords for v in cell.points]
      pd = matviz_vtk.fill_vtk_polydata(coords, cell.cells)
      appends.AddInputData(pd)
#       vertices.extend(verts)
#       faces.extend(cells)
#       offset = len(vertices)
      
#       for f in faces:
#         assert(type(f[0]) is int)
#         assert(type(f[1]) is int)
#         assert(type(f[2]) is int)
#       if short_e < short:
#         short = short_e
# 
#     assert(short > -1 and short < 99999)
#     points, cells = matviz_vtk.vtk_polydata_to_numpy(cleanFilter.GetOutput())
#     vertices, faces = matviz_vtk.vtk_polydata_to_numpy(pd)
#     print("shortest edge:",short)
#     for i, v in enumerate(vertices):
#       assert(v.idx == i)
    
#     vertices, faces = merge_duplicated_points(vertices,faces,short)
    
#     points, cells = resolve_hanging_nodes(points,cells,(dx,dy,dz),(nx,ny,nz))
    
#     coords = [ v.coords for v in vertices]
#     pd, _ = matviz_vtk.fill_vtk_polydata(coords,faces)
    
    appends.Update()
    
#     print("before cleaning:",pd.GetNumberOfCells())
    # clean coinciden points
    clean = vtk.vtkCleanPolyData()
#     clean.SetInputData(pd)
    clean.SetInputConnection(appends.GetOutputPort())
    clean.Update()
    
#     print("after cleaning:",clean.GetOutput().GetNumberOfCells())

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(clean.GetOutput())
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

# for a given list of basecells, make sure that all points
# at interfaces between two cells coincide
# @param basecells: list of Basecell() objects
# @param pos: list with respective (i,j,k) position in the grid
def align_cell_interfaces(nx,ny,nz,basecells,pos):
  # store base cells in a 3d array to have info about connectivity
  basecell_grid = np.empty((nx,ny,nz), dtype=object)
  
  for i,p in enumerate(pos):
    basecell_grid[p[0],p[1],p[2]] = basecells[i]
  
  new = []
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        print("i,j,k:",i,j,k)
        this = basecell_grid[i,j,k]
        assert(len(this.boundary_points) > 0)
        assert(this is not None)
        if i < nx-1:
#           print("x")
          right = basecell_grid[i+1,j,k]
          this_bp = this.boundary_points[Face_Name.XMAX.value]
          right_bp = right.boundary_points[Face_Name.XMIN.value]
          if len(this_bp) != len(right_bp):
#             print("\nlen(this):",len(this_bp))
            for p in this_bp:
              print(p.coords)
              
            print("\nlen(right):",len(right_bp))
            for p in right_bp:
              print(p.coords)

#           assert(len(this_bp) == len(right_bp))
          if len(this_bp) == len(right_bp) or len(this_bp) > len(right_bp):
            right.replace_boundary_points(this_bp,Face_Name.XMIN.value)
          else:
            this.replace_boundary_points(right_bp,Face_Name.XMAX.value)
        if j < ny-1:
          print("y")
          top = basecell_grid[i,j+1,k]
          this_bp = this.boundary_points[Face_Name.YMAX.value]
          top_bp = top.boundary_points[Face_Name.YMIN.value]
          if len(this_bp) != len(top_bp):
            print("len(this_bp):",len(this_bp)," len(top_bp):",len(top_bp))
            print("\nthis_bp:")
            for p in this_bp:
              print(p)
            print("\ntop_bp:")
            for p in top_bp:
              print(p)  
#           assert(len(this_bp) == len(top_bp))
          if len(this_bp) == len(top_bp) or len(this_bp) > len(top_bp):
            top.replace_boundary_points(this_bp,Face_Name.YMIN.value)
          else:
            this.replace_boundary_points(top_bp,Face_Name.YMAX.value) 
        if k < nz-1:
          print("z")
          front = basecell_grid[i,j,k+1]
          this_bp = this.boundary_points[Face_Name.ZMAX.value]
          front_bp = front.boundary_points[Face_Name.ZMIN.value]
#           assert(len(this_bp) == len(front_bp))
          if len(this_bp) == len(front_bp) or len(this_bp) < len(front_bp):
            front.replace_boundary_points(this_bp,Face_Name.ZMIN.value)
          else:
            this.replace_boundary_points(front_bp,Face_Name.ZMAX.value)
        new.append(this)  
  
  assert(len(new) == len(basecells))
  
  basecells = new
  
  return basecells      

# writes array with nondes to vtk file with extension *.vtr
def write_nondes_to_vtr_file(args,min_bb,max_bb,nondes):
  tmp = args.hom_samples.split(',')
  samples = [int(tmp[0]),int(tmp[1]),int(tmp[2])]
  resolution = np.asarray([args.bc_res*samples[0],args.bc_res*samples[1],args.bc_res*samples[2]])
  
  bounds = np.ones(6) * (-1)
  bounds[0:3] = min_bb[0:3]
  bounds[3:6] = max_bb[0:3]
  
  nondes_grid = np.full(resolution,0,dtype=int)
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  
  hx = width[0] / float(resolution[0])
  hy = width[1] / float(resolution[1])
  hz = width[2] / float(resolution[2])
  
  for p in nondes:
    i,j,k = draw_profile_functions.cartesian_to_voxel_coords(p,bounds[0],bounds[1],bounds[2],hx,hy,hz)
    nondes_grid[i,j,k] = 1
    
  x = np.arange(bounds[0],bounds[3]+hx,hx)
  y = np.arange(bounds[1],bounds[4]+hy,hy)
  z = np.arange(bounds[2],bounds[5]+hz,hz)
  
  from pyevtk.hl import gridToVTK
  gridToVTK("nondesign",x,y,z,cellData={"nondes":nondes_grid})
