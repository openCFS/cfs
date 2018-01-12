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

try:
  import basecell
  import draw_profile_functions
except:
  print("Warning: Could not load basecell and draw_profile_functions!")
  
    
# similar to create_3d_cross_ip; # without rotation and shearing
def create_3d_interpretation_ortho(args,coords,min_bb,max_bb,s1,s2,s3,scale,samples,grad,nondes=None):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ... 
  # coords, s1, s2, s3, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # min_bb/max_bb: bounding box of design regions
  # scale: parameter for scaling the cell size if necessary
  # nondes: (centers, min_bb, max_bb, elem_dim)
  # nondes[centers]: list of elements (corner vertices) that define non-design regions
  
  # writes array with nondes to vtk file with extension *.vtr
#   write_nondes_to_vtr_file(args,min_bb,max_bb,nondes)
  
  # MPI_Init() or MPI_Init_thread() is actually called when you import the MPI
  # use the standard communicator
  comm = MPI.COMM_WORLD
  
  # point coordinates from h5 file
  centers, _, _ = coords[0:3]
  
  # appending cells
  appends = vtk.vtkAppendPolyData()
  
  if scale <= 0:
    scale = 1.0

  # order: min_x,min_y,min_z,max_x,max_y,max_z
  design_bounds = np.ones(6) * (-1)
  design_bounds[0:3] = min_bb[0:3]
  design_bounds[3:6] = max_bb[0:3]
  
  # 0:xmin,1:ymin,2:zmin,3:xmax,4:ymax,5:zmax  
  h_des = (abs(design_bounds[3] - design_bounds[0]), abs(design_bounds[4] - design_bounds[1]), abs(design_bounds[5] - design_bounds[2]))
  
  # set size dx/dy/dz of one cell
  dx_des = h_des[0] / samples[0]
  dy_des = h_des[1] / samples[1]
  dz_des = h_des[2] / samples[2]
  
  min_thresh = 0.1
  max_thresh = 0.82
    
#   print("samples:",samples)  
#   print("h_des:",h_des)
#   print("design_bounds:",design_bounds)  
#   print("dx_des,dy_des,dz_des:",dx_des,dy_des,dz_des)
#   print("nx,ny,nz:",nx,ny,nz) 
  
  data_grid, data_grid_near, sample_coords= matviz_vtk.get_interpolation_row_major(coords, design_bounds, grad, s1, s2, s3, samples[0], samples[1], samples[2], dx_des, dy_des, dz_des)
  
  ################# voxel world #########################################
  nondes_coords = nondes[0]
  nondes_min = nondes[1]
  nondes_max = nondes[2]
  assert(len(nondes_min) == 3)
  assert(len(nondes_max) == 3) 
  # np.minimum gives elementwise min value
  bounds = [None] * 6
  bounds[0:3] = np.minimum(np.asarray(min_bb),np.asarray(nondes_min))
  bounds[3:6] = np.maximum(np.asarray(max_bb),np.asarray(nondes_max))
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  print("\nmin_bb design:",min_bb)
  print("max_bb design:",max_bb)
  print("min:",bounds[0:3])
  print("max:",bounds[3:6])
  
  my_grid = DistributedGrid(comm,samples,args.bc_res,bounds)
  my_grid.to_info()
  
#   nCells = nx*ny*nz
#   # distribute loop evenly over number of processes and take care
#   # of remainder -> difference between work chunks can be 1
#   # e.g. 10 loop runs, 4 processes, number of runs per process: [3,3,2,2]
#   # number of loop runs per process:
#   nRuns = [int(nCells/commSize) + (1 if p < nCells%commSize else 0) for p in range(0,commSize)]
#   start = int(np.sum(nRuns[0:rank]))
#   end = int(start + nRuns[my_grid.rank])
  
#   if my_grid.rank == 0: 
#     eps = 1e-6
#     x = np.arange(0,resolution[0]+1,1)
#     y = np.arange(0,resolution[1]+1,1)
#     z = np.arange(0,resolution[2]+1,1)
#     decomp = np.full(resolution,-1,dtype=int)
#     decomp[start:end+1,:,:] = rank
#     sys.stdout.flush()
#     finished_workers = 0
#     from pyevtk.hl import gridToVTK  
#     while finished_workers != commSize-1:
#       # a status object containing source, tag and size of a message
#       status = MPI.Status()
#       data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
#       
#       assert(data[0] >= 0 and data[1] <= resolution[0])
#       #get the message tag
#       tag = status.Get_tag()
#       source = status.Get_source()
#       decomp[data[0]:data[1],:,:] = source
#       comm.send(None,dest=source,tag=999)
#       
#       print(" source:",source, " data:",data)
#       sys.stdout.flush()
#       finished_workers += 1
#     
#     gridToVTK("decomp",x,y,z,cellData={"decomp":decomp})
#   else:
#     result = (start,end)
#     comm.send(result,dest=0,tag=1)
#     status = MPI.Status()
#     tmp = comm.recv(source=0,status=status)
#     # make sure communication wors
#     if status.Get_tag() == 999:
#       print("\n             rank ",rank," exiting now")
#       sys.stdout.flush()
#       sys.exit()
#       
#   sys.exit()
  
  local_id = 0
  for i in range(my_grid.start_x,my_grid.end_x):
    for j in range(samples[1]):
      for k in range(samples[2]):
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
        
        flags = None
        bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
        bc_input.eta = 0.7
        bc_input.stiffness_as_diameter = True
        cell_obj = basecell.Basecell(bc_input,id)
        # local i,j,k
        li,lj,lk = get_3d_grid_coords(local_id, my_grid.chunks, samples[1], samples[2])
        print("rank:",my_grid.rank," global i,j,k:",i,j,k," x1,x2,y1,y2,z1,z2:",x1,x2,y1,y2,z1,z2)
        my_grid.grid.data[li*args.bc_res:(li+1)*args.bc_res,lj*args.bc_res:(lj+1)*args.bc_res,lk*args.bc_res:(lk+1)*args.bc_res] = cell_obj.voxels
    
        local_id += 1
    
  eps = 1e-6
  
  borders = my_grid.communicate_edges()
  print("len(borders)",len(borders))

#   void_elems = []
#   # read void non-design manually  
#   import csv
#   with open('holes.txt') as csvfile:
#     readCSV = csv.reader(csvfile, delimiter=',')  
#     for row in readCSV:
#       # 4 elems with 3 coord components
#       assert(len(row) == 12)
#       void_elems.append([(float(row[0]),float(row[1]),float(row[2])), (float(row[3]),float(row[4]),float(row[5])), (float(row[6]),float(row[7]),float(row[8])), (float(row[9]),float(row[10]),float(row[11]))])
#     
#     assert(len(void_elems) > 0)
    
  x = np.arange(my_grid.bounds[0],my_grid.bounds[3]+my_grid.grid.hx-eps,my_grid.grid.hx)
  y = np.arange(my_grid.bounds[1],my_grid.bounds[4]+my_grid.grid.hy-eps,my_grid.grid.hy)
  z = np.arange(my_grid.bounds[2],my_grid.bounds[5]+my_grid.grid.hz-eps,my_grid.grid.hz)
  
  draw_non_design(nondes_coords, my_grid.grid.data, my_grid.bounds[0:3], (my_grid.grid.hx,my_grid.grid.hy,my_grid.grid.hz), 0)
#   draw_non_design(void_elems, my_grid.grid.data, my_grid.bounds[0:3], (my_grid.grid.hx,my_grid.grid.hy,my_grid.grid.hz), -1)    
  from pyevtk.hl import gridToVTK  
  gridToVTK("voxels"+str(my_grid.rank),x,y,z,cellData={"voxels":my_grid.grid.data})
  
  # binary helper array
  shape = np.asarray(my_grid.grid.data.shape[0:3]) + np.array((2,2,2))
  helper = np.zeros(shape,dtype=int)
  # use voxel info for Marching cubes algorithm
  # set voxels on boundary wit value 0
  # set voxels inside structure with value 1
  # voxels outside structure have value -1
  for i in range(shape[0]-2):
    for j in range(shape[1]-2):
      for k in range(shape[2]-2):
        if my_grid.grid.data[i,j,k] != -1: # valid voxel
          helper[i+1,j+1,k+1] = 1
  
  for b in borders:
    # b[0] stores direction oft cartesian comm
    if b[0] > my_grid.rank:
      if b[1] is not None:
        sys.stdout.flush()
        helper[shape[0]-1,1:helper.shape[1]-1,1:helper.shape[1]-1] = b[1] 
    else:
      assert(b[0] < my_grid.rank)
      if b[1] is not None:
        helper[0,1:helper.shape[1]-1,1:helper.shape[1]-1] = b[1]
        sys.stdout.flush()
  
  hx_help = width[0] / float(shape[0])
  hy_help = width[1] / float(shape[1])
  hz_help = width[2] / float(shape[2])
  
  print("hx_help,hy_help,hz_help:",hx_help,hy_help,hz_help)
  print("hx,hy,hz:",my_grid.grid.hx,my_grid.grid.hy,my_grid.grid.hz)
  
  # vtk rectgrid is node based
  x = np.arange(my_grid.bounds[0],my_grid.bounds[3]+1*hx_help,hx_help)
  y = np.arange(my_grid.bounds[1],my_grid.bounds[4]+1*hy_help,hy_help)
  z = np.arange(my_grid.bounds[2],my_grid.bounds[5]+1*hz_help,hz_help)
  from pyevtk.hl import gridToVTK  
  gridToVTK("helper"+str(my_grid.rank),x,y,z,cellData={"helper":helper}) 
 
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, normals, values = measure.marching_cubes(helper,spacing=(np.float32(my_grid.grid.hx),np.float32(my_grid.grid.hy),np.float32(my_grid.grid.hz)),allow_degenerate=False)
  verts = np.asarray(verts) + (my_grid.grid.hx/2.0,my_grid.grid.hx/2.0,my_grid.grid.hx/2.0)
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(my_grid.bounds[0],my_grid.bounds[1],my_grid.bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
  matviz_vtk.show_write_vtk(transformFilter.GetOutput(), 10, "marching"+str(my_grid.rank)+".vtp")
  
  connectivity = basecell.getConnectivity(verts,faces)
  verts = basecell.taubin_smoothing(verts,connectivity)
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(my_grid.bounds[0],my_grid.bounds[1],my_grid.bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
  
  normals = vtk.vtkPolyDataNormals()
  normals.SetInputData(transformFilter.GetOutput())
  normals.SetConsistency(1)
  normals.SetAutoOrientNormals(1)
  normals.Update()
  
  pd = normals.GetOutput()
  
  matviz_vtk.show_write_vtk(pd, 10, "smoothed_marching"+str(my_grid.rank)+".vtp")

  return pd

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

# draw 4 faces/triangles of a tetrahedron spanned by vertices A,B,C,D
#     C
#   /   \
#  /     \ 
# A-------B
# triangle plane equation: A + u * (B-A) + v * (C - A), u,v \in [0,1] and u+v <= 1
# tetrahedron consists of 4 triangles: ABC, ABD, BCD, ADC
# @param value: value to draw with
def draw_tetrahedron(A,B,C,D,grid,value=0):
  triangles = []
  
  # triangle ABC
  triangles.append((np.asarray(A),np.asarray(B),np.asarray(C)))
  # triangle ABD
  triangles.append((np.asarray(A),np.asarray(B),np.asarray(D)))
  # triangle BCD
  triangles.append((np.asarray(B),np.asarray(C),np.asarray(D)))
  # triangle ADC
  triangles.append((np.asarray(A),np.asarray(D),np.asarray(C)))
  
  eps = 1e-6
  
  for t in triangles:
    draw_triangle(t[0], t[1], t[2], grid, value)

# @param v0, v1, v2: triangle vertices        
# @param value: value to draw with        
def draw_triangle(v0,v1,v2,grid,value):
  # directional vectors
  v1v0 = v1 - v0
  v2v0 = v2 - v0
  
  # number of samples
  samples10 = abs(max(v1v0))
  samples20 = abs(max(v2v0))
  
  eps = 1e-6
  
  for u in np.linspace(0,1,num=samples10+10):
    for v in np.linspace(0,1,num=samples20+10):
      if u + v > 1.0:
        continue
      
      p0 = v0 + u * v1v0 + v * v2v0 + eps
      p = [int(v) for v in p0]
      
      if not idx_out_of_bounds(p, grid.shape):
        grid[tuple(p)] = value
          
        
# @param elem: list of lists, each entry contains 4 vertices (cartesian) of a tet
# @param bounds: list of bounds of cartesian world
# @param h: lattice spacings in 3 directions 
# @param grid: where to draw        
# @param value: value to draw with        
def draw_non_design(tets,grid,bounds,h,value):
  from scipy.ndimage import binary_fill_holes
  assert(len(bounds) >=3 and len(h) == 3)
  for e in tets:
    assert(len(e) == 4)
    
    # A: (i,j,k)
    A = np.asarray(draw_profile_functions.cartesian_to_voxel_coords(e[0],bounds[0],bounds[1],bounds[2],h[0],h[1],h[2]))
    B = np.asarray(draw_profile_functions.cartesian_to_voxel_coords(e[1],bounds[0],bounds[1],bounds[2],h[0],h[1],h[2]))
    C = np.asarray(draw_profile_functions.cartesian_to_voxel_coords(e[2],bounds[0],bounds[1],bounds[2],h[0],h[1],h[2]))
    D = np.asarray(draw_profile_functions.cartesian_to_voxel_coords(e[3],bounds[0],bounds[1],bounds[2],h[0],h[1],h[2]))
    
    #     print("A:",A)
    #     print("B:",B)
    #     print("C:",C)
    #     print("D:",D)
    #       sys.exit()
    
    # find out which vertex is inside and which partial triangle to draw
    vertices = [A,B,C,D]
    vertex_outside = [idx_out_of_bounds(v, grid.shape) for v in vertices]
    
    # whole tet lies outside local grid
    if np.all(vertex_outside):
      continue  
      # at least one vertex is inside domain
    elif np.any(np.invert(vertex_outside)):
      # list of vertices inside domain
      valid_verts_idx = np.where(np.invert(vertex_outside))[0]
      for idx in valid_verts_idx:
        grid[tuple(vertices[idx])] = value
      
      draw_tetrahedron(A, B, C, D, grid, value)
      
    else: # draw complete tet
      assert(not (idx_out_of_bounds(A, grid.shape) or idx_out_of_bounds(B, grid.shape) or idx_out_of_bounds(C, grid.shape)))
      
      # works only if A is a tuple
      grid[tuple(A)] = value 
      grid[tuple(B)] = value
      grid[tuple(C)] = value
      grid[tuple(D)] = value
      
      draw_tetrahedron(A,B,C,D,grid,value)
      
      # get bounding box of this tet and fill holes inside 4 drawn triangles
      mindim = [min(A[i],B[i],C[i],D[i]) for i in range(3)]
      maxdim = [max(A[i],B[i],C[i],D[i]) for i in range(3)] 
      
      subgrid = grid[int(mindim[0]):int(maxdim[0])+1,int(mindim[1]):int(maxdim[1])+1,int(mindim[2]):int(maxdim[2])+1]
      helper = np.ones(subgrid.shape) * 999
      # convert to binary 3dimage, 0 is hole and 1 structure
      for i in range(subgrid.shape[0]):
        for j in range(subgrid.shape[1]):
          for k in range(subgrid.shape[2]):
            # solid
            if subgrid[i,j,k] == value:
              helper[i,j,k] = 1
            else:
              #void
              assert(subgrid[i,j,k] != value)
              helper[i,j,k] = 0
       
      helper = binary_fill_holes(helper).astype(int)
       
      for i in range(helper.shape[0]):
        for j in range(helper.shape[1]):
          for k in range(helper.shape[2]):
            if helper[i,j,k] == 1:
              subgrid[i,j,k] = value
              
      grid[int(mindim[0]):int(maxdim[0])+1,int(mindim[1]):int(maxdim[1])+1,int(mindim[2]):int(maxdim[2])+1] = subgrid           

# check if given tuple of indices lies within given array bounds
# @param point (i,j,k)
# @param bounds: list of 3 ints    
def idx_out_of_bounds(point,bounds):
  if 0 <= point[0] < bounds[0] and 0 <= point[1] < bounds[1] and 0 <= point[2] < bounds[2]:
    return False
  else:
    return True
  
#  Setup 2D distributed grid with shared ghost boundaries
#  using mpi4py and the MPI Cartesian Communicator. 

class DistributedGrid():
  # total_samples: list with number of total samples in 3 directions
  # bc_res: resolution of one base cell (usually 40)
  # bounds: global (xmin,ymin,zmin,xmax,ymax,zmax)
  def __init__(self, comm, total_samples,bc_res,bounds):
    assert(len(total_samples) == 3)
    assert(len(bounds) == 6)
    
    self.rank, self.size = comm.Get_rank(), comm.Get_size()
    self.cart = comm.Create_cart(dims=(1,self.size))
    
    self.coords = self.cart.Get_coords(self.rank)
    
    # distribute along x-axis    
    # number of chunks for this rank
    # start_x,end_x: first and last x-slice for this rank
    self.chunks, self.start_x, self.end_x = self.calc_num_cell_slices(total_samples[0])
    self.bounds = bounds[:]
    # only in x-direction different from global ones
    self.bounds[0] = self.start_x * (bounds[3]-bounds[0]) / total_samples[0] + bounds[0]
    self.bounds[3] = self.end_x * (bounds[3]-bounds[0]) / total_samples[0] + bounds[0]
    self.grid = RectGrid(int(self.chunks*bc_res),int(total_samples[1]*bc_res),int(total_samples[2]*bc_res), self.bounds)
       
  # calculate number of slices, start and end idx along x-axis for this rank
  # @param total_slices: number of total slices
  def calc_num_cell_slices(self,total_slices):
    # 1D Slab Decomposition of global voxel grid: http://www.2decomp.org/decomp.html 
    # distribute number of cells along x-axis evenly over number of processes and take care
    # of remainder -> difference between work chunks can be 1
    chunks = [int(total_slices/self.size) + (1 if p < total_slices%self.size else 0) for p in range(0,self.size)]
    start_slice = int(np.sum(chunks[0:self.rank]))
    end_slice = int(start_slice + chunks[self.rank])
    
    return chunks[self.rank], start_slice, end_slice
  
  def to_info(self):
    print("---- mpi distr grid ----")  
    print("rank:",self.rank)
    print("bounds:",self.bounds)
    self.grid.to_info()
    sys.stdout.flush()
  
  # communicate ghost layers in both directions  
  def communicate_edges(self):
    recv = []
    for direction in [-1,1]:
      # coordinate dimension of shift is always 1, because we only have one dimension
      # coordinate dimension 1-based
      source, dest = self.cart.Shift(1,direction)
      print("rank ", self.rank, " sending slice ",self.grid.nx-1," to rank ",dest)
#       print("rank ", self.rank," direction:",direction," left:",source," right:",dest)
      sys.stdout.flush()

      sendbuf = None
      if direction == -1:
        sendbuf = np.copy(self.grid.data[0,:,:])
      else:            
        sendbuf = np.copy(self.grid.data[self.grid.nx-1,:,:])
      
      assert(sendbuf is not None)
      recvbuf = np.full_like(sendbuf,999)
      
      helper = np.zeros_like(sendbuf, dtype=int)
      for i in range(helper.shape[0]):
        for j in range(helper.shape[1]):
            if sendbuf[i,j] != -1:
              helper[i,j] = 1
              
      sendbuf = helper        
      
      self.cart.Sendrecv(sendbuf,dest=dest,source=source,recvbuf=recvbuf)
      
      if source == MPI.PROC_NULL:
        # got no data, cause neighbor does not exist
        recv.append((source,None))
      else:
        recv.append((source,recvbuf))
    
    return recv
      
class RectGrid():
  def __init__(self,nx,ny,nz,bounds):
    self.data = np.full((nx,ny,nz),999,dtype=int)

    self.nx = nx
    self.ny = ny
    self.nz = nz

    width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]    
    self.hx = width[0] / nx
    self.hy = width[1] / ny
    self.hz = width[2] / nz
  
  def to_info(self):  
    print("res:",self.nx,self.ny,self.nz)
    print("spacing:",self.hx,self.hy,self.hz)
    sys.stdout.flush()
    
    