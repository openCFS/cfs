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
  # nondes: store info on solid and void nondesign regions, has 2 entries: 0 -> solid non-design, 1 -> void non_design
  # nondes[0]: solid nondesign -> (centers, min_bb, max_bb, elem_dim)
  # nondes[0][centers]: list of elements (corner vertices) that define non-design regions
  
  # MPI_Init() or MPI_Init_thread() is actually called when you import the MPI
  # use the standard communicator
  comm = MPI.COMM_WORLD
  
  # point coordinates from h5 file
  centers, _, _ = coords[0:3]
  
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
  
  my_mpi_grid = MPI_Grid(comm)
  
  ################# voxel world #########################################
  nondes_min = 99999
  nondes_max = -99999
  nondes_coords = None
  holes = None
  design_elems = None
  if nondes and nondes[0] and nondes[1]: # assume holes are within bounding box of solid non-design and design domain
    nondes_coords = nondes[0][0]
    nondes_min = nondes[0][1]
    nondes_max = nondes[0][2]
    assert(len(nondes_min) == 3)
    assert(len(nondes_max) == 3)
    holes = nondes[1][0]
    holes_min = nondes[1][1]
    holes_max = nondes[1][2]
    assert(len(holes_min) == 3)
    assert(len(holes_max) == 3)
    design_elems = nondes[2][0]
    design_elems_min = nondes[2][1]
    design_elems_max = nondes[2][2]

  # broadcast all nondes to all ranks
  (nondes_coords,nondes_min,nondes_max) = my_mpi_grid.comm.bcast((nondes_coords,nondes_min,nondes_max),root=0)
#   print("rank ", my_mpi_grid.rank," nondes:",len(nondes_coords))
  # broadcast all holes and design elems to all ranks
  holes = my_mpi_grid.comm.bcast(holes,root=0)
  design_elems = my_mpi_grid.comm.bcast(design_elems,root=0)        
    
  # np.minimum gives elementwise min value
  bounds = [None] * 6
  bounds[0:3] = np.minimum(np.asarray(min_bb),np.asarray(nondes_min))
  bounds[3:6] = np.maximum(np.asarray(max_bb),np.asarray(nondes_max))
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  print("\nmin_bb design:",min_bb)
  print("max_bb design:",max_bb)
  print("min:",bounds[0:3])
  print("max:",bounds[3:6])
  
  my_mpi_grid.init_data_grid(samples, args.bc_res, bounds)
  my_mpi_grid.to_info()
  
  local_id = 0
  for k in range(samples[2]):
    for j in range(samples[1]):
      for i in range(my_mpi_grid.start_x,my_mpi_grid.end_x):
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

        li,lj,lk = get_3d_grid_coords(local_id, my_mpi_grid.chunks, samples[1], samples[2])        
        # bounds (voxel coords) of local base cell
        # xmin,ymin,zmin
        h = (my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz)
        left,lower,back = draw_profile_functions.voxel_to_cartesian_coords((li*args.bc_res,lj*args.bc_res,lk*args.bc_res), design_bounds[0:3], h)
        # xmax,ymax,zmax
        right,upper,front = draw_profile_functions.voxel_to_cartesian_coords(((li+1)*args.bc_res,(lj+1)*args.bc_res,(lk+1)*args.bc_res), design_bounds[0:3], h)
        
        left_lower_back = (left,lower,back)
        left_upper_back = (left,upper,back)
        left_lower_front = (left,lower,front)
        left_upper_front = (left,upper,front)
        
        right_lower_back = (right,lower,back)
        right_upper_back = (right,upper,back)
        right_lower_front = (right,lower,front)
        right_upper_front = (right,upper,front)
        
        # if base cell corners outside design domain, don't compute
        if out_of_bounds(left_lower_back, design_bounds) and out_of_bounds(left_upper_back, design_bounds)\
          and out_of_bounds(left_lower_front, design_bounds) and out_of_bounds(left_upper_front, design_bounds)\
          and out_of_bounds(right_lower_back, design_bounds) and out_of_bounds(right_upper_back, design_bounds)\
          and out_of_bounds(right_lower_front, design_bounds) and out_of_bounds(right_upper_front, design_bounds):
            continue
        
        flags = None
        bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
        bc_input.eta = 0.7
        bc_input.stiffness_as_diameter = True
        cell_obj = basecell.Basecell(bc_input,id)
        # local i,j,k
        print("rank:",my_mpi_grid.rank," global i,j,k:",i,j,k, " local:",li,lj,lk," x1,x2,y1,y2,z1,z2:",x1,x2,y1,y2,z1,z2)
        my_mpi_grid.grid.data[li*args.bc_res:(li+1)*args.bc_res,lj*args.bc_res:(lj+1)*args.bc_res,lk*args.bc_res:(lk+1)*args.bc_res] = cell_obj.voxels
         
        local_id += 1
    
  eps = 1e-6

  if nondes:
    hx = (my_mpi_grid.bounds[3]-my_mpi_grid.bounds[0])/my_mpi_grid.grid.nx
    hy = (my_mpi_grid.bounds[4]-my_mpi_grid.bounds[1])/my_mpi_grid.grid.ny
    hz = (my_mpi_grid.bounds[5]-my_mpi_grid.bounds[2])/my_mpi_grid.grid.nz
    #x = np.arange(my_mpi_grid.bounds[0],my_mpi_grid.bounds[3]+hx-eps,hx)
    #y = np.arange(my_mpi_grid.bounds[1],my_mpi_grid.bounds[4]+hy-eps,hy)
    #z = np.arange(my_mpi_grid.bounds[2],my_mpi_grid.bounds[5]+hz-eps,hz)
    
    tmp = np.zeros_like(my_mpi_grid.grid.data,dtype=bool)
    draw_non_design(design_elems, tmp, my_mpi_grid.bounds,(my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=True)
    my_mpi_grid.grid.data *= tmp
    
    design_elems = None
    draw_non_design(nondes_coords, my_mpi_grid.grid.data, my_mpi_grid.bounds,(my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=True)
#     draw_non_design(holes, my_mpi_grid.grid.data, my_mpi_grid.bounds, (my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=False)
  #   x = np.arange(0,my_mpi_grid.grid.data.shape[0]+1,1)
  #   y = np.arange(0,my_mpi_grid.grid.data.shape[1]+1,1)
  #   z = np.arange(0,my_mpi_grid.grid.data.shape[2]+1,1)
  #   from pyevtk.hl import gridToVTK
  #   gridToVTK("voxels"+str(my_mpi_grid.rank),x,y,z,cellData={"voxels":my_mpi_grid.grid.data.astype(int)})
  
    nondes_coords = None
    holes = None
    
  borders = my_mpi_grid.communicate_edges()
      
  # binary helper array
  shape = np.asarray(my_mpi_grid.grid.data.shape[0:3]) + np.array((2,2,2))
  helper = np.zeros(shape,dtype=bool)
  helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = my_mpi_grid.grid.data
  
  # hope for python's garbage collector to delete voxel array
  my_mpi_grid.grid.data = None
  
  for b in borders:
    # b[0] stores direction oft cartesian comm
    if b[0] > my_mpi_grid.rank:
      assert(b[1] is not None)
      helper[shape[0]-1,1:shape[1]-1,1:shape[2]-1] = b[1] 
    else:
      assert(b[0] < my_mpi_grid.rank)
      assert(b[1] is not None)
      helper[0,1:shape[1]-1,1:shape[2]-1] = b[1]
  
  print("hx,hy,hz:",my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz)
  hx = my_mpi_grid.grid.hx
  hy = my_mpi_grid.grid.hy
  hz = my_mpi_grid.grid.hz
  
  #print(helper.shape,np.sum(helper))
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, _, _ = measure.marching_cubes(helper,spacing=(np.float32(my_mpi_grid.grid.hx),np.float32(my_mpi_grid.grid.hy),np.float32(my_mpi_grid.grid.hz)),allow_degenerate=False,step_size=1)
  # translate from (0,0,0) to correct position
  # and marching cubes shift everything by 0.5 * hx/hy/hz
  shift = np.asarray(my_mpi_grid.bounds[0:3]) - np.asarray((my_mpi_grid.grid.hx/2.0,my_mpi_grid.grid.hy/2.0,my_mpi_grid.grid.hz/2.0))
  print("shift:",shift)
  verts = [p+shift for p in verts]
  
  # hope for python's garbage collector to delete voxel array
  helper = None
  
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  matviz_vtk.show_write_vtk(pd, 10, "marching"+str(my_mpi_grid.rank)+".vtp")
  
  pd = None
  
  my_mpi_grid.set_vertices_and_faces(list(verts),list(faces))
  
  my_mpi_grid.gather_data(append=True)
  
  print("after gather data rank:",my_mpi_grid.rank)
  sys.stdout.flush()
  data = None
  if my_mpi_grid.rank == 0:
    import pymesh
    print("before reducing: len(verts):",len(my_mpi_grid.vertices)," len(faces):",len(my_mpi_grid.faces))
    sys.stdout.flush()
    
    verts = np.asarray(my_mpi_grid.vertices)
    faces = np.asarray(my_mpi_grid.faces)
    
    verts, faces, info = pymesh.remove_duplicated_vertices_raw(verts,faces,1e-4) 
    verts, faces, info = pymesh.remove_duplicated_faces_raw(verts,faces)
    print("after reducing: len(verts):",len(verts)," len(faces):",len(faces))
    sys.stdout.flush()
    
    pd = matviz_vtk.fill_vtk_polydata(verts, faces)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(pd)
    normals.SetConsistency(1)
    normals.SetAutoOrientNormals(1)
    normals.Update()
    sys.stdout.flush()
    matviz_vtk.show_write_vtk(pd, 10, "marching_all.vtp")
    
    data = (verts,faces)
    
  # broadcast all verts to all ranks
  data = my_mpi_grid.comm.bcast(data,root=0)
  my_mpi_grid.set_vertices_and_faces(data[0],data[1])
  
  for v in my_mpi_grid.vertices:
    assert(v is not None)
  assert(np.amax(my_mpi_grid.faces) < len(my_mpi_grid.vertices))
  
  # update after changing vertices and faces
  my_mpi_grid.update_connectivity_verts_faces()
  
  xmin = min(my_mpi_grid.vertices, key=lambda t: t[0])[0]
  xmax = max(my_mpi_grid.vertices, key=lambda t: t[0])[0]
  ymin = min(my_mpi_grid.vertices, key=lambda t: t[1])[1]
  ymax = max(my_mpi_grid.vertices, key=lambda t: t[1])[1]
#   zmin = min(my_mpi_grid.vertices, key=lambda t: t[2])[2]
  zmax = max(my_mpi_grid.vertices, key=lambda t: t[2])[2]
  
  # do parallel smoothing here
  mpi_taubin_smoothing(my_mpi_grid,(xmin,ymin,0,xmax,ymax,zmax),niter=args.bc_smooth)
  
  # send smoothed data to root
  my_mpi_grid.gather_data(append=False, root=0)
  
  if my_mpi_grid.rank != 0:
    sys.exit()
  
  pd = matviz_vtk.fill_vtk_polydata(my_mpi_grid.vertices,my_mpi_grid.faces)
  matviz_vtk.show_write_vtk(pd, 10, "smoothed_marching"+str(my_mpi_grid.rank)+".vtp")
    
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
  
  nondes_grid = np.full(resolution,0,dtype=bool)
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
        
# @param tets: list of lists, each entry contains 4 vertices (cartesian) of a tet
# @param grid: where to draw
# @param bounds: list of bounds of cartesian world
# @param h: lattice spacings in 3 directions 
# @param solid or void non-design?        
def draw_non_design(tets,grid,bounds,h,solid=True):
  assert(len(bounds) >=3 and len(h) == 3)
  ug = vtk.vtkUnstructuredGrid()
  vtkpoints = vtk.vtkPoints()
  for e in tets:
    #print("e:",e)
    tetra = vtk.vtkTetra()
    for i,p in enumerate(e):
      id = vtkpoints.InsertNextPoint(p)
      tetra.GetPointIds().SetId(i,id)
    
    ug.InsertNextCell(tetra.GetCellType(),tetra.GetPointIds())
  
  ug.SetPoints(vtkpoints)
  
  vtkpoints = None
  
  resample = vtk.vtkResampleToImage() 
  resample.SetInputDataObject(ug)
  resample.UseInputBoundsOff()
  resample.SetSamplingBounds(bounds[0],bounds[3],bounds[1],bounds[4],bounds[2],bounds[5])
  resample.SetSamplingDimensions(grid.shape)
  resample.Update()
  
  ug = None
  
  itp = vtk.vtkImageDataToPointSet()
  itp.SetInputDataObject(resample.GetOutput())
  itp.Update()
  
  resample = None
  
  from vtk.util import numpy_support
  fact = numpy_support.vtk_to_numpy(itp.GetOutput().GetPointData().GetArray('vtkValidPointMask')).astype(bool)
  points = numpy_support.vtk_to_numpy(itp.GetOutput().GetPoints().GetData())
  
  prod = [p for i,p in enumerate(points) if fact[i]]
  assert(prod is not None)
  for p in prod:
    # transfer to voxel world
    i,j,k = draw_profile_functions.cartesian_to_voxel_coords(p,bounds[0],bounds[1],bounds[2],h[0],h[1],h[2])
    assert(not idx_out_of_bounds((i,j,k),grid.shape))
    grid[i,j,k] = solid

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
class MPI_Grid():
  # init cartesian mpi world
  def __init__(self,comm):
    self.comm = comm    
    self.rank, self.size = comm.Get_rank(), comm.Get_size()
    self.cart = comm.Create_cart(dims=(1,self.size))
    self.coords = self.cart.Get_coords(self.rank)
  
  # total_samples: list with number of total samples in 3 directions
  # bc_res: resolution of one base cell (usually 40)
  # bounds: global (xmin,ymin,zmin,xmax,ymax,zmax)  
  def init_data_grid(self,total_samples,bc_res,bounds):
    assert(len(total_samples) == 3)
    assert(len(bounds) == 6)  
    # distribute along x-axis    
    # number of chunks for this rank
    # start_x,end_x: first and last x-slice for this rank
    self.chunks, self.start_x, self.end_x = self.calc_num_cell_slices(total_samples[0])
    self.bounds = bounds[:]
    # only in x-direction different from global ones
    self.bounds[0] = self.start_x * (bounds[3]-bounds[0]) / total_samples[0] + bounds[0]
    self.bounds[3] = self.end_x * (bounds[3]-bounds[0]) / total_samples[0] + bounds[0]
    self.grid = RectGrid(int(self.chunks*bc_res),int(total_samples[1]*bc_res),int(total_samples[2]*bc_res), self.bounds)
    
  def set_vertices_and_faces(self,verts,faces):
    # contains ALL vertices
    self.vertices = verts
    self.faces = faces
    self.update_connectivity_verts_faces()
    self.chunks_vertices, self.start_verts_idx, self.end_verts_idx = self.calc_num_vertices()
    print("rank ",self.rank," start:",self.start_verts_idx," end:",self.end_verts_idx," chunks:",self.chunks_vertices," total:",len(self.vertices))
  
  def update_connectivity_verts_faces(self):
    self.connectivity = basecell.getConnectivity(self.vertices, self.faces)
    
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
  
  def calc_num_vertices(self):
    assert(self.vertices is not None and len(self.vertices) > 0)
    num_verts = len(self.vertices)
    chunks = [int(num_verts/self.size) + (1 if p < num_verts%self.size else 0) for p in range(0,self.size)]
    # start idx of 'vertices' list
    start_idx = int(np.sum(chunks[0:self.rank]))
    end_idx = int(start_idx + chunks[self.rank])
    
    return chunks[self.rank], start_idx, end_idx
    
  def to_info(self):
    print("---- mpi distr grid ----")  
    print("rank:",self.rank)
    print("bounds:",self.bounds)
    print("start:",self.start_x," end:",self.end_x)
    self.grid.to_info()
    sys.stdout.flush()
  
  # communicate ghost layers in both directions  
  def communicate_edges(self):
    recv = []
    for direction in [-1,1]:
      # coordinate dimension of shift is always 1, because we only have one dimension
      # coordinate dimension 1-based
      source, dest = self.cart.Shift(1,direction)

      sendbuf = None
      if direction == -1:
        sendbuf = np.copy(self.grid.data[0,:,:])
      else:            
        sendbuf = np.copy(self.grid.data[self.grid.nx-1,:,:])
      
      assert(sendbuf is not None)
      recvbuf = np.full_like(sendbuf,999)
      
      self.cart.Sendrecv(sendbuf,dest=dest,source=source,recvbuf=recvbuf)
      
      if source != MPI.PROC_NULL:
        # got no data, cause neighbor does not exist
        recv.append((source,recvbuf))

    return recv
  
  def communicate_vertices(self):
    assert(self.vertices is not None)
    recv = []
    for direction in [-1,1]:
      # coordinate dimension of shift is always 1, because we only have one dimension
      # coordinate dimension 1-based
      source, dest = self.cart.Shift(1,direction)
      
      sendbuf = (self.start_verts_idx, self.end_verts_idx,np.copy(self.vertices[self.start_verts_idx:self.end_verts_idx]))
      
      tmp = self.cart.sendrecv(sendbuf,dest=dest,source=source)
#       print("rank:",self.rank," got data from ",source)
      
      if source != MPI.PROC_NULL:
        recv.append(tmp)
    
#     print("rank:",self.rank," got ",len(recv)," data chunks")  
    return recv
      
  # collect all points and cells on a single rank (0)
  # if append==True, append data coming from other ranks
  # if append == False, update data coming from other ranks 
  def gather_data(self,append,root=0):
    if append:
      # find out number of vertices from each rank
      vertscount = self.comm.gather(3*len(self.vertices),root=root)
      facescount = self.comm.gather(3*len(self.faces),root=root)
     
      vertsbuf = None
      facesbuf = None
      if self.rank == root:
        vertsbuf = np.empty(int(np.sum(vertscount)))
        facesbuf = np.empty(int(np.sum(facescount)),dtype=int)
      self.comm.Gatherv(sendbuf=np.array(self.vertices),recvbuf=(vertsbuf,vertscount),root=root)
      self.comm.Gatherv(sendbuf=np.array(self.faces),recvbuf=(facesbuf,facescount),root=root)
      if self.rank == root:
        self.vertices = np.reshape(vertsbuf,(int(np.sum(vertscount)/3),3))
        self.faces = np.reshape(facesbuf,(int(np.sum(facescount)/3),3))
        # add offset
        offset = 0
        # cumulative sum
        cumsum_verts = [0]
        cumsum_verts.extend(np.cumsum(np.array(vertscount)/3.0))
        cumsum_verts = [int(id) for id in cumsum_verts]
        
        cumsum_faces = [0]
        cumsum_faces.extend(np.cumsum(np.array(facescount)/3.0))
        cumsum_faces = [int(id) for id in cumsum_faces]
          
        for r in range(len(cumsum_faces)-1):
          self.faces[cumsum_faces[r]:cumsum_faces[r+1]] += cumsum_verts[r]
        sys.stdout.flush()
    else: # update
      sendbuf = self.vertices[self.start_verts_idx:self.end_verts_idx]
      numverts = self.comm.gather(3*len(sendbuf),root=root)
      recvbuf = None
      if self.rank == root:
        recvbuf = np.empty(int(np.sum(numverts)))
      self.comm.Gatherv(sendbuf=np.array(sendbuf),recvbuf=(recvbuf,numverts),root=root)
      if self.rank == root:
        self.vertices = np.reshape(recvbuf,(int(np.sum(numverts)/3),3))  
        
  def update_vertices(self,recv):
    for r in recv:
      if r is not None:
        start = r[0]
        end = r[1]
        assert(r[2] is not None and r[2] is not None)
        self.vertices[start:end] = r[2]
        
class RectGrid():
  def __init__(self,nx,ny,nz,bounds):
    self.data = np.zeros((nx,ny,nz),dtype=bool)

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

def mpi_taubin_smoothing(mpi_grid,bounds=None,niter=None):
  assert(mpi_grid.vertices is not None)
  assert(mpi_grid.faces is not None)
  assert(mpi_grid.connectivity is not None)
  # smoothing parameter: p_i = p_i + lambda*L(p_{i,j})
  lamb = 0.4
  iter = 0
  res = 999
  old_points = None
  start = mpi_grid.start_verts_idx
  end = mpi_grid.end_verts_idx
  #while res > 1e-2:
  while iter < niter:
    sys.stdout.flush()
    old_points = mpi_grid.vertices.copy()
    assert(max(max(mpi_grid.connectivity)) < len(mpi_grid.vertices))
    # shrink
    mpi_grid.vertices = basecell.laplacian_smoothing(mpi_grid.vertices,mpi_grid.connectivity,lamb,start=start,end=end,rank=mpi_grid.rank,bounds=bounds)
    # communicate smoothed vertices to other ranks
    mpi_grid.update_vertices(mpi_grid.communicate_vertices())
    
    # expand
    mpi_grid.vertices = basecell.laplacian_smoothing(mpi_grid.vertices,mpi_grid.connectivity,-lamb-0.04,start=start,end=end,rank=mpi_grid.rank,bounds=bounds)
    assert(mpi_grid.vertices is not None)
    # communicate smoothed vertices to other ranks
    mpi_grid.update_vertices(mpi_grid.communicate_vertices())
    
    res = basecell.residual(old_points[start:end], mpi_grid.vertices[start:end])
    print("rank:", mpi_grid.rank, " iter:", iter, " residual:", res)
    iter += 1
    
  print("Taubin smoothing with ", iter, " iterations and res=",res)  

# @param point: list with 3 components
# @param bounds: list with 6 entries (xmin,ymin,zmin,xmax,ymax,zmax)
def out_of_bounds(point,bounds):
  eps = 1e-6
  if bounds[0]-eps <= point[0] <= bounds[3]+eps and bounds[1]-eps <= point[1] <= bounds[4]+eps and bounds[2]-eps <= point[2] <= bounds[5]+eps:
    return False
  else:
    print("point ",point, " is out of bounds ",design_bounds)
    return True