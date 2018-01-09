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
  rank = comm.Get_rank()
  commSize = comm.Get_size()
  
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
    
  # where we want nodes
  nx = int(h_des[0] / dx_des)
  ny = int(h_des[1] / dy_des)
  nz = int(h_des[2] / dz_des)

#   print("samples:",samples)  
#   print("h_des:",h_des)
#   print("design_bounds:",design_bounds)  
#   print("dx_des,dy_des,dz_des:",dx_des,dy_des,dz_des)
#   print("nx,ny,nz:",nx,ny,nz) 
  
  data_grid, data_grid_near, sample_coords= matviz_vtk.get_interpolation_row_major(coords, design_bounds, grad, s1, s2, s3, nx, ny, nz, dx_des, dy_des, dz_des)
  
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
  print("min_bb design:",min_bb)
  print("max_bb design:",max_bb)
  print("min:",bounds[0:3])
  print("max:",bounds[3:6])
  resolution = (int(samples[0]*args.bc_res),int(samples[1]*args.bc_res),int(samples[2]*args.bc_res))
  print("voxel grid res:",resolution)
  voxel_grid  = np.full(resolution,999,dtype=int)
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  
  hx = width[0] / float(resolution[0])
  hy = width[1] / float(resolution[1])
  hz = width[2] / float(resolution[2])
  
  nCells = nx*ny*nz
  # distribute loop evenly over number of processes and take care
  # of remainder -> difference between work chunks can be 1
  # e.g. 10 loop runs, 4 processes, number of runs per process: [3,3,2,2]
  # number of loop runs per process:
  nRuns = [int(nCells/commSize) + (1 if p < nCells%commSize else 0) for p in range(0,commSize)]
  start = int(np.sum(nRuns[0:rank]))
  end = int(start + nRuns[rank])
  
#   # 1D Slab Decomposition of global voxel grid: http://www.2decomp.org/decomp.html 
#   # distribute slices along x-axis evenly over number of processes and take care
#   # of remainder -> difference between work chunks can be 1
#   chunks = [int(resolution[0]/commSize) + (1 if p < resolution[0]%commSize else 0) for p in range(0,commSize)]
#   start = int(np.sum(chunks[0:rank]))
#   end = int(start + chunks[rank])
  
#   print("chunks:",chunks)
#   print("rank:",rank," start:",start," end:",end)
#   sys.stdout.flush()

#   if rank == 0: 
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
    
    print("rank:",rank," i,j,k:",i,j,k, " before basecell;",x1,x2,y1,y2,z1,z2)
    sys.stdout.flush()
    
    bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
    bc_input.eta = 0.7
    bc_input.stiffness_as_diameter = True
    cell_obj = basecell.Basecell(bc_input,id)
#     cell_obj = basecell.Basecell(bc_input,id,nondes,(left_front_corner,right_back_corner))
#     cell_obj.scale(dx, dy, dz)
#     cell_obj.translate(left_front_corner[0],left_front_corner[1],left_front_corner[2])
#     cell_obj.update()
#     if x1 <= 0.076 and x2 <= 0.076 and y1 <= 0.076 and y2 <= 0.076 and z1 <= 0.076 and z2 <= 0.076:
#       continue
#     cell_center = np.asarray(left_front_corner) + np.asarray([dx/2,dy/2,dz/2])
#     cell_obj.center = cell_center
    voxel_grid[i*args.bc_res:(i+1)*args.bc_res,j*args.bc_res:(j+1)*args.bc_res,k*args.bc_res:(k+1)*args.bc_res] = cell_obj.voxels
    
  eps = 1e-6
  x = np.arange(bounds[0],bounds[3]+hx-eps,hx)
  y = np.arange(bounds[1],bounds[4]+hy-eps,hy)
  z = np.arange(bounds[2],bounds[5]+hz-eps,hz)

  void_elems = []
  # read void non-design manually  
  import csv
  with open('holes.txt') as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')  
    for row in readCSV:
      # 4 elems with 3 coord components
      assert(len(row) == 12)
      void_elems.append([(float(row[0]),float(row[1]),float(row[2])), (float(row[3]),float(row[4]),float(row[5])), (float(row[6]),float(row[7]),float(row[8])), (float(row[9]),float(row[10]),float(row[11]))])
    
    assert(len(void_elems) > 0)
  
  draw_non_design(nondes_coords, voxel_grid, bounds, (hx,hy,hz), 0)
  draw_non_design(void_elems, voxel_grid, bounds, (hx,hy,hz), -1)    
  from pyevtk.hl import gridToVTK  
  gridToVTK("voxels",x,y,z,cellData={"voxels":voxel_grid})
  
  # binary helper array
  shape = np.asarray(voxel_grid.shape[0:3]) + np.array((2,2,2))
  helper = np.zeros(shape,dtype=int)
  # use voxel info for Marching cubes algorithm
  # set voxels on boundary wit value 0
  # set voxels inside structure with value 1
  # voxels outside structure have value -1
  for i in range(shape[0]-2):
    for j in range(shape[1]-2):
      for k in range(shape[2]-2):
        if voxel_grid[i,j,k] != -1: # valid voxel
          helper[i+1,j+1,k+1] = 1
          
  hx_help = width[0] / float(shape[0])
  hy_help = width[1] / float(shape[1])
  hz_help = width[2] / float(shape[2])
  
  print("hx_help,hy_help,hz_help:",hx_help,hy_help,hz_help)
  print("hx,hy,hz:",hx,hy,hz)
  
  # vtk rectgrid is node based
  x = np.arange(bounds[0],bounds[3]+1*hx_help,hx_help)
  y = np.arange(bounds[1],bounds[4]+1*hy_help,hy_help)
  z = np.arange(bounds[2],bounds[5]+1*hz_help,hz_help)
  from pyevtk.hl import gridToVTK  
  gridToVTK("helper",x,y,z,cellData={"helper":helper}) 
 
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, normals, values = measure.marching_cubes(helper,spacing=(np.float32(hx),np.float32(hy),np.float32(hz)),allow_degenerate=False)
  verts = np.asarray(verts) + (hx/2.0,hy/2.0,hz/2.0)
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(bounds[0],bounds[1],bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
  matviz_vtk.show_write_vtk(transformFilter.GetOutput(), 10, "marching.vtp")
  
  connectivity = basecell.getConnectivity(verts,faces)
  verts = basecell.taubin_smoothing(verts,connectivity,30)
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  translation = vtk.vtkTransform()
  translation.Translate(bounds[0],bounds[1],bounds[2])
  transformFilter = vtk.vtkTransformPolyDataFilter()
  transformFilter.SetInputData(pd)
  transformFilter.SetTransform(translation)
  transformFilter.Update()
 
  matviz_vtk.show_write_vtk(transformFilter.GetOutput(), 10, "smooothed_marching.vtp")
   
  normals = vtk.vtkPolyDataNormals()
  normals.SetInputData(transformFilter.GetOutput())
  normals.SetConsistency(1)
  normals.SetAutoOrientNormals(1)
  normals.Update()
  
  return normals.GetOutput()

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
    v0 = t[0]
    v1 = t[1]
    v2 = t[2]
    
    # directional vectors
    v1v0 = v1 - v0
    v2v0 = v2 - v0
    
    # number of samples
    samples10 = abs(max(v1v0))
    samples20 = abs(max(v2v0))
    
    for u in np.linspace(0,1,num=samples10+10):
      for v in np.linspace(0,1,num=samples20+10):
        if u + v > 1.0:
          continue
        
        p0 = v0 + u * v1v0 + v * v2v0 + eps
        p = [int(v) for v in p0]
        
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
