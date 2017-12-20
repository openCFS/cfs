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
  
  basecells = []
  nProblem = nx*ny*nz

  # distribute loop evenly over number of processes and take care
  # of remainder -> difference between work chunks can be 1
  # e.g. 10 loop runs, 4 processes, number of runs per process: [3,3,2,2] 
  # number of loop runs per process:
  nRuns = [int(nProblem/commSize) + (1 if p < nProblem%commSize else 0) for p in range(0,commSize)]
  start = int(np.sum(nRuns[0:rank]))
  end = int(start + nRuns[rank])
  
  print("rank:",rank," start:",start," end:",end)
  sys.stdout.flush()
  
  
  ################# voxel world #########################################
  nondes_coords = nondes[0]
  nondes_min = nondes[1]
  nondes_max = nondes[2]
  assert(len(nondes_min) == 3)
  assert(len(nondes_max) == 3) 
  # np.minimum gives elementwise min value
  bounds[0:3] = np.minimum(np.asarray(min_bb),np.asarray(nondes_min))
  bounds[3:6] = np.maximum(np.asarray(max_bb),np.asarray(nondes_max))
  print("min_bb:",min_bb)
  print("max_bb:",max_bb)
  print("min_dim:",bounds[0:3])
  print("max_dim:",bounds[3:6])
  resolution = (int(samples[0]*args.bc_res),int(samples[1]*args.bc_res),int(samples[2]*args.bc_res))
  print("resolution:",resolution)
  voxel_grid  = np.full(resolution,-1,dtype=int)
  width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
  
  hx_old = width[0] / float(resolution[0])
  hy_old = width[1] / float(resolution[1])
  hz_old = width[2] / float(resolution[2])
  
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
    
    print("rank:",rank," i,j,k:",i,j,k, " before basecell;",x1,x2,y1,y2,z1,z2)
    sys.stdout.flush()
    
    bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="surface_mesh",bc_flags=flags)
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
    voxel_grid[i*args.bc_res:(i+1)*args.bc_res,j*args.bc_res:(j+1)*args.bc_res,k*args.bc_res:(k+1)*args.bc_res] = cell_obj.voxels
    
    basecells.append((cell_obj,(i,j,k)))
    print("appended ",i,j,k,left_front_corner,x1,x2,y1,y2,z1,z2)
  
  eps = 1e-6
  x = np.arange(bounds[0],bounds[3]+hx_old-eps,hx_old)
  y = np.arange(bounds[1],bounds[4]+hy_old-eps,hy_old)
  z = np.arange(bounds[2],bounds[5]+hz_old-eps,hz_old)
  
  for p in nondes_coords:
    i,j,k = draw_profile_functions.cartesian_to_voxel_coords(p,bounds[0],bounds[1],bounds[2],hx_old,hy_old,hz_old)
    voxel_grid[i,j,k] = -2
  
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
          
#   helper[0,:,:] = helper[1,:,:]
#   helper[:,0,:] = helper[:,1,:]
#   helper[:,:,0] = helper[:,:,1]
#   helper[shape[0]-1,:,:] = helper[shape[0]-2,:,:]
#   helper[:,shape[1]-1,:] = helper[:,shape[1]-2,:]
#   helper[:,:,shape[2]-1] = helper[:,:,shape[2]-2]
  
  hx = width[0] / float(shape[0])
  hy = width[1] / float(shape[1])
  hz = width[2] / float(shape[2])
  
  print("hx,hy,hz:",hx,hy,hz)
  print("old hx,hy,hz:",hx_old,hy_old,hz_old)
  
  x = np.arange(bounds[0],bounds[3]+1*hx,hx)
  y = np.arange(bounds[1],bounds[4]+1*hy,hy)
  z = np.arange(bounds[2],bounds[5]+1*hz,hz)
  from pyevtk.hl import gridToVTK  
  gridToVTK("helper",x,y,z,cellData={"helper":helper}) 
 
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, normals, values = measure.marching_cubes(helper,spacing=(np.float32(hx_old),np.float32(hy_old),np.float32(hz_old)),allow_degenerate=False)
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
      
    import time   
    start = time.time()
    # for appending base cells
    vertices = []   
    faces = []
    # appending cells
    appends = vtk.vtkAppendPolyData()
    # each bc (entry of basecells list) stores the basecell object and lists with boundary circle meshes      
    for i,cell in enumerate(basecells):
      if (i%100==0):
        end_vtk_iteration=time.time()
        print ("Time for vtk "+str(i)+" iterations" ,end_vtk_iteration-start ," s "   )
      
      pd = matviz_vtk.fill_vtk_polydata(cell.points, cell.cells)
      appends.AddInputData(pd)
    
    appends.Update()
    
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
