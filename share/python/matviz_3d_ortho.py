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
  
  my_mpi_grid = MPI_Grid(comm,samples,args.bc_res,bounds)
  #my_mpi_grid.to_info()
  
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
        
        flags = None
        bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
        bc_input.eta = 0.7
        bc_input.stiffness_as_diameter = True
        cell_obj = basecell.Basecell(bc_input,id)
        # local i,j,k
        li,lj,lk = get_3d_grid_coords(local_id, my_mpi_grid.chunks, samples[1], samples[2])
        print("rank:",my_mpi_grid.rank," global i,j,k:",i,j,k, " local:",li,lj,lk," x1,x2,y1,y2,z1,z2:",x1,x2,y1,y2,z1,z2)
        my_mpi_grid.grid.data[li*args.bc_res:(li+1)*args.bc_res,lj*args.bc_res:(lj+1)*args.bc_res,lk*args.bc_res:(lk+1)*args.bc_res] = cell_obj.voxels
    
        local_id += 1
    
  eps = 1e-6
  
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
    
  x = np.arange(my_mpi_grid.bounds[0],my_mpi_grid.bounds[3]+my_mpi_grid.grid.hx-eps,my_mpi_grid.grid.hx)
  y = np.arange(my_mpi_grid.bounds[1],my_mpi_grid.bounds[4]+my_mpi_grid.grid.hy-eps,my_mpi_grid.grid.hy)
  z = np.arange(my_mpi_grid.bounds[2],my_mpi_grid.bounds[5]+my_mpi_grid.grid.hz-eps,my_mpi_grid.grid.hz)
  
 # import binvox_rw
  #with open('skin_160.binvox', 'rb') as f:
 #   model = binvox_rw.read_as_3d_array(f)
  # model.data has boolean 3d array
  #range0 = my_mpi_grid.end_x*args.bc_res - my_mpi_grid.start_x*args.bc_res 
  #my_mpi_grid.grid.data[0:range0,0:samples[1]*args.bc_res,0:samples[2]*args.bc_res] = np.logical_or( model.data[my_mpi_grid.start_x*args.bc_res:my_mpi_grid.end_x*args.bc_res,0:samples[1]*args.bc_res,0:samples[2]*args.bc_res],my_mpi_grid.grid.data[0:range0+1,0:samples[1]*args.bc_res,0:samples[2]*args.bc_res] )
  draw_non_design(nondes_coords, my_mpi_grid.grid.data, my_mpi_grid.bounds[0:3], (my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz), 1)
#   draw_non_design(void_elems, my_mpi_grid.grid.data, my_mpi_grid.bounds[0:3], (my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz), -1)
  from scipy.ndimage import binary_fill_holes
  my_mpi_grid.grid.data = binary_fill_holes(my_mpi_grid.grid.data).astype(bool)
  from pyevtk.hl import gridToVTK  
  gridToVTK("voxels"+str(my_mpi_grid.rank),x,y,z,cellData={"voxels":my_mpi_grid.grid.data.astype(int)})
  
  borders = my_mpi_grid.communicate_edges()
    
  # binary helper array
  shape = np.asarray(my_mpi_grid.grid.data.shape[0:3]) + np.array((2,2,2))
  helper = np.zeros(shape,dtype=bool)
  helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = my_mpi_grid.grid.data
  
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
  # vtk rectgrid is node based
  x = np.arange(my_mpi_grid.bounds[0]-hx,my_mpi_grid.bounds[3]+2*hx,hx)
  y = np.arange(my_mpi_grid.bounds[1]-hy,my_mpi_grid.bounds[4]+2*hy,hy)
  z = np.arange(my_mpi_grid.bounds[2]-hz,my_mpi_grid.bounds[5]+2*hz,hz)
  from pyevtk.hl import gridToVTK  
  gridToVTK("helper"+str(my_mpi_grid.rank),x,y,z,cellData={"helper":helper.astype(int)}) 
 
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  verts, faces, normals, values = measure.marching_cubes(helper,spacing=(np.float32(my_mpi_grid.grid.hx),np.float32(my_mpi_grid.grid.hy),np.float32(my_mpi_grid.grid.hz)),allow_degenerate=False)
  verts = np.asarray(verts) + (my_mpi_grid.grid.hx/2.0,my_mpi_grid.grid.hx/2.0,my_mpi_grid.grid.hx/2.0)
  # translate from (0,0,0) to correct position
  shift = np.asarray(my_mpi_grid.bounds[0:3])
  verts = [p+shift for p in verts]
  
  pd = matviz_vtk.fill_vtk_polydata(verts, faces)
  matviz_vtk.show_write_vtk(pd, 10, "marching"+str(my_mpi_grid.rank)+".vtp")
  
  my_mpi_grid.set_vertices_and_faces(list(verts),list(faces))
  
  my_mpi_grid.gather_data(append=True)
  
  data = None
  if my_mpi_grid.rank == 0:
    import pymesh
    print("before reducing: len(verts):",len(my_mpi_grid.vertices)," len(faces):",len(my_mpi_grid.faces))
    verts = np.asarray(my_mpi_grid.vertices)
    faces = np.asarray(my_mpi_grid.faces)
    verts, faces, info = pymesh.remove_duplicated_vertices_raw(verts,faces,1e-4) 
    verts, faces, info = pymesh.remove_duplicated_faces_raw(verts,faces)
    print("after reducing: len(verts):",len(verts)," len(faces):",len(faces))
    
    pd = matviz_vtk.fill_vtk_polydata(verts, faces)
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(pd)
    clean.Update()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(pd)
    normals.SetConsistency(1)
    normals.SetAutoOrientNormals(1)
    normals.Update()
    matviz_vtk.show_write_vtk(normals.GetOutput(), 10, "marching_all.vtp")
    
    data = (verts,faces)
    
  #sys.exit()  
  
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
  mpi_taubin_smoothing(my_mpi_grid,(xmin,ymin,0,xmax,ymax,zmax))
  
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
    #draw_triangle(t[0], t[1], t[2], grid, value)
    draw_triangle_bresenham(t[0], t[1], t[2], grid, value)

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
  
  for u in np.linspace(0,1,num=np.sqrt(2)*samples10+10):
    for v in np.linspace(0,1,num=np.sqrt(2)*samples20+10):
      if u + v > 1.0 + eps:
        continue
      
      p0 = v0 + u * v1v0 + v * v2v0 + eps
      p = [int(v) for v in p0]
      
      if not idx_out_of_bounds(p, grid.shape):
        grid[tuple(p)] = value

def draw_triangle_bresenham(v0,v1,v2,grid,value=1):
  for trip in [(v0,v1,v2),(v1,v0,v2),(v2,v1,v0)]:
    res01 = bresenham_3d(trip[0],trip[1],grid,value)
    res02 = bresenham_3d(trip[0],trip[2],grid,value)
    
    left = None
    right = None
    if len(res01) > 0 and len(res02) > 0:
      for i in range(max(len(res01),len(res02))):
        if i < len(res01):
          left = res01[i]
        if i < len(res02):
          right = res02[i]
        if left is not None and right is not None:  
          bresenham_3d(left,right,grid,value)
        else:
          print("left right is None:",len(res01),len(res02))

def bresenham_3d(p0,p1,array,value=1):
  x0 = p0[0]
  y0 = p0[1]
  z0 = p0[2]
  
  x1 = p1[0]
  y1 = p1[1]
  z1 = p1[2]
  
  dx = np.abs(x1 - x0)
  sx = 1 if x0 < x1 else -1

  dy = np.abs(y1 - y0) 
  sy = 1 if y0 < y1 else -1

  dz = np.abs(z1 - z0) 
  sz = 1 if z0 < z1 else -1

  # max diff
  dm = np.max([dx,dy,dz])
  i = dm
  x1 = y1 = z1 = int(dm/2)
  
  res = []
  
  while True:
    if not idx_out_of_bounds((x0,y0,z0), array.shape):
      array[x0,y0,z0] = value
    res.append((x0,y0,z0))
    if (i == 0):
      break
    i -= 1
    x1 -= dx
    if x1 < 0:
      x1 += dm
      x0 += sx
      
    y1 -= dy
    if y1 < 0:
      y1 += dm
      y0 +=sy
    
    z1 -= dz
    if z1 <0:
      z1 += dm
      z0 +=sz
      
  return res  
          
        
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
      
      # get bounding box of this tet and fill holes inside 4 drawn triangles
      mindim = [min(0,A[i],B[i],C[i],D[i]) for i in range(3)]
      maxdim = [min(grid.shape[i],max(A[i],B[i],C[i],D[i])) for i in range(3)] 
      
      #subgrid = grid[int(mindim[0]):int(maxdim[0])+1,int(mindim[1]):int(maxdim[1])+1,int(mindim[2]):int(maxdim[2])+1]
      
      ## fill holes does not work, if we have a partial tet at the boundary
      ## solution: add ghost layers with value 1 to 5 faces and keep one to 0
      ## do all combinations and add pictures with logical or
      
      ## need bigger array with ghost layer
      #shape = np.asarray(subgrid.shape[0:3]) + np.array((2,2,2))
      #print("shape:",shape)
      #collect = []
      
      ###### exclude xmin #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmax
      #helper[shape[0]-1,1:shape[1],1:shape[2]] = 1
      ## layer at ymin
      #helper[1:shape[0],0,1:shape[2]] = 1
      ## layer at ymax
      #helper[1:shape[0],shape[1]-1,1:shape[2]] = 1
      ## layer at zmin
      #helper[1:shape[0],1:shape[1],0] = 1
      ## layer at zmax
      #helper[1:shape[0],1:shape[1],shape[2]-1] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
      
      ###### exclude xmax #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmin
      #helper[0,1:shape[1],1:shape[2]] = 1
      ## layer at ymin
      #helper[1:shape[0],0,1:shape[2]] = 1
      ## layer at ymax
      #helper[1:shape[0],shape[1]-1,1:shape[2]] = 1
      ## layer at zmin
      #helper[1:shape[0],1:shape[1],0] = 1
      ## layer at zmax
      #helper[1:shape[0],1:shape[1],shape[2]-1] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
      
      ####### exclude ymin #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmin
      #helper[0,1:shape[1],1:shape[2]] = 1
      ## layer at xmax
      #helper[shape[0]-1,1:shape[1],1:shape[2]] = 1
      ## layer at ymax
      #helper[1:shape[0],shape[1]-1,1:shape[2]] = 1
      ## layer at zmin
      #helper[1:shape[0],1:shape[1],0] = 1
      ## layer at zmax
      #helper[1:shape[0],1:shape[1],shape[2]-1] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
      
      ###### exclude ymax #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmin
      #helper[0,1:shape[1],1:shape[2]] = 1
      ## layer at xmax
      #helper[shape[0]-1,1:shape[1],1:shape[2]] = 1
      ## layer at ymin
      #helper[1:shape[0],0,1:shape[2]] = 1
      ## layer at zmin
      #helper[1:shape[0],1:shape[1],0] = 1
      ## layer at zmax
      #helper[1:shape[0],1:shape[1],shape[2]-1] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
      
      ###### exclude zmin #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmin
      #helper[0,1:shape[1],1:shape[2]] = 1
      ## layer at xmax
      #helper[shape[0]-1,1:shape[1],1:shape[2]] = 1
      ## layer at ymin
      #helper[1:shape[0],0,1:shape[2]] = 1
      ## layer at ymax
      #helper[1:shape[0],shape[1]-1,1:shape[2]] = 1
      ## layer at zmax
      #helper[1:shape[0],1:shape[1],shape[2]-1] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
      
      ###### exclude zmax #########
      #helper = np.zeros(shape,dtype=bool)
      #helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = subgrid
      ## fill layer at xmin
      #helper[0,1:shape[1],1:shape[2]] = 1
      ## layer at xmax
      #helper[shape[0]-1,1:shape[1],1:shape[2]] = 1
      ## layer at ymin
      #helper[1:shape[0],0,1:shape[2]] = 1
      ## layer at ymax
      #helper[1:shape[0],shape[1]-1,1:shape[2]] = 1
      ## layer at zmin
      #helper[1:shape[0],1:shape[1],0] = 1
      #collect.append(binary_fill_holes(helper).astype(bool))
       
      #new = np.zeros(shape,dtype=bool)
      #for c in collect:
        #new = np.logical_or(new,c)
        
      #grid[int(mindim[0]):int(maxdim[0])+1,int(mindim[1]):int(maxdim[1])+1,int(mindim[2]):int(maxdim[2])+1] = new[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1]  
      
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
      subgrid = binary_fill_holes(subgrid).astype(int)
      
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
class MPI_Grid():
  # total_samples: list with number of total samples in 3 directions
  # bc_res: resolution of one base cell (usually 40)
  # bounds: global (xmin,ymin,zmin,xmax,ymax,zmax)
  def __init__(self, comm, total_samples,bc_res,bounds):
    assert(len(total_samples) == 3)
    assert(len(bounds) == 6)

    self.comm = comm    
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
      data = self.comm.gather((self.vertices, self.faces),root=root)
      if self.rank == root:
        for d in data:
          offset = len(self.vertices)
          self.vertices.extend(d[0])
          self.faces.extend(list(np.asarray(d[1])+offset))
    else: # update
      data = self.comm.gather((self.start_verts_idx, self.end_verts_idx,self.vertices[self.start_verts_idx:self.end_verts_idx]),root=root)
      #print("rank:", self.rank," len(vertices):", len(self.vertices))
      if self.rank == root:
        self.update_vertices(data)
        
  def update_vertices(self,recv):
    for r in recv:
      if r is not None:
        start = r[0]
        end = r[1]
        assert(r[2] is not None and r[2] is not None)
        self.vertices[start:end] = r[2]
        
class RectGrid():
  def __init__(self,nx,ny,nz,bounds):
    self.data = np.full((nx,ny,nz),999,dtype=bool)

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

def mpi_taubin_smoothing(mpi_grid,bounds=None):
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
  while iter < 20:
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
    

