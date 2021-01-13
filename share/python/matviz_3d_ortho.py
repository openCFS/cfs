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

found_pymesh = False

# similar to create_3d_cross_ip; # without rotation and shearing
def create_3d_interpretation_ortho(args,reg_info,barycenters,min_bb,max_bb,design,scale,samples,grad,nondes=None):
  # args: options for basecell, e.g. voxel resolution for local microstructure, interpolation type, beta, eta, ...
  # coords, design, angles: element center coordinates and design values s1,s2,s3,angle per finite element
  # min_bb/max_bb: bounding box of design regions
  # scale: parameter for scaling the cell size if necessary
  # nondes: store info on solid and void nondesign regions, has 2 entries: 0 -> solid non-design, 1 -> void non_design
  # nondes[0]: solid nondesign -> (centers, min_bb, max_bb, elem_dim)
  # nondes[0][centers]: list of elements (corner vertices) that define non-design regions

  s1 = design['s1']
  s2 = design['s2']
  s3 = design['s3']


  # MPI_Init() or MPI_Init_thread() is actually called when you import the MPI
  # use the standard communicator
  comm = MPI.COMM_WORLD

  # point coordinates from h5 file
  centers, _, _ = barycenters[0:3]

  if scale <= 0:
    scale = 1.0

  # order: min_x,min_y,min_z,max_x,max_y,max_z
  design_bounds = np.ones(6) * (-1)
  design_bounds[0:3] = min_bb[0:3]
  design_bounds[3:6] = max_bb[0:3]

#   print("design_bounds:",design_bounds)
  # 0:xmin,1:ymin,2:zmin,3:xmax,4:ymax,5:zmax
  h_des = (abs(design_bounds[3] - design_bounds[0]), abs(design_bounds[4] - design_bounds[1]), abs(design_bounds[5] - design_bounds[2]))

  # set size dx/dy/dz of one cell
  dx_des = h_des[0] / samples[0]
  dy_des = h_des[1] / samples[1]
  dz_des = h_des[2] / samples[2]

  min_thresh = 0.153
  max_thresh = 0.98

  thresh = None
  if args.bc_thresh is not None:
    thresh = [float(t) for t in args.bc_thresh.split(',')]
    assert(len(thresh) == 2)

#   print("samples:",samples)
#   print("h_des:",h_des)
#   print("design_bounds:",design_bounds)
#   print("dx_des,dy_des,dz_des:",dx_des,dy_des,dz_des)
#   print("nx,ny,nz:",nx,ny,nz)

  data_grid, sample_coords = matviz_vtk.get_interpolation_natural_neighbor(barycenters, s1, s2, s3, dx_des, dy_des, dz_des)

  s1 = comm.bcast(s1,root=0)
  s2 = comm.bcast(s2,root=0)
  s3 = comm.bcast(s3,root=0)

  my_mpi_grid = MPI_Grid(comm)
  
  if my_mpi_grid.rank == 0:
    try:
      import pymesh
      global found_pymesh
      found_pymesh = True
    except:
      print("pymesh not found!")

  ################# voxel world #########################################
  nondes_min = 99999
  nondes_max = -99999
  nondes_coords = None
  holes = None
  design_elems = None
  if my_mpi_grid.rank == 0 and  nondes is not None: # assume holes are within bounding box of solid non-design and design domain
    if nondes[0] is not None:
      nondes_coords = nondes[0][0]
      nondes_min = nondes[0][1]
      nondes_max = nondes[0][2]
      assert(len(nondes_min) == 3)
      assert(len(nondes_max) == 3)
    if nondes[1] is not None:
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

  print("min_bb design:" + str(min_bb))
  print("max_bb design:" + str(max_bb))
  print("min:"+ str(bounds[0:3]))
  print("max:"+ str(bounds[3:6]))

  my_mpi_grid.init_data_grid(samples, args.bc_res, bounds)
  my_mpi_grid.to_info()

  for k in range(samples[2]):
    for j in range(samples[1]):
      for i in range(my_mpi_grid.start_x,my_mpi_grid.end_x):
        li = i - my_mpi_grid.start_x

        # data grid contains interpolated data for three faces of a cube: y-z, x-z and x-y face
        x1 = data_grid[0][i,j,k][0]
        x2 = data_grid[0][i+1,j,k][0]
        y1 = data_grid[1][i,j,k][1]
        y2 = data_grid[1][i,j+1,k][1]
        z1 = data_grid[2][i,j,k][2]
        z2 = data_grid[2][i,j,k+1][2]

        all_values = [x1,x2,y1,y2,z1,z2]

        # handle void and solid cells
        if thresh is not None:
          if (np.array(all_values) > thresh[1]).sum() >= 3:
#             print("found solid cell: ","rank:",my_mpi_grid.rank," global i,j,k:",i,j,k, " local:",li,j,k," x1,x2,y1,y2,z1,z2:",x1,x2,y1,y2,z1,z2," coords: ",sample_coords[i,j,k])
            print("found solid cell: rank:" + str(my_mpi_grid.rank) + " global i,j,k:" + str([i,j,k]) + " x-direction: " + str(li*args.bc_res) + ":" + str((li+1)*args.bc_res) + " x1,x2,y1,y2,z1,z2:" + str([x1,x2,y1,y2,z1,z2]))
            my_mpi_grid.grid.data[li*args.bc_res:(li+1)*args.bc_res,j*args.bc_res:(j+1)*args.bc_res,k*args.bc_res:(k+1)*args.bc_res] = 1
            continue
          # check if we have an interface between lattice and void
          # in this case, t < thresh[0] is True but we don't want to cut off this cell
          # but round up to thresh[0]
          elif any(t < thresh[0] for t in all_values):
            if (np.array(all_values) < thresh[0]).sum() > 3:
  #             print("found void cell: ","rank:",my_mpi_grid.rank," global i,j,k:",i,j,k, " local:",li,j,k," x1,x2,y1,y2,z1,z2:",x1,x2,y1,y2,z1,z2," coords: ",sample_coords[i,j,k])
              print("found void cell: rank:" + str(my_mpi_grid.rank) +" global i,j,k:" + str([i,j,k]) + " x1,x2,y1,y2,z1,z2:" + str([x1,x2,y1,y2,z1,z2]))
              continue # skip void cell
            else:
              print("detected interface s/v rank:" + str(my_mpi_grid.rank) +" global i,j,k:" + str([i,j,k]) + " x1,x2,y1,y2,z1,z2:" + str([x1,x2,y1,y2,z1,z2]))
              all_values = [np.maximum(t,thresh[0]) for t in all_values]
              x1, x2, y1, y2, z1, z2 = all_values
              print("rounded to ", all_values)


        flags = None
        bc_input  = basecell.Basecell_Data(args.bc_res,args.bc_bend,x1,x2,y1,y2,z1,z2,args.bc_interpolation,args.bc_beta,args.bc_eta,target="volume_mesh",bc_flags=flags)
        bc_input.eta = 0.7
        bc_input.stiffness_as_diameter = True
        cell_obj = basecell.Basecell(bc_input,id)
        print("rank:"+ str(my_mpi_grid.rank) + " global i,j,k:" + str([i,j,k]) + " x1,x2,y1,y2,z1,z2:" + str([x1,x2,y1,y2,z1,z2]))
        my_mpi_grid.grid.data[li*args.bc_res:(li+1)*args.bc_res,j*args.bc_res:(j+1)*args.bc_res,k*args.bc_res:(k+1)*args.bc_res] = cell_obj.voxels
        sys.stdout.flush()

  eps = 1e-6

  
  hx = (my_mpi_grid.bounds[3]-my_mpi_grid.bounds[0])/my_mpi_grid.grid.nx
  hy = (my_mpi_grid.bounds[4]-my_mpi_grid.bounds[1])/my_mpi_grid.grid.ny
  hz = (my_mpi_grid.bounds[5]-my_mpi_grid.bounds[2])/my_mpi_grid.grid.nz

  tmp = np.zeros_like(my_mpi_grid.grid.data,dtype=int)
  draw_non_design(design_elems, tmp, my_mpi_grid.bounds,(my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=True)
  my_mpi_grid.grid.data *= tmp
  
  netskin = compute_net_skin((my_mpi_grid.grid.nx,my_mpi_grid.grid.ny,my_mpi_grid.grid.nz), tmp, my_mpi_grid.grid.data)
  
  # debugging
  x = np.arange(0,my_mpi_grid.grid.nx+1,1)
  y = np.arange(0,my_mpi_grid.grid.ny+1,1)
  z = np.arange(0,my_mpi_grid.grid.nz+1,1)

  from pyevtk.hl import gridToVTK
  gridToVTK("netskin",x,y,z,cellData={"netskin":netskin})
  
  my_mpi_grid.grid.data = np.clip(netskin+my_mpi_grid.grid.data,0,1)
  
  vol = comm.gather(np.sum(my_mpi_grid.grid.data),root=0)
  if my_mpi_grid.rank == 0:
    print("design volume:" + str(np.sum(vol) / (samples[0]*samples[1]*samples[2]*args.bc_res**3)))
      
  if nondes:    
    design_elems = None
    if nondes_coords is not None:
      draw_non_design(nondes_coords, my_mpi_grid.grid.data, my_mpi_grid.bounds,(my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=True)
    if holes is not None:
      draw_non_design(holes, my_mpi_grid.grid.data, my_mpi_grid.bounds, (my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz),solid=False)


    vol = comm.gather(np.sum(my_mpi_grid.grid.data),root=0)
    if my_mpi_grid.rank == 0:
      print("volume with nondes:" + str(np.sum(vol) / (samples[0]*samples[1]*samples[2]*args.bc_res**3)))
    nondes_coords = None
    holes = None
  else:
    vol = comm.gather(np.sum(my_mpi_grid.grid.data),root=0)
    if my_mpi_grid.rank == 0:
      print("volume:" + str(np.sum(vol) / (samples[0]*samples[1]*samples[2]*args.bc_res**3)))


  borders = my_mpi_grid.communicate_edges()

  # binary helper array
  shape = np.asarray(my_mpi_grid.grid.data.shape[0:3]) + np.array((2,2,2))
  helper = np.zeros(shape,dtype=bool)
  helper[1:shape[0]-1,1:shape[1]-1,1:shape[2]-1] = my_mpi_grid.grid.data

#   x = np.arange(bounds[0],bounds[3]+hx,hx)
#   y = np.arange(bounds[1],bounds[4]+hy,hy)
#   z = np.arange(bounds[2],bounds[5]+hz,hz)
#
#   from pyevtk.hl import gridToVTK
#   gridToVTK("image",x,y,z,cellData={"data":my_mpi_grid.grid.data})

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

  print("hx,hy,hz:" + str([my_mpi_grid.grid.hx,my_mpi_grid.grid.hy,my_mpi_grid.grid.hz]))
  hx = my_mpi_grid.grid.hx
  hy = my_mpi_grid.grid.hy
  hz = my_mpi_grid.grid.hz

  #print(helper.shape,np.sum(helper))
  from skimage import measure
  # coords of vertices lie in [0,1-h]
  normals = []
  my_mpi_grid.vertices = []
  my_mpi_grid.faces = []
#   marching_cubes(helper,(np.float32(my_mpi_grid.grid.hx),np.float32(my_mpi_grid.grid.hy),np.float32(my_mpi_grid.grid.hz)),my_mpi_grid.vertices, my_mpi_grid.faces,normals)
  my_mpi_grid.vertices, my_mpi_grid.faces, _, _ = measure.marching_cubes_lewiner(helper,spacing=(np.float32(my_mpi_grid.grid.hx),np.float32(my_mpi_grid.grid.hy),np.float32(my_mpi_grid.grid.hz)),allow_degenerate=False,step_size=1)
  # translate from (0,0,0) to correct position
  # and marching cubes shift everything by 0.5 * hx/hy/hz
  shift = np.asarray(my_mpi_grid.bounds[0:3]) - np.asarray((my_mpi_grid.grid.hx/2.0,my_mpi_grid.grid.hy/2.0,my_mpi_grid.grid.hz/2.0))
  print("shift:" + str(shift))
  my_mpi_grid.vertices = [p+shift for p in my_mpi_grid.vertices]

  # hope for python's garbage collector to delete voxel array
  helper = None

#   pd = matviz_vtk.fill_vtk_polydata(my_mpi_grid.vertices, my_mpi_grid.faces)
#   matviz_vtk.show_write_vtk(pd, 10, "marching"+str(my_mpi_grid.rank)+".vtp")

  pd = None

  # communicate all pieces to root for global vertices list
  my_mpi_grid.gather_data(append=True)

  sys.stdout.flush()
  data = None
  if my_mpi_grid.rank == 0:
#     if found_pymesh:
#       print("before reducing: len(verts):" + str(len(my_mpi_grid.vertices)) + " len(faces):" + str(len(my_mpi_grid.faces)))
#       sys.stdout.flush()
#    
#       my_mpi_grid.vertices, my_mpi_grid.faces, info = pymesh.remove_duplicated_vertices_raw(np.asarray(my_mpi_grid.vertices),np.asarray(my_mpi_grid.faces),1e-4)
#       my_mpi_grid.vertices, my_mpi_grid.faces, info = pymesh.remove_duplicated_faces_raw(np.asarray(my_mpi_grid.vertices),np.asarray(my_mpi_grid.faces))
#       print("after reducing: len(verts):" + str(len(my_mpi_grid.vertices)) + " len(faces):" + str(len(my_mpi_grid.faces)))
#       sys.stdout.flush()

    # add vertex id to list of vertices, need it later on when scattering
    # data for smoothing

    pd = matviz_vtk.fill_vtk_polydata(my_mpi_grid.vertices, my_mpi_grid.faces)
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(pd)
    normals.SetConsistency(1)
    normals.SetAutoOrientNormals(1)
    normals.Update()
    sys.stdout.flush()
    vtpName = args.save+".vtp" if args.save is not None else "interpreted.vtp"
    matviz_vtk.show_write_vtk(pd, 10, vtpName)

#     data = (verts,faces)
  if args.bc_smooth > 0:
  # #   # broadcast all verts to all ranks
  #   data = my_mpi_grid.comm.bcast(data,root=0)
  #   my_mpi_grid.set_vertices_and_faces(data[0],data[1])
    overlap_offsets = None
    # on which chunks each rank should work on
    chunks = None
    offsets = None
    # how many elements each rank gets
    rankcounts = None
    displ = None
    recvbuf = None

    if my_mpi_grid.rank == 0:
      my_mpi_grid.update_connectivity_verts_faces()

      chunks, offsets = my_mpi_grid.calc_vertices_chunks()
      print(" rank 0 offsets:" + str(offsets))
      new_v = []
      for i,v in enumerate(my_mpi_grid.vertices):
        new_v.append((v[0],v[1],v[2],i))

      # create dict to make vertex id independent from list index
      my_mpi_grid.vertices = np.array(new_v)

  #     print("\n offsets")
  #     for p in range(my_mpi_grid.size):
  #       print("rank: start=",vertices_offsets[p][0], " end=",vertices_offsets[p][1])
      # broadcast chunks of vertices to ranks
      # make sure ranks overlap to ensure smoothing at the rank boundaries
      # add number of voxels of two voxel layers in +x/-x direction to chunks
      # e.g. one voxel layer has 1 x 14 x 9 = 126 voxels, add 126 vertices from each side
      # (left and right) to the rank's chunk

      add = 1.2 * samples[1] * samples[2] * args.bc_res**2
      rankcounts = [None] * len(chunks)
      displ = []
      for p in range(my_mpi_grid.size):
        # first is start_idx, second is end_idx
        start = offsets[p][0]
        end =  offsets[p][1]
        if  p != 0:
          start -= add
          print("rank: " + str(my_mpi_grid.rank) + " start: " + str(start) + " end: " + str(end) + " offsets: " + str(offsets))
          assert(start >= 0)
        if p != my_mpi_grid.size - 1:
          end += add

        displ.append(start)
        rankcounts[p] = end-start
  #       print("\noverlapping offsets:")
  #       print(" rank ", p, ": start=",start," end=",end)

      # need to flatten array with vertices for communication
      displ = np.array(displ,dtype=int) * 4
      rankcounts = np.array(rankcounts,dtype=int) * 4

      str(("displ:" + str(np.array(displ)/4)))

    # stuff we're working on
    offsets = my_mpi_grid.comm.bcast(offsets,root=0)
    chunks = my_mpi_grid.comm.bcast(chunks,root=0)

    # extended stuff that we need to work
    displ = my_mpi_grid.comm.bcast(displ,root=0)
    rankcounts = my_mpi_grid.comm.bcast(rankcounts,root=0)

    recvbuf = np.empty(int(rankcounts[my_mpi_grid.rank]))
    my_mpi_grid.comm.Scatterv(sendbuf=[my_mpi_grid.vertices,rankcounts,displ,MPI.DOUBLE],recvbuf=recvbuf,root=0)
    my_mpi_grid.connectivity = my_mpi_grid.comm.bcast(my_mpi_grid.connectivity,root=0)


    # convert to dict with vertex id as key and coordinates as value
    my_mpi_grid.vertices = convert_verts_list_to_dict(np.reshape(recvbuf, (int(rankcounts[my_mpi_grid.rank]/4),4)))
    #if my_mpi_grid.rank == 1:
      #print("rank ", my_mpi_grid.rank, " verts:",np.reshape(recvbuf, (int(rankcounts[my_mpi_grid.rank]/4),4)))
      #print(len(recvbuf)/4)
      #print("verts:",my_mpi_grid.vertices)
    #sys.exit()
    my_mpi_grid.start_verts_idx = int(offsets[my_mpi_grid.rank][0])
    my_mpi_grid.end_verts_idx = int(offsets[my_mpi_grid.rank][1])

  #   print("\nrank ", my_mpi_grid.rank, " recieved ", len(my_mpi_grid.vertices)," overlapping data from id ", displ[my_mpi_grid.rank]/4," to" , displ[my_mpi_grid.rank]/4+rankcounts[my_mpi_grid.rank]/4, " working on data starting at", offsets[my_mpi_grid.rank], " chunks:",chunks)

    for v in my_mpi_grid.vertices:
      assert(v is not None)

    xmin = min(my_mpi_grid.vertices.values(), key=lambda t: t[0])[0]
    xmax = max(my_mpi_grid.vertices.values(), key=lambda t: t[0])[0]
    ymin = min(my_mpi_grid.vertices.values(), key=lambda t: t[1])[1]
    ymax = max(my_mpi_grid.vertices.values(), key=lambda t: t[1])[1]
  #   zmin = min(my_mpi_grid.vertices, key=lambda t: t[2])[2]
    zmax = max(my_mpi_grid.vertices.values(), key=lambda t: t[2])[2]


    # do parallel smoothing here
    mpi_taubin_smoothing(my_mpi_grid,(xmin,ymin,0,xmax,ymax,zmax),niter=args.bc_smooth)

    # send smoothed data to root
    my_mpi_grid.gather_data(append=False, root=0)

    if my_mpi_grid.rank != 0:
      sys.exit()

    pd = matviz_vtk.fill_vtk_polydata(list(my_mpi_grid.vertices.values()),my_mpi_grid.faces,pointIds=list(my_mpi_grid.vertices.keys()))
    matviz_vtk.show_write_vtk(pd, 10, "smoothed_marching"+str(my_mpi_grid.rank)+".vtp")

  if my_mpi_grid.rank != 0:
      sys.exit()

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
  nx, ny, nz, _, _ = array.shape

  if idx[0] >= nx or idx[0] < 0  or idx[1] >= ny or idx[1] < 0  or idx[2] >= nz  or idx[2] < 0 :
    return None
  else:
    if array[idx][0] == -1 or array[idx][1] == -1 or array[idx][2] == -1:
      print("fallback")
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

# @param shape: extents of 3d image
# @param trimmer: domain used for trimming the design
# @param trimmed: result of trimming operation
def compute_net_skin(shape, trimmer, trimmed):
  from scipy import ndimage
  from skimage.segmentation import find_boundaries
  assert(len(shape) == 3)
  rows, cols, depth = shape
  
  import cv2
  #erode solid part
  kernel = np.ones((31,31),dtype=np.uint8)
  eroded = cv2.erode(np.uint8(trimmer),kernel,iterations = 1)
  #eroded = ndimage.binary_erosion(trimmer,structure=np.ones((int(0.5*rows), int(0.5*cols), int(0.5*depth))))
  inv_eroded = (1 - eroded).astype(bool)
  # throw away pixels that are not on boundary of trimmed structure
  candidates = (trimmed * inv_eroded).astype(bool)
  # boundaries of trimmer
  perim = find_boundaries(trimmer)
  
  directions = [np.array([1,0,0]), np.array([-1,0,0]), np.array([0,1,0]), np.array([0,-1,0]), np.array([0,0,1]), np.array([0,0,-1])]
  #directions = [np.array([1,0]), np.array([0,1])]
  n_max = 10
  
  projected = []
  for vec in directions:
    for r in range(rows):
      for c in range(cols):
        for d in range(depth):
          for n in range(1,n_max+1):
            # project solid voxels along unit vectors
            if candidates[r,c,d] == 1:
              projected.append(np.array([r,c,d]) + n*vec)
              #print("n,[r,c],vec,projected[-1]:",n,np.array([r,c]),vec,projected[-1])
              #sys.exit()
              
  netskin = np.zeros(shape,dtype=np.uint8)
  
  for p in projected:
    t = tuple(p)
    if p[0] >= trimmed.shape[0] or p[1] >= trimmed.shape[1] or p[2] >= trimmed.shape[2]:
      continue
    if perim[t] == 1:
      #print("trimmed[tuple(p)]:",trimmed[tuple(p)])
      netskin[t] = 1
      
  return netskin       

   
# @param elems: list of lists, each entry contains 4 vertices (cartesian) of a tetra-/hexahedron
# @param grid: where to draw
# @param bounds: list of bounds of cartesian world
# @param h: lattice spacings in 3 directions
# @param solid or void non-design?
# @param toarray: return np array with non-design?
def draw_non_design(elems,grid,bounds,h,solid=True,toarray=False):
  assert(len(bounds) >=3 and len(h) == 3)
  assert(elems is not None)
  ug = vtk.vtkUnstructuredGrid()
  vtkpoints = vtk.vtkPoints()
  for e in elems:
    #print("e:",e)
    assert(len(e) == 4 or len(e) == 8)
    vtkObj = vtk.vtkTetra() if len(e) == 4 else vtk.vtkHexahedron()
    for i,p in enumerate(e):
      id = vtkpoints.InsertNextPoint(p)
      vtkObj.GetPointIds().SetId(i,id)

    ug.InsertNextCell(vtkObj.GetCellType(),vtkObj.GetPointIds())

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
  res = np.zeros(grid.shape) if toarray else None
  for p in prod:
    # transfer to voxel world
    i,j,k = draw_profile_functions.cartesian_to_voxel_coords(p,bounds[0],bounds[1],bounds[2],h[0],h[1],h[2])
    assert(not idx_out_of_bounds((i,j,k),grid.shape))
    grid[i,j,k] = solid
    if toarray:
      res[i,j,k] = solid 
    
  if toarray:
    assert(res is not None)
    return res  
    

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
    self.connectivity = None
    self.start_verts_idx = None
    self.end_verts_idx = None

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

  # compute which rank gets how many vertices for smoothing
  # and what the offset of the global vertices array is
  def calc_vertices_chunks(self):
    assert(self.rank == 0)
    assert(self.vertices is not None and len(self.vertices) > 0)
    num_verts = len(self.vertices)
    chunks = [int(num_verts/self.size) + (1 if p < num_verts%self.size else 0) for p in range(0,self.size)]
    offsets = []
    for p in range(self.size):
      # start idx of 'vertices' list
      start_idx = int(np.sum(chunks[0:p]))
      end_idx = int(start_idx + chunks[p])
      offsets.append((start_idx,end_idx))

    return chunks, offsets

  def to_info(self):
    print("---- mpi distr grid ----")
    print("rank:" + str(self.rank))
    print("bounds:" + str(self.bounds))
    print("start:" + str(self.start_x) + " end:" + str(self.end_x))
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
    recvkey = []
    recvval = []
    for direction in [-1,1]:
      # coordinate dimension of shift is always 1, because we only have one dimension
      # coordinate dimension 1-based
      source, dest = self.cart.Shift(1,direction)

      key = list(self.vertices.keys()) # vertex ids
      values = list(self.vertices.values()) # vertex coords

      tmp_key = self.cart.sendrecv(key,dest=dest,source=source)
      tmp_values = self.cart.sendrecv(values,dest=dest,source=source)
#       print("rank:",self.rank," got data from ",source)

      if source != MPI.PROC_NULL:
        recvkey.extend(tmp_key)
        recvval.extend(tmp_values)

#     print("rank:",self.rank," got ",len(recv)," data chunks")
    return [recvkey,recvval]

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
      keysbuf = None
      coordsbuf = None
      keys = list(self.vertices.keys())
      coords = list(self.vertices.values())
      # find out where to truncate the list
      start_key_idx = keys.index(self.start_verts_idx)
      keys = keys[start_key_idx:]
      coords = coords[start_key_idx:]

      print("rank " + str(self.rank) + " sends " + str(len(keys)) + " vertices to root \n")
      lenkeys = self.comm.gather(len(keys),root=root)
      lencoords = self.comm.gather(len(coords)*3,root=root)

      if self.rank == root:
        keysbuf = np.empty(int(np.sum(lenkeys)),dtype=int)
        coordsbuf = np.empty(int(np.sum(lencoords)))

      self.comm.Gatherv(sendbuf=np.array(keys),recvbuf=(keysbuf,lenkeys),root=root)
      self.comm.Gatherv(sendbuf=np.array(coords),recvbuf=(coordsbuf,lencoords),root=root)

      if self.rank == root:
        print("lenkeys:" + str(lenkeys))
        print("lencoords:" + str(lencoords))
        print("int(np.sum(lencoords)/3):"+ str(int(np.sum(lencoords)/3)))
        coordsbuf = np.reshape(coordsbuf,(int(np.sum(lencoords)/3),3))
        print("root recieved " + str(np.sum(lenkeys))," vertex ids \n")
        if coordsbuf.shape[0] != np.sum(lenkeys):
          print("coordsbuf.shape[0] != np.sum(lenkeys)\n")
          print(str(coordsbuf.shape[0]) + "," + str(np.sum(lenkeys)))
          print("lenkeys:"+ str(lenkeys))
          sys.exit()

        for i,k in enumerate(keysbuf):
          self.vertices[k] = coordsbuf[i]

  def update_vertices(self,recv):
    key, values = recv
    for i,k in enumerate(key):
      if k in self.vertices:
        self.vertices[k] = values[i]

class RectGrid():
  def __init__(self,nx,ny,nz,bounds):
    self.data = np.zeros((nx,ny,nz),dtype=int)

    self.nx = nx
    self.ny = ny
    self.nz = nz

    width = [bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]]
    self.hx = width[0] / nx
    self.hy = width[1] / ny
    self.hz = width[2] / nz

  def to_info(self):
    print("res:" + str([self.nx,self.ny,self.nz]))
    print("spacing:" + str([self.hx,self.hy,self.hz]))
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

#   print("rank ", mpi_grid.rank," keys:",mpi_grid.vertices.keys())
#   sys.exit()
  print("niter:" + str(niter))

  while iter < niter:
    sys.stdout.flush()
    old_points = mpi_grid.vertices.copy()
    assert(max(max(mpi_grid.connectivity)) < len(mpi_grid.vertices))
    #if mpi_grid.rank == 1:
      #print("\n before smooth verts:",mpi_grid.vertices)
    # shrink
    mpi_grid.vertices = laplacian_smoothing_dict(mpi_grid.vertices,mpi_grid.connectivity,lamb,start=start,end=end,rank=mpi_grid.rank)
    # communicate smoothed vertices to other ranks
    mpi_grid.update_vertices(mpi_grid.communicate_vertices())

    # expand
    mpi_grid.vertices = laplacian_smoothing_dict(mpi_grid.vertices,mpi_grid.connectivity,-lamb-0.04,start=start,end=end,rank=mpi_grid.rank,bounds=bounds)
    assert(mpi_grid.vertices is not None)
    # communicate smoothed vertices to other ranks
    mpi_grid.update_vertices(mpi_grid.communicate_vertices())

    #res = basecell.residual(old_points[start:end], mpi_grid.vertices[start:end])
    #print("rank:", mpi_grid.rank, " iter:", iter, " residual:", res)
    print("rank:" + str(mpi_grid.rank) + " iter:" + str(iter))
    iter += 1

  print("Taubin smoothing with " + str(iter) + " iterations and res="  + str(res))

# @param point: list with 3 components
# @param bounds: list with 6 entries (xmin,ymin,zmin,xmax,ymax,zmax)
def out_of_bounds(point,bounds):
  eps = 1e-6
  if bounds[0]-eps <= point[0] <= bounds[3]+eps and bounds[1]-eps <= point[1] <= bounds[4]+eps and bounds[2]-eps <= point[2] <= bounds[5]+eps:
    return False
  else:
    print("point "  + str(point) + " is out of bounds "  + str(bounds))
    return True

# convert list to dict with list idx as key and list content as values
def convert_verts_list_to_dict(ls):
  return {int(v[3]): v[0:3] for v in ls}

# same as laplacian_smoothing from basecell.py, but vertices are stored in a dict
# and vertex ids are not contiguous
# laplacian smoothing: p_i = p_i + \lambda * L(p_i)
# using weighted average: L(p_i) = (w_ij*p_j + w_ik*p_k) / (w_ik+w_ik) - p_i, assuming neighbors are p_j,p_k
# bounds: tuple/list with 6 entries - points on these boundaries are not smoothed
# bounds order: xmin,ymin,zmin,xmax,ymax,zmax
def laplacian_smoothing_dict(points,connectivity,lamb,start=0,end=None,bounds=None,rank=None):
  if end == None:
    end = len(points)
  new_points = dict()
  if bounds is None:
    bounds = [-9999,-9999,-9999,9999,9999,9999] # default unit cube

  #print("rank:",rank," smoothing start:",start," end:",end, " lambda:",lamb)
  for id, coords in points.items():
    p = coords

    if id < start or id >= end:
      new_points[id] = p
      continue
    if np.isclose(p[0], bounds[0]) or np.isclose(p[1], bounds[1]) or np.isclose(p[2], bounds[2]) or np.isclose(p[0], bounds[3]) or np.isclose(p[1], bounds[4]) or np.isclose(p[2], bounds[5]):
      new_points[id] = p
      assert(new_points[id] is not None)
    else:
      # calculate L(p_i)
      # w_ij*p_j + w_ik*p_k
      L = 0
      # nid is id of a neighbor node
      neighborhood = connectivity[id]
      #print("neighborhood: ",neighborhood)
      for nid in neighborhood:
        if nid not in points:
          print("rank :" + str(rank) + " point "+ str(id) + " has neighbors " + str(neighborhood) + " vertex with id " + str(nid) + " is not in points list!\n")#keys:",points.keys())
          if rank == 1:
            print("number of verts:" + str(len(points)))
            print(points)
          sys.exit()

        # n are coords of neighbor with id nid
        n = np.asarray(points[nid])
        assert(n is not None)
        # distance between neighbor and this node
        w = 1.0 / np.linalg.norm(len(neighborhood))
        assert(p is not None)
        L += w * (n-p)
        #print("w:",w," n:",n," p:",p)
        #print("n-p:",n-p," w*(n-p):",w*(n-p)," L:",L)

      new_points[id] = p + lamb * L
      #print(p,"+",lamb,"*",L)
      #print("old:",points[i]," new:",new_points[i])
      assert(new_points[id] is not None)

  assert(len(new_points) == len(points))
#   for i,p in enumerate(new_points):
#     if p is None:
#       print("Rank:",rank," i:",i, " has None")
#     assert(p is not None)

  #for i in range(len(points)):
  #  print("old:",points[i]," new:",new_points[i])

  return new_points

