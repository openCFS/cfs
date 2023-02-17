#!/usr/bin/env python

# while mesh_tool.py collects base features for creating meshes, and create_mesh.py provides genertion of basic/generic meshes, 
# this create_special_mesh.py creates common special meshes. Note that purely project related meshes shall not go here!

# you can also import this file and use the functions
from mesh_tool import *
import argparse

# sets the outer elements to the provided region name. This is a layer of one element
# works for 2D and 3D
def set_ghost_region(mesh, name = 'ghost'):
  for i, e in enumerate(mesh.elements):
    x, y, z = mesh.element_idx_by_pos(i, -42)
    if x == 0 or y == 0 or z == 0:
      e.region = name
    if x == mesh.nx-1 or y == mesh.ny-1 or z == mesh.nz-1:
      e.region = name

# set the region for a centered rect or circ inclusion.
# works for 2D and 3D!
# @param type is 'rect' for box or cube or 'circ' for circle or ball
# @param size is in meters
# returns the number of changed elements
def set_centered_inclusion(mesh, type, size, name = 'inclusion'):
  assert type in ['rect','circ']
  assert size is not None
  center = calc_barycenter(mesh) 

  # the difference is just max-norm or Euclidian norm
  order = np.inf if type == 'rect' else 2 

  cnt = 0
  
  for e in mesh.elements:
    test = mesh.calc_barycenter(e)
    if np.linalg.norm(test-center,ord=order) <= size / 2:
      e.region = name
      cnt += 1   
  
  return cnt  

# set a top layer of given thickness in meter to given name.
# works in 2D and 3D where in both cases top is for the y coordinate
# can be used as reinforcment
# @name note that the name must not be 'top' as this is already used for nodes and
#       cfs cannot handle same name for different entity types 
def set_top_layer(mesh, thickness, name = 'layer'):
  
  max_y = np.max(mesh.nodes[:,1]) # 2D and 3D
  
  for e in mesh.elements:
    test = mesh.calc_barycenter(e)
    if test[1] >= max_y - thickness:
      e.region = name

# set a rectangular region by min and max coordinates. For 2D and 3D
# coordinates might be out of the domain
# returns the number of changed elements
def set_rect_region(mesh, min, max, name = 'box'):
  assert len(mesh.nodes[0]) == len(min) == len(max)
  
  for e in mesh.elements:
    test = mesh.calc_barycenter(e)
    if (test >= min).all() and (test <= max).all(): # the comparson is [True, False] for each component
      e.region = name

# create 2d LBM named elements and boundary region.
# LBM pipe_bend and two_inlet_one_outlet example as used by Pingen et al. 2007
# @param mesh assumes the inclusion to be already set
def enrich_lbm2d(mesh, type):

  set_ghost_region(mesh, 'boundary')

  eps = 1e-4
  nx = mesh.nx 
  ny = mesh.ny
  nxx = nx+1 # for boundary conditions
  nyy = ny+1 

  if type == 'pipe_bend':
    mesh.ne.append(('inlet', list(range(int(0.1 * nx * ny + eps), int(0.3 * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int((ny - 1) * nx + 0.7 * nx - eps), int((ny - 1) * nx + 0.9 * nx + eps)))))
  elif type == 'two_inlet_one_outlet':
    mesh.ne.append(('inlet', list(range(int((0.25 - 1. / 16) * nx * ny + eps), int((0.25 + 1. / 16) * nx * ny + nx + eps), nx))))
    mesh.ne.append(('inlet', list(range(int((0.75 - 1. / 16) * nx * ny + eps), int((0.75 + 1. / 16) * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int(0.375 * nx * ny - 1 + eps), int(0.625 * nx * ny - 1 + eps), nx))))
  elif type == 'two_inlet_two_outlet':
    inletLength = 0.15 * ny
    mesh.ne.append(('inlet',list(range(int(0.2*nx*ny+eps), int(0.2*nx*ny+eps + inletLength*nx),nx)) ))
    mesh.ne.append(('inlet',list(range(int(0.8*nx*ny+eps-inletLength*nx), int(0.8*nx*ny+eps),nx)) ))
    mesh.ne.append(('outlet',list(range(int(0.2*nx*ny+eps + nx-1), int(0.2*nx*ny+eps + inletLength*nx+ nx-1),nx)) ))
    mesh.ne.append(('outlet',list(range(int(0.8*nx*ny+eps-inletLength*nx+ nx-1), int(0.8*nx*ny+eps+ nx-1),nx)) ))
  elif type == "pipe":
    mesh.ne.append(('inlet',list(range(nx,nx*(ny-1),nx))))
    mesh.ne.append(('outlet',list(range(2*nx-1,nx*ny-1,nx))))
  elif type == "diffuser":
    mesh.ne.append(('inlet',list(range(nx,nx*(ny-1),nx))))
    mesh.ne.append(('outlet',list(range(int(0.3*nx*ny+nx-1),int(0.7*nx*ny+nx-1),nx))))
  elif type == "low_in_high_out":
    mesh.ne.append(('inlet', list(range(int(0.1 * nx * ny + eps), int(0.3 * nx * ny + nx + eps), nx))))
    mesh.ne.append(('outlet', list(range(int(0.7 * nx * ny - 1 + eps), int(0.9 * nx * ny - 1 + eps), nx))))
  else:
    print("unkwnon lbm2d type '" + type + "'")
    sys.exit(-1)

  return mesh

# create a lvm 2d solar heater with chimney
# @param chimney height in m
def enrich_solar_heater_2d(mesh, chimney):

  eps = 1e-4
  nx = mesh.nx 
  ny = mesh.ny
  nxx = nx+1 # for boundary conditions
  nyy = ny+1 

  # rename to design to not take over legacy mech name for lbm
  mesh.rename_region('mech', 'design', 1)
  set_ghost_region(mesh, 'boundary')
  
  # ----------------------|outlet|
  # |  obstacle           |      | chimney height chimney
  # |                     |      | 
  # -----------------------      |
  # i                            |
  # n                            |
  # l   design                   |
  # e                            |
  # t                            |
  # ------------------------------
  
  inlet = []
  outlet = []
  inlet_nodes = [] # need to be made unique, but a set is not ordered
  outlet_nodes = []
  
  surf = []
  
  assert mesh.dy is not None
  height = mesh.dy * ny
  assert height > chimney + .05 # we have 10 for the chimney
  
  total_ny = ny
  wall_ny  = int(chimney / mesh.dy)
  coll_ny  = total_ny - wall_ny # collector including boundary
  print('total_ny',total_ny,'wall_ny',wall_ny,'coll_ny',coll_ny)
  
  n_outlet = .1 * nx 
  assert .8 * coll_ny <= coll_ny - 2 # ensure space for boundary
  
  for i, e in enumerate(mesh.elements):
    x, y = mesh.element_idx_by_pos(i)
    if x == 0:
      if y >= .2 * coll_ny and y <= .8 * coll_ny:
        inlet.append(i)
        assert len(e.nodes) == 4
        inlet_nodes.append(e.nodes[0]) # left lower
        inlet_nodes.append(e.nodes[3]) # left upper
      
        se = Element('inlet_surf', Ansys.LINE, 1)
        se.nodes = [e.nodes[0], e.nodes[3]]
        surf.append(se) 
    if y >= coll_ny and y != total_ny -1 and x > 0 and x < nx - n_outlet -1:
      e.region = 'obstacle'
    if y == ny-1 and x >= nx - n_outlet -2:  
      if x > nx - n_outlet -2 and x < nx-1:
        outlet.append(i)
        outlet_nodes.append(e.nodes[2])
        outlet_nodes.append(e.nodes[3])

        se = Element('outlet_surf', Ansys.LINE, 1)
        se.nodes = [e.nodes[2], e.nodes[3]]
        surf.append(se) 

  outlet_nodes.sort()
  inlet_nodes.sort()

  mesh.elements.extend(surf)
  mesh.ne.append(('inlet', inlet))
  mesh.ne.append(('outlet', outlet))
  mesh.bc.append(('inlet_nodes', list(dict.fromkeys(inlet_nodes))))
  mesh.bc.append(("outlet_nodes", list(dict.fromkeys(outlet_nodes))))
  return mesh

# create a lvm 2d solar heater with chimney
# @param chimney height in m
def enrich_solar_heater_2d_old(mesh, chimney = .1):

  eps = 1e-4
  nx = mesh.nx 
  ny = mesh.ny
  nxx = nx+1 # for boundary conditions
  nyy = ny+1 

  # rename to design to not take over legacy mech name for lbm
  mesh.rename_region('mech', 'design', 1)
  
  #                       |outlet|
  #    wall               |      | chimney height chimney
  #                       |      | 
  # -----------------------      |
  # i                            |
  # n                            |
  # l   design                   |
  # e                            |
  # t                            |
  # ------------------------------
  
  inlet = []
  outlet = []
  inlet_nodes = [] # need to be made unique, but a set is not ordered
  outlet_nodes = []
  
  surf = []
  
  assert mesh.dy is not None
  height = mesh.dy * ny
  assert height > chimney + .05 # we have 10 for the chimney
  
  total_ny = ny
  wall_ny  = int(chimney / mesh.dy)
  coll_ny  = total_ny - wall_ny # collector including boundary
  print('total_ny',total_ny,'wall_ny',wall_ny,'coll_ny',coll_ny)
  
  n_outlet = .1 * nx 
  assert .8 * coll_ny <= coll_ny - 2 # ensure space for boundary
  
  for i, e in enumerate(mesh.elements):
    x, y = mesh.element_idx_by_pos(i)
    if x == 0:
      if y <= coll_ny:
        e.region = 'boundary'
      if y >= .2 * coll_ny and y < .8 * coll_ny:
        inlet.append(i)
        assert len(e.nodes) == 4
        inlet_nodes.append(e.nodes[0]) # left lower
        inlet_nodes.append(e.nodes[3]) # left upper
      
        se = Element('inlet_surf', Ansys.LINE, 1)
        se.nodes = [e.nodes[0], e.nodes[3]]
        surf.append(se) 
    if x == nx-1:        
      e.region = 'boundary'
    if x == nx-n_outlet-2 and y >= coll_ny-1:
      e.region = 'boundary'
      
    if y == 0:
      e.region = 'boundary'
    if y == coll_ny-1 and x < nx - n_outlet -2:
      e.region = 'boundary'
    if y > coll_ny-1 and x < nx - n_outlet -2:
      e.region = 'obstacle'
    if y == ny-1 and x >= nx - n_outlet -2:  
      e.region = 'boundary'
      if x > nx - n_outlet -2 and x < nx-1:
        outlet.append(i)
        outlet_nodes.append(e.nodes[2])
        outlet_nodes.append(e.nodes[3])

        se = Element('outlet_surf', Ansys.LINE, 1)
        se.nodes = [e.nodes[2], e.nodes[3]]
        surf.append(se) 

  mesh.elements.extend(surf)

  mesh.ne.append(('inlet', inlet))
  mesh.ne.append(('outlet', outlet))
  mesh.bc.append(('inlet_nodes', list(dict.fromkeys(inlet_nodes))))
  mesh.bc.append(("outlet_nodes", list(dict.fromkeys(outlet_nodes))))
  return mesh



# see enrich_lbm2d
def enrich_lbm3d(mesh, type):

  set_ghost_region(mesh, 'boundary')

  eps = 1e-4
  
  nx = mesh.nx 
  ny = mesh.ny 
  nz = mesh.nz

  nnx = nx+1
  nny = ny+1
  nnz = nz+1

  mix, miy, miz, max, may, maz = calc_min_max_coords(mesh)
  assert mix == miy == miz == 0.0
  dx = max/nx
  dy = may/ny 
  dz = maz/nz 

  side = (("heat_bottom", []))
  mesh.bc.append(side)
  for z in range(1, nnz-1):
    for x in range(1,nnx-1):
      side[1].append((z*nny+1)*nnx+x)

  side = (("heat_top", []))
  mesh.bc.append(side)
  for z in range(1, nnz-1):
    for x in range(1,nnx-1):
      side[1].append((z*nny+ny-1)*nnx+x)

  side = (("heat_back", []))
  mesh.bc.append(side)
  for y in range(1,nny-1):
    for x in range(nnx*nnz + 1 , nnx*nnz + nnx-1,1):
      side[1].append(x+nnx*y)

  side = (("heat_front", []))
  mesh.bc.append(side)
  for y in range(1,nny-1):
    for x in range(nnx*nny*(nnz-2) + 1 , nnx*nny*(nnz-2) + nnx-1,1):
      side[1].append(x+nnx*y)

  if type == 'pipe_bend':
    area = 0.04
    in_x = int(np.sqrt(area / (dx * dx))+eps) ## find out how many elements in one direction
    dist_wall = 0.04
    x_wall = int(np.sqrt(dist_wall / (dx * dx))+eps) #find out distance from wall
    for i in range(nx*ny*nz-(x_wall+in_x)*nx*ny,nx*ny*nz-x_wall*nx*ny,nx*ny):
      mesh.ne.append(('inlet',list(range(i+x_wall*nx,i+(x_wall+in_x)*nx,nx))))
    for i in range(0,in_x*nx*ny,nx*ny):
      mesh.ne.append(('outlet', list(range((x_wall+1)*nx*ny-in_x-x_wall+i,(x_wall+1)*nx*ny-x_wall+i,1))))
  elif type == 'pipe':
    for i in range (1,nz-1):
      mesh.ne.append(('inlet', list(range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx))))
      mesh.ne.append(('outlet', list(range(int(i*nx*ny + 2*nx-1),int(i*nx*ny + 2*nx-1+nx*(ny-2)),nx))))
  elif type == 'diffuser':
    for i in range (1,nz-1):
      mesh.ne.append(('inlet', list(range(i*nx*ny+nx,(i+1)*nx*ny-nx,nx))))
    for i in range (int(0.3*nz),int(0.7*nz)):
      mesh.ne.append(('outlet', list(range(int(i*nx*ny + 0.3*nx*ny + nx-1),int(i*nx*ny + 0.7*nx*ny+ nx-1),nx))))
  elif type == 'distributor':
    center_x = nx / 2.0
    center_y = ny / 2.0
    center_z = nz / 2.0
    width_x = 0.1 * nx
    width_y = 0.1 * ny
    width_z = 0.1 * nz
    for i in range(int(round(center_z-width_z/2.0)),int(round(center_z+width_z/2.0))):
      mesh.ne.append(('outlet',list(range(int(i*nx*nz-center_x+2.0*width_x),int(i*nx*nz-center_x+3.0*width_x),1)))) # top face
      mesh.ne.append(('outlet',list(range(int((i-1)*nx*nz+center_x+2.0*width_x),int((i-1)*nx*nz+center_x+3.0*width_x),1)))) # bottom face
    for i in range(int(center_z-int(width_z)),int(center_z+int(width_z))):
      mesh.ne.append(('inlet',list(range(int(i*nx*ny+nx*int(center_y-width_y)),int(i*nx*ny+nx*int(center_y+width_y)),nx)))) #left face
    for i in range(int(round(center_y-width_y/2.0)),int(round(center_y+width_y/2.0)),1):
      mesh.ne.append(('outlet',list(range(int(i*nx+center_x+2.0*width_x),int(i*nx+center_x+3.0*width_x),1)))) #back face
      mesh.ne.append(('outlet',list(range(int(nx*ny*(nz-1)+i*nx+center_x+2.0*width_x),int(nx*ny*(nz-1)+i*nx+center_x+3.0*width_x),1)))) #front face

  return mesh

# helper wich prints info for a list of tuples (name, int) or (name, list)
def print_list(list):
  s = ''
  for n, v in list: 
    s += ' ' + n + '='
    s += str(v if type(v) is int else len(v))
  return s
  
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="This are special purpose meshes beyond create_mesh.py")
  parser.add_argument("--res", help="width resolution. Only quadratic/cubic elements", type=int, required = True )
  parser.add_argument("--width", help="in meter", type=float, default = 1.0)
  parser.add_argument("--height", help="in meter", type=float, default = 1.0)
  parser.add_argument("--depth", help="in meter", type=float, default = 1.0)

  parser.add_argument("--type", help="if not lbm, select here and configure with height and depth", choices=['bulk2d', 'bulk3d'])  
  parser.add_argument("--lbm2d", help="2D flow settings, inclusion works", choices=['pipe_bend','two_inlet_one_outlet','two_inlet_two_outlet','pipe','diffuser','low_in_high_out'])
  parser.add_argument("--lbm3d", help="3D flow settings, inclusion works", choices=['pipe','pipe_bend','diffuser','distributor'])

  parser.add_argument("--solar2d", help="2D solar heater", action='store_true')
  parser.add_argument("--solar_wall", help="wall size im m for solar2d", type=float, default = 0.2)
  
  parser.add_argument("--ghost", help="set a one element ghost region. 2D/3D", action='store_true')
  
  parser.add_argument('--inclusion', help="centered/top inclusion for 2D/3D, also for lbm", choices=["rect", "circ", "top_panel"])
  parser.add_argument('--inclusion_size', help="inclusion size in m, mandatory", type=float)
  parser.add_argument('--inclusion_name', help="if not given, defaults are used depending on type")
  
  parser.add_argument('--box', nargs='+', help="add a rectanular region 'box' by 4/6 corner coordingates in m.", type=float)
  
  parser.add_argument('--file', help="optional give output file name. ")
  
  args = parser.parse_args()
  mesh = None
  
  name = None
  
  if args.lbm2d and args.lbm3d:
    print('Error: decide which lbm to use :)')
    sys.exit(-1)
    
  if args.type and (args.lbm2d or args.lbm3d):
    print('Error: either type or lbm :)')
    sys.exit(-1)
      
  if args.type == 'bulk2d' or args.lbm2d:
    mesh = create_2d_mesh(args.res, width = args.width, height = args.height)
    if args.lbm2d:
      enrich_lbm2d(mesh, args.lbm2d)
      name = args.lbm2d
    else:
      name = args.type    
    # inclusion comes later    

  if args.type == 'bulk3d' or args.lbm3d:
    mesh = create_3d_mesh(args.res, width = args.width, height = args.height, depth=args.depth)
    if args.lbm3d:
      enrich_lbm3d(mesh, args.lbm3d)  
      name = args.lbm3d
    else:
      name = args.type    

  if args.solar2d:
    mesh = create_2d_mesh(args.res, width = args.width, height = args.height)
    enrich_solar_heater_2d(mesh, args.solar_wall)
    name = 'solar-wall_' + str(args.solar_wall) + '-nx' # chimney

  if mesh == None:
    print('Error: give either type or a lbm')
    sys.exit(-1)

  name += '_' + str(args.res)
  if args.width != 1.0:
    name += '-w_' + str(args.width)
  if args.height != 1.0:
    name += '-h_' + str(args.height)
  if args.depth != 1.0:
    name += '-d_' + str(args.depth)

  # the ghost layer is already done with name 'boundary' for lbm
  if args.ghost:
    set_ghost_region(mesh, 'ghost')
    name += '_ghost'

  # inclusions go for 2D and 3D
  if args.inclusion:
    if not args.inclusion_size:
      print('Error: inclusion_size is mandatory for inclusion')
      sys.exit(-1)
    inc_name = args.inclusion_name if args.inclusion_name  else 'inclusion' if args.type else 'obstacle'
    if args.inclusion in ['rect', 'circ']:
      set_centered_inclusion(mesh, args.inclusion, args.inclusion_size, inc_name)
    if args.inclusion == 'top_panel':
      set_top_layer(mesh, args.inclusion_size, inc_name)
    name += '_' + args.inclusion + '_' + str(args.inclusion_size)    
      
  if args.box:
    if mesh.is2d():
      if len(args.box) != 4:
        print('Error: expect for --box x_min y_min x_max y_max')
        set_rect_region(mesh, min=[args.box[0],args.box[1]],max=[args.box[2],args.box[3]],name='box')
    else:
      if len(args.box) != 6:
        print('Error: expect for --box x_min y_min z_min x_max y_max z_max')
        set_rect_region(mesh, min=[args.box[0],args.box[1],args.box[2]],max=[args.box[3],args.box[4],args.box[5]],name='box')
    name += '_box'          
      
  file = name + '.mesh' if args.file == None else args.file 
  
  write_ansys_mesh(mesh, file)
  
  print('regions:' + print_list(mesh.count_regions()))
  if len(mesh.ne) > 0:
    print('named elements:' + print_list(mesh.ne))
  print('named nodes:' + print_list(mesh.bc))
  print("created file '" + file + "' with " + str(len(mesh.elements)) + " elements")
