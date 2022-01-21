import vtk
import pickle
#from vtk.util.colors import *
from numpy import *
from matviz_rot import *
from matviz_vtk import *
# polyhedron lib from: http://cens.ioc.ee/projects/polyhedron/
#from polyhedron import Vrep, Hrep
import scipy.interpolate as ip



def create_3d_mesh_unstructured(coords, nondes_coords, nodes_force, nodes_support, design, ip_nx, ip_ny, ip_nz, grad, scale):

  st1 = design['s1']
  st2 = design['s2']
  st3 = design['s3']
  angle = design['angle']

  mesh = Mesh()
  centers, min, max = coords[0:3]  # design nodes
  nondes_centers, nondes_min, nondes_max = nondes_coords[0:3]  # nondesign nodes
  if scale <= 0:
    scale = 1.0
  # calculate size of an element dx,dy,dz
  delta = (abs(max[0] - min[0]), abs(max[1] - min[1]), abs(max[2] - min[2]))
  dx = float(delta[0] / ip_nx)
  dy = float(delta[1] / ip_ny)
  dz = float(delta[2] / ip_nz)
  # number of elements nx,ny,nz in each direction
  nx = ip_nx
  ny = ip_ny
  nz = ip_nz

  # non-design region points are added to design region points in order to create a structured mesh
  centers_full = pickle.loads(pickle.dumps(centers)) # pickle is faster than deepcopy
  for i in range(len(nondes_centers)):
    centers_full.append((nondes_centers[i]))
  coords_full = (centers_full, (numpy.min((min[0], nondes_min[0])), numpy.min((min[1], nondes_min[1])), numpy.min((min[2], nondes_min[2]))), (numpy.max((max[0], nondes_max[0])), numpy.max((max[1], nondes_max[1])), numpy.max((max[2], nondes_max[2]))))
  centers_full, min_full, max_full = coords_full[0:3]

  delta_full = (abs(max_full[0] - min_full[0]), abs(max_full[1] - min_full[1]), abs(max_full[2] - min_full[2]))
  nx_full = int(delta_full[0] / dx)
  ny_full = int(delta_full[1] / dy)
  nz_full = int(delta_full[2] / dz)

  # calculate the upper and lower bound of the discretization with respect to min,max of design nodes
  overhead_x_l = int(round(abs(min_full[0] - min[0]) / dx))
  overhead_y_l = int(round(abs(min_full[1] - min[1]) / dy))
  overhead_z_l = int(round(abs(min_full[2] - min[2]) / dz))

  overhead_x_u = int(round(abs(max_full[0] - max[0]) / dx))
  overhead_y_u = int(round(abs(max_full[1] - max[1]) / dy))
  overhead_z_u = int(round(abs(max_full[2] - max[2]) / dz))
  # add regular grid nodes to mesh
  for z in range(-overhead_z_l, nz + 1 + overhead_z_u):
    for y in range(-overhead_y_l, ny + 1 + overhead_y_u):
      for x in range(-overhead_x_l, nx + 1 + overhead_x_u):
        mesh.nodes.append([min[0] + float(x) * dx, min[1] + float(y) * dy, min[2] + float(z) * dz])
  
  mesh.nodes = np.array(mesh.nodes) # make the expected np.array
        
  # initialize pseudo-designvariables for design and nondesign region
  st1_full = zeros((len(centers_full) , 1))
  st2_full = zeros((len(centers_full), 1))
  st3_full = zeros((len(centers_full), 1))
  idx = 0
  for i in range(len(st1_full)):
    if idx < len(st1):
      st1_full[i, 0] = st1[idx, 0]
      st2_full[i, 0] = st2[idx, 0]
      st3_full[i, 0] = st3[idx, 0]
      idx += 1
    else:
      st1_full[i, 0] = 1.
      st2_full[i, 0] = 1.
      st3_full[i, 0] = 1.

  # add pseudo angles for non-design region
  if angles != None:
    ang_full = zeros((len(centers_full), 3))
    for i in range((ang_full.shape)[0]):
      if i < (angles.shape)[0]:
        ang_full[i, :] = angles[i, :]
      else:
        ang_full[i, :] = zeros((1, 3))
  # interpolate design values with nearest neighbor interpolation for design region
  ip_data, ip_near, out, ndim = get_interpolation_unstructured(coords, grad, st1, st2, st3, 0, 0, 0, nx, ny, nz, dx, dy, dz, None, None, angles)
  # interpolate design values by linear interpolation for non-design and design region
  ip_data_full, ip_near_full, out_full, ndim_full = get_interpolation_unstructured(coords_full, grad, st1_full, st2_full, st3_full, -overhead_x_l, -overhead_y_l, - -overhead_z_l, nx + overhead_x_u, ny + overhead_y_u, nz + overhead_z_u, dx, dy, dz, min, max, ang_full)

  i = 0
  idx = 0
  within = 0
  invalid = 0
  nnx = nx + overhead_x_u + overhead_x_l
  nny = ny + overhead_y_u + overhead_y_l
  nnz = nz + overhead_z_u + overhead_z_l
  for z in range(-overhead_z_l, nz + overhead_z_u):
    for y in range(-overhead_y_l, ny + overhead_y_u):
      for x in range(-overhead_x_l, nx + overhead_x_u):
        coord = out_full[idx]
        if x >= 0 and x < nx and y >= 0 and y < ny and z >= 0 and z < nz:
          s1, s2, s3 = ip_data[i, 0:3]
          angle = None if angles is None else ip_data[i, 3:6]
          i += 1
        else:
          s1, s2, s3 = [-1., -1., -1.]
          angle = zeros((3, 1))

        s1_full, s2_full, s3_full = ip_data_full[idx, 0:3]
        angle_full = None if angles is None else ip_data_full[idx, 3:6]
        e = Element()
        # if s1 < 0 point is out of the convex hull
        # print 's1_schleife = ' + repr(s1)
        if s1 > 0.0:
          within += 1
          # remove holes of the robot arm
          if not valid_position(coord, coords):
            invalid += 1
            e.region = 'void'
          else:
            e.region = 'design'
        else:
          if s1_full > 0.:
            # remove holes of the robot arm
            if not valid_position(coord, coords):
              invalid += 1
              e.region = 'void'
            else:
              e.region = 'nondesign'
          else:
            e.region = 'void'
        e.type = HEXA8
        e.stiff1 = s1
        e.stiff2 = s2
        e.stiff3 = s3
        e.rotX = None if angles is None else angle[0]
        e.rotY = None if angles is None else angle[1]
        e.rotZ = None if angles is None else angle[2]
        ll = (nnx + 1) * (nny + 1) * (z + overhead_z_l) + (nnx + 1) * (y + overhead_y_l) + (x + overhead_x_l)  # lowerleft
        e.nodes = ((ll + (nnx + 1) * (nny + 1), ll + (nnx + 1) * (nny + 1) + nnx + 1, ll + (nnx + 1) * (nny + 1) + nnx + 1 + 1, ll + (nnx + 1) * (nny + 1) + 1, ll, ll + nnx + 1, ll + nnx + 1 + 1, ll + 1))
        mesh.elements.append(e)
        idx += 1
  # add support and load nodes to the mesh
  min_mod = (min[0] - overhead_x_l * dx, min[1] - overhead_y_l * dy, min[2] - overhead_z_l * dz)
  max_mod = (min[0] + (overhead_x_u + nx) * dx, min[1] + (ny + overhead_y_u) * dy, min[2] + (nz + overhead_z_u) * dz)
  mesh = map_bc_to_node(mesh, nodes_force, min_mod, max_mod, 'load', dx, dy, dz, nnx + 1, nny + 1, nnz + 1, -1)
  mesh = map_bc_to_node(mesh, nodes_support, min_mod, max_mod, 'support', dx, dy, dz, nnx + 1, nny + 1, nnz + 1, 1)

  if grad != 'nearest':
    print(str(within) + ' elements out of ' + str(len(out)) + ' in convex hull')
  if invalid > 0:
    print(str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)')
  return mesh







def map_bc_to_node(mesh, nodes, min, max, name, dx, dy, dz, nnx, nny, nnz, flag):
  # map boundary nodes from unstructured grid to structured grid
  elem = []
  for i in range(len(nodes)):
    value = nodes[i]
    ind_x = abs(value[0] - min[0]) / dx
    ind_y = abs(value[1] - min[1]) / dy
    ind_z = abs(value[2] - min[2]) / dz
    index = int(nny * nnx * int(ind_z) + nnx * int(ind_y) + int(ind_x))
    out_ind = 0.
    if ind_x >= round(ind_x):
      if ind_y >= round(ind_y):
        if ind_z >= round(ind_z):
          out_ind = index
        else:
          out_ind = index + (nnx) * (nny)
      else:
        if ind_z >= round(ind_z):
          out_ind = index + nnx
        else:
          out_ind = index + (nnx) * (nny) + nnx
    else:
      if ind_y >= round(ind_y):
        if ind_z >= round(ind_z):
          out_ind = index + 1
        else:
          out_ind = index + (nnx) * (nny) + 1
      else:
        if ind_z >= round(ind_z):
          out_ind = index + nnx + 1
        else:
          out_ind = index + (nnx) * (nny) + nnx + 1
    if flag > 0:
      elem.append(out_ind + 2 * nnx)
    else:
      elem.append(out_ind - 2 * nnx)
  mesh.bc.append((name, elem))
  return mesh


def get_interpolation_unstructured(coords, grad, s1, s2, s3, lx, ly, lz, ux, uy, uz, dx, dy, dz, mi=None, ma=None, angle=None):
  # we make our own elem
  if mi != None and max != None:
    centers = coords[0]  # skip elem
  else:
    centers, mi, ma = coords[0:3]
  if ux - lx == 0 or uy - ly == 0 or uz - lz == 0:
    print('chose a higher hom_samples such that also the smallest side gets discretized')
    exit()

  out = numpy.zeros(((ux - lx) * (uy - ly) * (uz - lz), 3))
  idx = 0
  for z in range(lz, uz):
    for y in range(ly, uy):
      for x in range(lx, ux):
        out[idx] = ((mi[0] + 0.5 * dx + float(x) * dx, mi[1] + 0.5 * dy + float(y) * dy, mi[2] + 0.5 * dz + float(z) * dz))
        idx += 1

  v = numpy.zeros((len(s1), 3 if angle == None else 6))

  v[:, 0] = s1[:, 0]
  v[:, 1] = s2[:, 0]
  v[:, 2] = s3[:, 0]
  v[:, 3:6] = angle[:, :]
  print('interpolation start')
  ip_data = ip.griddata(centers, v, out, grad, -500.)
  # any interpolation, ie. linear interpolation can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation which can also interpolate values outside the convex hull
  ip_near = ip.griddata(centers, v, out, 'nearest') if grad != 'nearest' else None

  print('interpolation end')

  return ip_data, ip_near, out, (ux - lx, uy - ly, uz - lz)


# # for the robot arm we have check for the two nondesign holes as they are within the
# # convex hull of the design :(
def valid_position(pos, coords):
  mi, ma = coords[1:3]
  delta = (abs(ma[0] - mi[0]), abs(ma[1] - mi[1]), abs(ma[2] - mi[2]))
  # if int(delta[0]) == 508 and int(delta[2]) == 126:
  if (pos[0] + 147.4) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center -147, 0, 0
    return False
  if (pos[0] - 250.0) ** 2 + pos[2] ** 2 < 30.0 ** 2:  # center 250, 0, 0
    return False


  return True


#def mkhull(points):
#    p = Vrep (points)
#    return Hrep (p.A, p.b)

def point_inside_polygon(point, poly):
  #     -----------------------------------------------
  # C-Library cddlib (version 0.91) README FILE
  # -----------------------------------------------
  # 1. The C-library  cddlib is a C implementation of the Double Description
  # Method of Motzkin et al. for generating all vertices (i.e. extreme points)
  # and extreme rays of a general convex polyhedron in R^d given by a system
  # of linear inequalities:
  #
  #    P = { x=(x1, ..., xd)^T :  b - A  x  >= 0 }
  #
  # where  A  is a given m x d real matrix, b is a given m-vector
  # and 0 is the m-vector of all zeros.
  #     return inside
  if all(dot(p.A, point) <= p.b):
    return True
  else:
    return False

def write_node_design(name, mesh, coords, design):

  st1 = design['s1']
  st2 = design['s2']
  st3 = design['s3']
  angles = design['angle']

  # maps element density to element corner nodes
  nodes_st1 = zeros((len(mesh.nodes), 1))
  nodes_st2 = zeros((len(mesh.nodes), 1))
  nodes_st3 = zeros((len(mesh.nodes), 1))
  if angles != None:
    nodes_angles = zeros((len(mesh.nodes), 3))
  else:
    nodes_angles = None
  nb_des_count = zeros((len(mesh.nodes), 1))
  nb_ndes_count = zeros((len(mesh.nodes), 1))

  count = 0
  for i in range(len(mesh.elements)):
    if mesh.elements[i].region == 'design':
      for j in range(len(mesh.elements[i].nodes)):
        nodes_st1[mesh.elements[i].nodes[j]] += st1[count]
        nodes_st2[mesh.elements[i].nodes[j]] += st2[count]
        nodes_st3[mesh.elements[i].nodes[j]] += st3[count]
        if nodes_angles != None:
          nodes_angles[mesh.elements[i].nodes[j], :] += angles[count, :]
        nb_des_count[mesh.elements[i].nodes[j]] += 1
      count += 1
    elif mesh.elements[i].region == 'nondesign':
      for j in range(len(mesh.elements[i].nodes)):
        nb_ndes_count[mesh.elements[i].nodes[j]] += 1
  output_file = open(name, "w")
  output_file.write(str(len(nodes_st1)) + '  ' + str(6 if angles != None else 3) + ' \n')
  for i in range(len(mesh.nodes)):
    if nb_des_count[i] != 0:
      nodes_st1[i] /= nb_des_count[i]
      nodes_st2[i] /= nb_des_count[i]
      nodes_st3[i] /= nb_des_count[i]
      nodes_angles[i][:] /= nb_des_count[i]
    elif nb_ndes_count[i] == 0 and nb_des_count[i] == 0:
      # void node
      nodes_st1[i] = -10
      nodes_st2[i] = -10
      nodes_st3[i] = -10
    elif nb_des_count[i] == 0:
      # nondesign node
      nodes_st1[i] = -1
      nodes_st2[i] = -1
      nodes_st3[i] = -1
    if angles != None:
      output_file.write(str(nodes_st1[i]) + '\t ' + str(nodes_st2[i]) + '\t ' + str(nodes_st3[i]) + '\t ' + str(nodes_angles[i, :]) + ' \n')
    else:
      output_file.write(str(nodes_st1[i]) + '\t ' + str(nodes_st2[i]) + '\t ' + str(nodes_st3[i]) + ' \n')
  output_file.close()


def post_process_mesh(mesh, coords, nondes_coords, ip_nx, ip_ny, ip_nz):
  centers, min, max = coords[0:3]  # design nodes
  nondes_centers, nondes_min, nondes_max = nondes_coords[0:3]  # design nodes
  # calculate size of an element dx,dy,dz
  delta = (abs(max[0] - min[0]), abs(max[1] - min[1]), abs(max[2] - min[2]))
  dx = float(delta[0] / ip_nx)
  dy = float(delta[1] / ip_ny)
  dz = float(delta[2] / ip_nz)
#   print 'dx = ' + repr(dx)
  # number of elements nx,ny,nz in each direction
  nx = ip_nx
  ny = ip_ny
  nz = ip_nz

  # non-design region points are added to design region points in order to create a structured mesh
  centers_full = pickle.loads(pickle.dumps(centers)) # pickle is faster than deepcopy
  for i in range(len(nondes_centers)):
    centers_full.append((nondes_centers[i]))
  coords_full = (centers_full, (numpy.min((min[0], nondes_min[0])), numpy.min((min[1], nondes_min[1])), numpy.min((min[2], nondes_min[2]))), (numpy.max((max[0], nondes_max[0])), numpy.max((max[1], nondes_max[1])), numpy.max((max[2], nondes_max[2]))))
  centers_full, min_full, max_full = coords_full[0:3]

  delta_full = (abs(max_full[0] - min_full[0]), abs(max_full[1] - min_full[1]), abs(max_full[2] - min_full[2]))
  nx_full = int(delta_full[0] / dx)
  ny_full = int(delta_full[1] / dy)
  nz_full = int(delta_full[2] / dz)

  overhead_x_l = int(round(abs(min_full[0] - min[0]) / dx))
  overhead_y_l = int(round(abs(min_full[1] - min[1]) / dy))
  overhead_z_l = int(round(abs(min_full[2] - min[2]) / dz))

  overhead_x_u = int(round(abs(max_full[0] - max[0]) / dx))
  overhead_y_u = int(round(abs(max_full[1] - max[1]) / dy))
  overhead_z_u = int(round(abs(max_full[2] - max[2]) / dz))

  nnx = nx + overhead_x_u + overhead_x_l
  nny = ny + overhead_y_u + overhead_y_l
  nnz = nz + overhead_z_u + overhead_z_l
  idx = 0
  i = 0
  print('mesh.elements = ' + repr(len(mesh.elements)))
  print('nnx = ' + repr(nnx))
  print('nny = ' + repr(nny))
  print('nnz = ' + repr(nnz))
  yrange = list(range(-overhead_x_l, nx + overhead_y_u))
  for z in range(nnz):
    for x in range(nnx):
      idx = z * (nnx) * (nny) + x
      idx2 = idx + nnx
      idx3 = idx + 2 * nnx
      idx4 = idx + 3 * (nnx)
      tmp = mesh.elements[idx2]
      tmp2 = mesh.elements[idx3]
      tmp3 = mesh.elements[idx4]
#       print 'tmp = ' + repr(tmp.region)
#       print 'tm2 = ' + repr(tmp2.region)
      if tmp.region == 'nondesign' and tmp2.region == 'nondesign' and tmp3.region == 'nondesign':
#         print 'drin'
        i = 4
        while i < nny - 1:
          idx4 = idx + i * (nnx)
          tmp4 = mesh.elements[idx4]
          tmp4.region = 'nondesign'
          i += 1
#   print 'idx ende ' + repr(idx)
  return mesh

# def create_3d_mesh_unstructured_5(coords, nondes_coords, nodes_force, nodes_support, st1, st2, st3, angles, ip_nx, ip_ny, ip_nz, grad, scale):
#   mesh = Mesh()
#   centers, min, max = coords[0:3]  # design nodes
# #   print 'coords_len = ' + repr(len(coords[0]))
# #   print 'min = ' + repr(min)
# #   print 'max = ' + repr(max)
#   nondes_centers, nondes_min, nondes_max = nondes_coords[0:3]  # nondesign nodes
#   if scale <= 0:
#     scale = 1.0
#   # calculate size of an element dx,dy,dz
#   delta = (abs(max[0] - min[0]), abs(max[1] - min[1]), abs(max[2] - min[2]))
#   dx = float(delta[0] / ip_nx)
#   dy = float(delta[1] / ip_ny)
#   dz = float(delta[2] / ip_nz)
# #   print 'dx = ' + repr(dx)
#   # number of elements nx,ny,nz in each direction
#   nx = ip_nx
#   ny = ip_ny
#   nz = ip_nz
#
#   # non-design region points are added to design region points in order to create a structured mesh
#   centers_full = pickle.loads(pickle.dumps(centers)) # pickle is faster than deepcopy
#   for i in range(len(nondes_centers)):
#     centers_full.append((nondes_centers[i]))
#   coords_full = (centers_full, (numpy.min((min[0], nondes_min[0])), numpy.min((min[1], nondes_min[1])), numpy.min((min[2], nondes_min[2]))), (numpy.max((max[0], nondes_max[0])), numpy.max((max[1], nondes_max[1])), numpy.max((max[2], nondes_max[2]))))
#   centers_full, min_full, max_full = coords_full[0:3]
#   print 'min = ' + repr(min)
#   print 'max = ' + repr(max)
#   print 'min_full = ' + repr(min_full)
#   print 'max_full = ' + repr(max_full)
#
#
#   delta_full = (abs(max_full[0] - min_full[0]), abs(max_full[1] - min_full[1]), abs(max_full[2] - min_full[2]))
#   nx_full = int(delta_full[0] / dx)
#   ny_full = int(delta_full[1] / dy)
#   nz_full = int(delta_full[2] / dz)
#
# #   overhead_x = round((nx_full - nx) * 0.5)
# #   overhead_y = round((ny_full - ny) * 0.5)
# #   overhead_z = round((nz_full - nz) * 0.5)
#
#   overhead_x_l = int(round(abs(min_full[0] - min[0]) / dx))
#   overhead_y_l = int(round(abs(min_full[1] - min[1]) / dy))
#   overhead_z_l = int(round(abs(min_full[2] - min[2]) / dz))
#
#   overhead_x_u = int(round(abs(max_full[0] - max[0]) / dx))
#   overhead_y_u = int(round(abs(max_full[1] - max[1]) / dy))
#   overhead_z_u = int(round(abs(max_full[2] - max[2]) / dz))
#
# #   overhead_x_l = 0
# #   overhead_y_l = 0
# #   overhead_z_l = 0
# #
# #   overhead_x_u = 0
# #   overhead_y_u = 0
# #   overhead_z_u = 0
# #   nx2 = nx + 1 + 2 * overhead_x
# #   ny2 = ny + 1 + 2 * overhead_y
# #   nz2 = nz + 1 + 2 * overhead_z
# #
# #   num = (nx2 + 1) * (ny2 + 1) * (nz2 + 1)
# #   st1_full = zeros((num, 1))
# #   st2_full = zeros((num, 1))
# #   st3_full = zeros((num, 1))
# #   ang_full = zeros((num, 3))
# #   idx = 0
# #   i = 0
# #   for z in range(-overhead_z, nz + 1 + overhead_z):
# #     for y in range(-overhead_y, ny + 1 + overhead_y):
# #       for x in range(-overhead_x, nx + 1 + overhead_x):
# #         mesh.nodes.append((min[0] + float(x) * dx, min[1] + float(y) * dy, min[2] + float(z) * dz))
#   idx = 0
#   for z in range(-overhead_z_l, nz + 1 + overhead_z_u):
#     for y in range(-overhead_y_l, ny + 1 + overhead_y_u):
#       for x in range(-overhead_x_l, nx + 1 + overhead_x_u):
#         mesh.nodes.append((min[0] + float(x) * dx, min[1] + float(y) * dy, min[2] + float(z) * dz))
#         idx += 1
# #         if x >= 0 and x <= nx and y >= 0 and y <= ny and z >= 0 and z <= nz:
# #           st1_full[idx, 0] = st1[i]
# #           st2_full[idx, 0] = st2[i]
# #           st3_full[idx, 0] = st3[i]
# #           if angles <> None:
# #             ang_full[idx, 3:6] = angles[i, 3:6]
# #           i += 1
# #         idx += 1
#   st1_full = zeros((len(centers_full) , 1))
#   st2_full = zeros((len(centers_full), 1))
#   st3_full = zeros((len(centers_full), 1))
#   idx = 0
#   for i in range(len(st1_full)):
#     if idx < len(st1):
#       st1_full[i, 0] = 1.
#       st2_full[i, 0] = 1.
#       st3_full[i, 0] = 1.
#       idx += 1
#     else:
#       st1_full[i, 0] = -1.
#       st2_full[i, 0] = -1.
#       st3_full[i, 0] = -1.
#
# #   print 'shape angle = ' + repr(len(angles))
# #   print 'angles = ' + repr(angles[0])
#
#   # add pseudo angles for non-design region
#   if angles <> None:
# #     print 'angles = ' + repr(angles.shape)
#     ang_full = zeros((len(centers_full), 3))
#     for i in range((ang_full.shape)[0]):
#       if i < (angles.shape)[0]:
#         ang_full[i, :] = angles[i, :]
#       else:
#         ang_full[i, :] = zeros((1, 3))
#
#   print 'st1.shape = ' + repr(st1.shape)
#   print 'coord.shape = ' + repr(len(coords[0]))
#   print 'nx =' + repr(nx)
#   print 'ny =' + repr(ny)
#   print 'nz =' + repr(nz)
#   print 'nx_full =' + repr(nx_full)
#   print 'ny_full =' + repr(ny_full)
#   print 'nz_full =' + repr(nz_full)
#   print 'overhead_x_l =' + repr(overhead_x_l)
#   print 'overhead_y_l =' + repr(overhead_y_l)
#   print 'overhead_z_l =' + repr(overhead_z_l)
#
#   print 'overhead_x_u =' + repr(overhead_x_u)
#   print 'overhead_y_u =' + repr(overhead_y_u)
#   print 'overhead_z_u =' + repr(overhead_z_u)
#   # ip_data, ip_near, out, ndim = get_interpolation_unstructured(coords, grad, st1, st2, st3, 0, 0, 0, nx, ny, nz, dx, dy, dz, None, None, angles)
#   print 'st1_full.shape = ' + repr(st1_full.shape)
#   print 'coord_full.shape = ' + repr(len(coords_full[0]))
#   print 'type nondes_centers = ' + repr(type(nondes_centers))
# #   non = []
# #   min = [10000., 10000., 10000]
# #   max = [-10000., -10000., -10000.]
# #   for i in range(len(nondes_centers)):
# #     if nondes_centers[i][1] <= -27.:
# #       for j in range(3):
# #         if min[j] > nondes_centers[i][j]:
# #           min[j] = nondes_centers[i][j]
# #         if max[j] < nondes_centers[i][j]:
# #           max[j] = nondes_centers[i][j]
# #       non.append(nondes_centers[i])
# #
# #   nonn = (non, min, max)
# #   print 'test min' + repr(min)
# #   print 'test max' + repr(max)
# #   st1 = ones((len(non), 1))
# #   st2 = ones((len(non), 1))
# #   st3 = ones((len(non), 1))
# #   angles = zeros((len(non), 3))
#   print 'len of centers = ' + repr(len(centers))
#   print 'len of centers_full = ' + repr(len(centers_full))
#   for i in range(len(centers_full)):
#     if centers_full[i][0] >= -118. and centers_full[i][0] < -112. and centers_full[i][1] >= -17. and centers_full[i][1] < -11. and centers_full[i][2] >= 30. and centers_full[i][2] < 38.:
#       print 'Nummer = ' + repr(i)
#       print 'center_full = ' + repr(centers_full[i])
#       print 'stiff1 = ' + repr(st1_full[i, 0])
#   print 'center_nondes[75] = ' + repr(nondes_centers[75])
#   print 'center_nondes[76] = ' + repr(nondes_centers[76])
#   ip_data_full, ip_near_full, out_full, ndim_full = get_interpolation_unstructured(coords_full, 'nearest', st1_full, st2_full, st3_full, -overhead_x_l, -overhead_y_l, -overhead_z_l, nx + overhead_x_u, ny + overhead_y_u, nz + overhead_z_u, dx, dy, dz, min, max, ang_full)
#   ip_data_full2, ip_near_full2, out_full2, ndim_full2 = get_interpolation_unstructured(coords_full, grad, st1_full, st2_full, st3_full, -overhead_x_l, -overhead_y_l, -overhead_z_l, nx + overhead_x_u, ny + overhead_y_u, nz + overhead_z_u, dx, dy, dz, min, max, ang_full)
#   ip_data_full3, ip_near_full3, out_full3, ndim_full3 = get_interpolation_unstructured(coords, grad, st1, st2, st3, overhead_x_l, -overhead_y_l, -overhead_z_l, nx + overhead_x_u, ny + overhead_y_u, nz + overhead_z_u, dx, dy, dz, min, max, angles)
#
#
#   i = 0
#   idx = 0
#   within = 0
#   invalid = 0
#   nnx = nx + overhead_x_u + overhead_x_l
#   nny = ny + overhead_y_u + overhead_y_l
#   nnz = nz + overhead_z_u + overhead_z_l
#   for z in range(-overhead_z_l, nz + overhead_z_u):
#     for y in range(-overhead_y_l, ny + overhead_y_u):
#       for x in range(-overhead_x_l, nx + overhead_x_u):
#         coord = out_full[idx]
# #         s1_full, s2_full, s3_full = ip_data_full[idx, 0:3]
# #         s1_full2, s2_full2, s3_full2 = ip_data_full2[idx, 0:3]
# #         s1_full3, s2_full3, s3_full3 = ip_data_full3[idx, 0:3]
# #
# #         # assert(s1_full == -1. or s1_full == 1.)
# #         angle_full = None if angles is None else ip_data_full[idx, 3:6]
#         s1_full = ip_data_full[idx, 0]
#         s1_full2 = ip_data_full2[idx, 0]
#         s1_full3 = ip_data_full3[idx, 0]
#         if idx == 13847:
#           print 's1_full_nearest = ' + repr(s1_full)
#           print 's1_full_linear = ' + repr(s1_full2)
#           print 's1_linear = ' + repr(s1_full3)
#           print 'coord = ' + repr(coord)
#         e = Element()
#         # if s1 < 0 point is out of the convex hull
#         # print 's1_schleife = ' + repr(s1)
#         if s1_full2 > 0.:
#           within += 1
#           if not valid_position(coord, coords):
#             invalid += 1
#             e.region = 'void'
#           else:
#             if s1_full == 1.:
#               if s1_full3 > 0.:
#                 e.region = 'design'
#               else:
#                 e.region = 'nondesign'
#             else:
#               e.region = 'nondesign'
#         else:
#           e.region = 'void'
#         e.type = HEXA8
#         e.stiff1 = 1.
#         e.stiff2 = 1.
#         e.stiff3 = 1.
# #         e.rotX = None if angles is None else 0.
# #         e.rotY = None if angles is None else 0.
# #         e.rotZ = None if angles is None else 0.
#         e.rotX = 0.
#         e.rotY = 0.
#         e.rotZ = 0.
#         ll = (nnx + 1) * (nny + 1) * (z + overhead_z_l) + (nnx + 1) * (y + overhead_y_l) + (x + overhead_x_l)  # lowerleft
#         e.nodes = ((ll + (nnx + 1) * (nny + 1), ll + (nnx + 1) * (nny + 1) + nnx + 1, ll + (nnx + 1) * (nny + 1) + nnx + 1 + 1, ll + (nnx + 1) * (nny + 1) + 1, ll, ll + nnx + 1, ll + nnx + 1 + 1, ll + 1))
#         mesh.elements.append(e)
#         idx += 1
#   if grad <> 'nearest':
#     print str(within) + ' elements out of ' + str(len(out_full)) + ' in convex hull'
#   if invalid > 0:
#     print str(invalid) + ' elements out of ' + str(len(out_full)) + ' are considered invalid (shall not be visualized)'
#   return mesh
# def create_3d_mesh_unstructured_4(coords, nondes_coords, nodes_force, nodes_support, st1, st2, st3, angles, ip_nx, ip_ny, ip_nz, grad, scale):
#   mesh = Mesh()
#   centers, min, max = coords[0:3]  # design nodes
#   nondes_centers, nondes_min, nondes_max = nondes_coords[0:3]  # nondesign nodes
#   if scale <= 0:
#     scale = 1.0
#   # calculate size of an element dx,dy,dz
#   delta = (abs(max[0] - min[0]), abs(max[1] - min[1]), abs(max[2] - min[2]))
#   dx = float(delta[0] / ip_nx)
#   dy = float(delta[1] / ip_ny)
#   dz = float(delta[2] / ip_nz)
#   # number of elements nx,ny,nz in each direction
#   nx = ip_nx
#   ny = ip_ny
#   nz = ip_nz
#     # non-design region points are added to design region points in order to create a structured mesh
#   centers_full = pickle.loads(pickle.dumps(centers)) # pickle is faster than deepcopy
#   for i in range(len(nondes_centers)):
#     centers_full.append((nondes_centers[i]))
#   coords_full = (centers_full, (numpy.min((min[0], nondes_min[0])), numpy.min((min[1], nondes_min[1])), numpy.min((min[2], nondes_min[2]))), (numpy.max((max[0], nondes_max[0])), numpy.max((max[1], nondes_max[1])), numpy.max((max[2], nondes_max[2]))))
#   centers_full, min_full, max_full = coords_full[0:3]
#   print 'min = ' + repr(min)
#   print 'max = ' + repr(max)
#   print 'min_full = ' + repr(min_full)
#   print 'max_full = ' + repr(max_full)
#
#
#   delta_full = (abs(max_full[0] - min_full[0]), abs(max_full[1] - min_full[1]), abs(max_full[2] - min_full[2]))
#   nx_full = int(delta_full[0] / dx)
#   ny_full = int(delta_full[1] / dy)
#   nz_full = int(delta_full[2] / dz)
#
#   overhead_x_l = int(round(abs(min_full[0] - min[0]) / dx))
#   overhead_y_l = int(round(abs(min_full[1] - min[1]) / dy))
#   overhead_z_l = int(round(abs(min_full[2] - min[2]) / dz))
#
#   overhead_x_u = int(round(abs(max_full[0] - max[0]) / dx))
#   overhead_y_u = int(round(abs(max_full[1] - max[1]) / dy))
#   overhead_z_u = int(round(abs(max_full[2] - max[2]) / dz))
#
#   nnx = nx + overhead_x_u + overhead_x_l
#   nny = ny + overhead_y_u + overhead_y_l
#   nnz = nz + overhead_z_u + overhead_z_l
#
#
#   out = numpy.zeros(((nnx + 1) * (nny + 1) * (nnz + 1), 3))
#   for z in range(-overhead_z_l, nz + 1 + overhead_z_u):
#     for y in range(-overhead_y_l, ny + 1 + overhead_y_u):
#       for x in range(-overhead_x_l, nx + 1 + overhead_x_u):
#         mesh.nodes.append((min[0] + float(x) * dx, min[1] + float(y) * dy, min[2] + float(z) * dz))
#
#   idx = 0
#   for z in range(-overhead_z_l, nz + overhead_z_u):
#     for y in range(-overhead_y_l, ny + overhead_y_u):
#       for x in range(-overhead_x_l, nx + overhead_x_u):
#         out[idx] = ((min[0] + 0.5 * dx + float(x) * dx, min[1] + 0.5 * dy + float(y) * dy, min[2] + 0.5 * dz + float(z) * dz))
#         idx += 1
#
#   i = 0
#   idx = 0
#   within = 0
#   invalid = 0
#   part1 = array([[-147.4, -2.842, 60.5], [-147.4, -2.842, -60.5], [-147.4, -51.16, 60.5], [-147.4, -51.16, -60.5], [250., -2.842, 45.], [250., -2.842, -45.], [250., -51.16, 45.], [250., -51.16, -45.]])
#   p1 = mkhull(part1)
#   for z in range(-overhead_z_l, nz + overhead_z_u):
#     for y in range(-overhead_y_l, ny + overhead_y_u):
#       for x in range(-overhead_x_l, nx + overhead_x_u):
#         coord = out[idx]
#         # test vs box
#         design = all(dot (p1.A, coord) <= p1.b)
#         nondesign = False
#         e = Element()
#         # if s1 < 0 point is out of the convex hull
#         # print 's1_schleife = ' + repr(s1)
#         if design:
#           within += 1
#           if not valid_position(coord, coords):
#             invalid += 1
#             e.region = 'void'
#           else:
#             e.region = 'design'
#         else:
# #           if nondesign:
# #             if not valid_position(coord, coords):
# #               invalid += 1
# #               e.region = 'void'
# #             else:
# #               e.region = 'nondesign'
# #           else:
#           e.region = 'void'
#         e.type = HEXA8
#         e.stiff1 = 1.
#         e.stiff2 = 1.
#         e.stiff3 = 1.
#         e.rotX = None if angles is None else 0.
#         e.rotY = None if angles is None else 0.
#         e.rotZ = None if angles is None else 0.
#         ll = (nnx + 1) * (nny + 1) * (z + overhead_z_l) + (nnx + 1) * (y + overhead_y_l) + (x + overhead_x_l)  # lowerleft
#         e.nodes = ((ll + (nnx + 1) * (nny + 1), ll + (nnx + 1) * (nny + 1) + nnx + 1, ll + (nnx + 1) * (nny + 1) + nnx + 1 + 1, ll + (nnx + 1) * (nny + 1) + 1, ll, ll + nnx + 1, ll + nnx + 1 + 1, ll + 1))
#         mesh.elements.append(e)
#         idx += 1
#
#   if grad <> 'nearest':
#     print str(within) + ' elements out of ' + str(len(out)) + ' in convex hull'
#   if invalid > 0:
#     print str(invalid) + ' elements out of ' + str(len(out)) + ' are considered invalid (shall not be visualized)'
#   return mesh


