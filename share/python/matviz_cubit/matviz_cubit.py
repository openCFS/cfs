#!/usr/bin/env python3

# This is a collection of functions to create unit cells (e.g. 2D triangles or 3D crosses)
# or dehomogenize macroscopic results with Coreform Cubit and mesh these geometries.
# The file provides similar functionality as matviz_2d.py.
# However, the functions have to be called from Cubit's Python command line.
# For examples, see dehomogenize_triangle.py and dehomogenize_cross_3d.py.


import cubit
import time
import numpy as np
import matviz_2d

# if param is given, create a parallelogram with two triangular holes and repeat in each dimension
# or dehomogenize a macroscopic result
# @param coords - (elem centers, min_bb, max_bb, elem_dim)
# @param design - design to dehomogenize else None
# @param grad - interpolation method (linear|nearest)
# @param samples - number of samples in each dimension for dehomogenization
# @param thres - threshold for design, currently not implemented
# @param equilateral - if triangles should be equilateral (base 1, height sqrt(3)/2 or isosceles (base 1, height 1)
# @param radius - radius for rounded corners (can be None)
# @param savefile - name of file to save geometry to
# @param param - volume for generation of unit cell (None for dehomogenization)
# @param repetitions - number of repetitions of unit cell
def show_triangle_grad(coords, design, grad, samples, thres, equilateral=True, radius=None, savefile=None, param=None, repetitions=1):

  # gives the angles to draw a circle segment
  def get_angles_for_chord(x, y, idx, angle, boundary):
    mid_circ_angle = (90 - angle*180/np.pi / 2) * 2
    # angles for arc
    if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
      if idx == 0:
        start, end = 90, 90 + mid_circ_angle
      elif idx == 1:
        start, end = 270 - mid_circ_angle/2, 270 + mid_circ_angle/2
      elif idx == 2:
        start, end = 90 - mid_circ_angle, 90
    else:
      if idx == 0:
        start, end = 270 - mid_circ_angle, 270
      elif idx == 1:
        start, end = 90 - mid_circ_angle/2, 90 + mid_circ_angle/2
      elif idx == 2:
        start, end = 270, 270 + mid_circ_angle
    if boundary:
      if x == 0:
        if y % 2 == 0:
          start, end = 90 - mid_circ_angle, 90
        else:
          start, end = 270, mid_circ_angle - 90
      else:
        if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
          start, end = 270 - mid_circ_angle, 270
        else:
          start, end = 90, 90 + mid_circ_angle
    return start, end

  # draws circle segments
  def draw_chords(x, y, tupels, tupel_ids, radius, innerTriangle):
    triaIds = np.zeros((3,))
    circIds = np.zeros((3,))
    for idx in range(3):
      p = tupels[idx]
      p1 = tupels[(idx+1)%3]
      p2 = tupels[(idx+2)%3]

      angle = np.pi/3.0

      start, end = get_angles_for_chord(x, y, idx, angle, False)

      # midpoint of arc
      vec = (p1+p2) / 2 - p
      mid_circ = p + vec/np.linalg.norm(vec) * radius / np.sin(angle/2)

      # arc start and end points
      a1 = np.zeros((2,))
      a2 = np.zeros((2,))
      a1[0] = mid_circ[0] + radius * np.cos(start/180*np.pi)
      a1[1] = mid_circ[1] - radius * np.sin(start/180*np.pi)
      a2[0] = mid_circ[0] + radius * np.cos(end/180*np.pi)
      a2[1] = mid_circ[1] - radius * np.sin(end/180*np.pi)

      cubit.silent_cmd('create vertex {} {} 0'.format( *a1 ))
      n4 = cubit.get_last_id("vertex")
      cubit.silent_cmd('create vertex {} {} 0'.format( *a2 ))
      n5 = cubit.get_last_id("vertex")
      cubit.silent_cmd('create surface vertex {} {} {}'.format(tupel_ids[idx], n4, n5))
      triaIds[idx] = cubit.get_last_id("surface")

      cubit.silent_cmd('create vertex {} {} 0'.format( *mid_circ ))
      n6 = cubit.get_last_id("vertex")
      cubit.silent_cmd('create surface circle radius {} zplane '.format(radius))
      circIds[idx] = cubit.get_last_id("surface")
      cubit.silent_cmd('move surface {} location vertex {} include_merged '.format(circIds[idx], n6))
      cubit.silent_cmd('delete vertex {} '.format(n6))

    cubit.silent_cmd('subtract surface {} {} {} from surface {} '.format(*triaIds, innerTriangle))
    cubit.silent_cmd('unite surface {} {} {} {}'.format(cubit.get_last_id("surface"), *circIds))

    return cubit.get_last_id("surface")


  start = time.time()
  
  assert(param is not None or design is not None)

  # not tested otherwise
  assert(equilateral == True)

  if param is None:
    centers, min_bb, max_bb, elem = coords

    # set size dx/dy/dz of one cell
    if samples is not None:
      dx = (max_bb[0] - min_bb[0]) / samples[0]
      dy = (max_bb[1] - min_bb[1]) / samples[1]
      elem = np.array([dx, dy, 0])

    elem = np.array([elem[0], elem[0] * np.sqrt(3)/2, 0])

    nx = int((max_bb[0]-min_bb[0]) / elem[0])
    ny = int(np.round((max_bb[1]-min_bb[1]) / elem[1]))
    # recalculate elem size, such that we always draw complete triangles in y-direction
    elem[1] = (max_bb[1]-min_bb[1]) / ny
  else:
    nx = repetitions
    ny = repetitions

  # edge length of outer triangle
  length = 1/nx

  # same conversion as in matviz_2d::show_triangle_grad
  rad = radius * length/np.sqrt(3)/4

  # this is a correction to match the radius from cubit_catalogue
  # 5 was number of cells in homogenization (nx)
  rad = rad / 8 * 5

  # same rounding as in matviz_2d::show_triangle_grad
  rad = np.round(rad * 1000) / 1000

  if param is None:
    print('Not cutting hole, if relative volume above {:.6f}.'.format(1 - np.pi*rad**2 / (np.sqrt(3)/4*length**2)))

  cubit.silent_cmd('reset')

  # draw outer shape
  if param is None:
    # draw rectangle
    width, height = max_bb[0]-min_bb[0], max_bb[1]-min_bb[1]
    cubit.silent_cmd('create surface rectangle width {} height {} zplane'.format(width, height))
    cubit.silent_cmd('move Surface {} x {} y {} z 0 include_merged '.format(cubit.get_last_id("surface"), width/2, height/2))
  else:
    nodes = np.zeros((3,2))
    nodes[0] = (0,0)
    nodes[1] = (1,0)
    nodes[2] = (0.5,np.sqrt(3)/2)

    node_ids = np.zeros((3,))
    cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[0] ))
    node_ids[0] = cubit.get_last_id("vertex")
    cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[1] ))
    node_ids[1] = cubit.get_last_id("vertex")
    cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[2] ))
    node_ids[2] = cubit.get_last_id("vertex")
    cubit.silent_cmd('create surface parallelogram vertex {} {} {}'.format(node_ids[2], node_ids[0], node_ids[1]))
  outerShape = cubit.get_last_id("surface")

  nx *= 2
  # draw additional triangles
  if param is not None:
    # for parallelogram we need a lot more.
    # however we will skip the additional ones
    nx += ny-1
  else:
    # two more wich will be cut later
    nx += 2

  triangles = np.zeros((nx*ny,))
  triangleidx = 0

  if param is None:
    # calculate new element centers
    centers_new = []
    for y in range(ny):
      for x in range(nx):
        # -length/2 will add one column too much. we need this to get
        # the correct orientation of triangles
        center = np.array([-length/2 + x * length/2, np.sqrt(3)/4 * length + y * np.sqrt(3)/2 * length])
        centers_new.append(center)
        # kind of integration points -> should be removed and graded triangle bars should be used
        if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
          #    2        4--2           2--5
          #   / \       |   \         /   |
          #  /   \      |    \       /    |
          # 1-----3     1-----3     1-----3
          centers_new.append(center + [-elem[0]/4.0, -elem[1]/4.0])
          centers_new.append(center + [+elem[0]/4.0, -elem[1]/4.0])
          centers_new.append(center + [         0.0, +elem[1]/4.0])
        else:
          # upside down
          # 1-----3     1-----3     1-----3
          #  \   /      |    /       \    |
          #   \ /       |   /         \   |
          #    2        4--2           2--5
          centers_new.append(center + [-elem[0]/4.0, +elem[1]/4.0])
          centers_new.append(center + [+elem[0]/4.0, +elem[1]/4.0])
          centers_new.append(center + [         0.0, -elem[1]/4.0])

    # interpolate original data at new centers
    c, v = matviz_2d.convert_single_data_interpolation_input(centers, design, None)
    ip_data = matviz_2d.ip.griddata(c, v, centers_new, grad, -1.0)
    # any interpolation but nearest neighbor can only interpolate in the convex hull
    ip_near = matviz_2d.ip.griddata(c, v, centers_new, 'nearest') if grad != 'nearest' else None

  for y in range(ny):
    for x in range(nx):
      if param is not None:
        area = np.min((1.0, np.max((0.0, param))))
        # skip to get parallelogram
        if x < y or x > nx - ny + y:
          continue
        # center of triangle
        mid = np.array([length/2 + x * length/2, np.sqrt(3)/4 * length + y * np.sqrt(3)/2 * length])
      else:
        if x == 0:
          # to get the right orientation, we added one column too much in centers_new
          continue

        # if the value is -1 we use the nearest interpolation
        # mid IS NOT THE MIDPOINT OF THE TRIANGLE BUT OF elem (= enclosing rectangle)
        mid, area = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+0)
        _, area1 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+1)
        _, area2 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+2)
        _, area3 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+3)
        # kind of integration -> should be removed and graded triangle bars should be used
        area = (3*area + area1+area2+area3)/6
        area = area[0]

      area *= np.sqrt(3)/4*length**2

      if area < np.sqrt(3)/4*length**2 - np.pi*rad**2:
        # when we want a parallelogram for homogenization, the design is constant and thus all holes are the same.
        # we only draw one hole and copy + rotate + translate it
        if param is None or (x == 0 and y == 0):
          # without round corners
          # area = 3*param*(g - sqrt(3)*param)
          # param = (np.sqrt(3)*length - np.sqrt( 3*length**2 - 4*np.sqrt(3)*area )) / 6

          # with round corners
          # area = 3*param*(g - sqrt(3)*param) + (3*sqrt(3) - pi)*rad^2
          p = (np.sqrt(3)*length - np.sqrt( 3*length**2 + 4*(9-np.sqrt(3)*np.pi)*rad**2 - 4*np.sqrt(3)*area )) / 6

          # nodes of inner triangle
          nodes = np.zeros((3,2))
          nodes[0,0] = mid[0] - length/2 + np.sqrt(3)*p
          nodes[1,0] = mid[0]
          nodes[2,0] = mid[0] + length/2 - np.sqrt(3)*p
          if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
            #    2        4--2           2--5
            #   / \       |   \         /   |
            #  /   \      |    \       /    |
            # 1-----3     1-----3     1-----3
            nodes[0,1] = mid[1] - np.sqrt(3)/4*length + p
            nodes[1,1] = mid[1] + np.sqrt(3)/4*length - 2*p
          else:
            # upside down
            # 1-----3     1-----3     1-----3
            #  \   /      |    /       \    |
            #   \ /       |   /         \   |
            #    2        4--2           2--5
            nodes[0,1] = mid[1] + np.sqrt(3)/4*length - p
            nodes[1,1] = mid[1] - np.sqrt(3)/4*length + 2*p
          nodes[2,1] = nodes[0,1]

          # The triangle approach can lead to very short edges, if the rounded 
          # corners almost form a circle. This causes problems when meshing.
          # Thus, we draw a circle instead, if the area of the hole is smaller
          # than 3.5 times the area of a circle with radius rad.
          threshold = 3.5
          if np.sqrt(3)/4*length**2 - area < np.pi*rad**2*threshold:
            # midpoint of incircle of triangle
            edgelengths = np.linalg.norm(nodes, axis=1)
            incirc_mid = np.dot(np.transpose(nodes), edgelengths) / sum(edgelengths)
            
            # radius of incircle
            incirc_rad = np.sqrt((np.sqrt(3)/4*length**2 - area) / np.pi)
            
            # draw circle and move to the correct position
            cubit.silent_cmd('create surface circle radius {} zplane'.format(incirc_rad))
            innerCut = cubit.get_last_id("surface")
            cubit.silent_cmd('move Surface {} location {} {} 0 include_merged'.format(innerCut, *incirc_mid))
  
          else: # draw a triangle with rounded corners
            # draw inner triangle
            node_ids = np.zeros((3,))
            cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[0] ))
            node_ids[0] = cubit.get_last_id("vertex")
            cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[1] ))
            node_ids[1] = cubit.get_last_id("vertex")
            cubit.silent_cmd('create vertex {} {} 0'.format( *nodes[2] ))
            node_ids[2] = cubit.get_last_id("vertex")
            cubit.silent_cmd('create surface vertex {} {} {}'.format(node_ids[0], node_ids[1], node_ids[2]))
            innerTriangle = cubit.get_last_id("surface")
  
            # draw rounded corners
            if rad > 0:
              innerCut = draw_chords(x+1, y, nodes, node_ids, rad, innerTriangle)
            else:
              innerCut = innerTriangle

          triangles[triangleidx] = innerCut
          triangleidx += 1
        else:
          # we copy + rotate + translate the first hole
          if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
            rotation_angle = 0
            translation_x = mid[0] - length/2
            translation_y = mid[1] - np.sqrt(3)/4 * length
          else:
            rotation_angle = 180
            translation_x = mid[0] + length/2
            translation_y = mid[1] + np.sqrt(3)/4 * length
          cubit.silent_cmd('surface {} copy rotate {} about z nomesh'.format(triangles[0], rotation_angle))
          id = cubit.get_last_id("surface")
          cubit.silent_cmd('move Surface {} x {} y {} z 0 include_merged '.format(id, translation_x, translation_y))
          triangles[triangleidx] = cubit.get_last_id("surface")
          triangleidx += 1
      else:
        pass

  mystring = ''
  for idx in triangles:
    if idx > 0:
      mystring += '{:.0f} '.format(idx)

  # cut inner triangles from outerShape
  cubit.silent_cmd('subtract surface {} from surface {} imprint '.format(mystring, outerShape))

  if savefile is not None:
    cubit.cmd('save cub5 "{}.cub5" overwrite journal'.format(savefile))

  shape = cubit.get_last_id("surface")

  name_regions_and_nodes(shape, param)

  end = time.time()
  print('Time for creating geometry: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))))

  return shape


# names regions and nodes of a geometry
# @param shape - id of geometry
# @param param - if None, macroscopic naming is applied, if not None, unit cell will be named
def name_regions_and_nodes(shape, param):
  # name region
  cubit.silent_cmd('block 1 add surface {} '.format(shape))
  cubit.silent_cmd('block 1 name "mech"')

  # name boundary nodes
  if param is not None:
    cubit.silent_cmd('nodeset 1 add curve 1 ')
    cubit.silent_cmd('nodeset 2 add curve 2 ')
    cubit.silent_cmd('nodeset 3 add curve 3 ')
    cubit.silent_cmd('nodeset 4 add curve 4 ')
  else:
    bb = cubit.get_bounding_box("surface", shape)

    all_curves = cubit.get_entities("curve")

    # find curves belonging to boundaries
    left = []
    right = []
    bottom = []
    top = []
    for curve in all_curves:
      start_point = cubit.curve(curve).position_from_fraction(0)
      end_point = cubit.curve(curve).position_from_fraction(1)
      if abs(start_point[0]-bb[0]) < 1e-12 and abs(end_point[0]-bb[0]) < 1e-12:
        left.append(curve)
      if abs(start_point[0]-bb[1]) < 1e-12 and abs(end_point[0]-bb[1]) < 1e-12:
        right.append(curve)
      if abs(start_point[1]-bb[3]) < 1e-12 and abs(end_point[1]-bb[3]) < 1e-12:
        bottom.append(curve)
      if abs(start_point[1]-bb[4]) < 1e-12 and abs(end_point[1]-bb[4]) < 1e-12:
        top.append(curve)

    # generate string from list
    left_ids = ' '.join('{:d}'.format(x) for x in left)
    right_ids = ' '.join('{:d}'.format(x) for x in right)
    bottom_ids = ' '.join('{:d}'.format(x) for x in bottom)
    top_ids = ' '.join('{:d}'.format(x) for x in top)

    print('Adding {} curve(s) to westy.'.format(len(left)))
    print('Adding {} curve(s) to southy.'.format(len(bottom)))
    print('Adding {} curve(s) to easty.'.format(len(right)))
    print('Adding {} curve(s) to northy.'.format(len(top)))

    cubit.silent_cmd('nodeset 1 add curve {} '.format(left_ids))
    cubit.silent_cmd('nodeset 2 add curve {} '.format(bottom_ids))
    cubit.silent_cmd('nodeset 3 add curve {} '.format(right_ids))
    cubit.silent_cmd('nodeset 4 add curve {} '.format(top_ids))

  # west, south, east and north are cubit identifiers, thus we use -y
  cubit.silent_cmd('nodeset 1 name "westy"')
  cubit.silent_cmd('nodeset 2 name "southy"')
  cubit.silent_cmd('nodeset 3 name "easty"')
  cubit.silent_cmd('nodeset 4 name "northy"')

  cubit.silent_cmd('nodeset 5 add vertex 1')
  cubit.silent_cmd('nodeset 5 name "north_east"')
  cubit.silent_cmd('nodeset 6 add vertex 2')
  cubit.silent_cmd('nodeset 6 name "north_west"')
  cubit.silent_cmd('nodeset 7 add vertex 3')
  cubit.silent_cmd('nodeset 7 name "south_west"')
  cubit.silent_cmd('nodeset 8 add vertex 4')
  cubit.silent_cmd('nodeset 8 name "south_east"')


# mesh a given geometry
# @param shape - id of geometry
# @param meshsize - approximate size of FE
# @param filename - filename of meshfile
def mesh_shape(shape, meshsize, filename):
  start = time.time()

  block_id = cubit.get_block_id("surface", shape)
  cubit.silent_cmd('block {} element type TRI3'.format(block_id))

  # mesh shape
  cubit.silent_cmd('surface {} size {} '.format(shape, meshsize))
  cubit.silent_cmd('surface {} scheme trimesh'.format(shape))
  cubit.silent_cmd('Trimesher geometry sizing on')
  cubit.silent_cmd('mesh surface {}'.format(shape))

  cubit.cmd('export ansys cdb "{}.cdb"  overwrite '.format(filename))

  bb = cubit.get_bounding_box("surface", shape)
  rel_volume = cubit.get_surface_area(shape) / bb[2] / bb[5]
  rel_meshed_volume = cubit.get_meshed_volume_or_area("surface",[shape]) / bb[2] / bb[5] 

  end = time.time()
  print('Time for meshing: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))))
  print('Relative Surface Area: {}'.format(rel_volume))
  print('Relative Meshed Area: {} ({} triangles)\n'.format(rel_meshed_volume, cubit.get_tri_count()))

  with open('{}.log'.format(filename), 'a+') as fid:
    print(filename, file=fid)
    print('Relative Surface Area: {}'.format(rel_volume), file=fid)
    print('Relative Meshed Area: {} ({} triangles)\n'.format(rel_meshed_volume, cubit.get_tri_count()), file=fid)
    print('Time for meshing: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))), file=fid)

    print('Number of nodes in westy:  {}'.format(cubit.get_nodeset_node_count(1)), file=fid)
    print('Number of nodes in southy: {}'.format(cubit.get_nodeset_node_count(2)), file=fid)
    print('Number of nodes in easty:  {}'.format(cubit.get_nodeset_node_count(3)), file=fid)
    print('Number of nodes in northy: {}'.format(cubit.get_nodeset_node_count(4)), file=fid)

if __name__ == "__main__":
  print('Hello!')
