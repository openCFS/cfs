#!/usr/bin/env python3

# This is a collection of functions to create unit cells (e.g. 2D triangles)
# or dehomogenize macroscopic results with Coreform Cubit and mesh these geometries.
# The file provides similar functionality as matviz_2d.py.
# However, the functions have to be called from Cubit's Python command line.
# For examples, see dehomogenize_triangle.py and dehomogenize_cross_3d.py.


import os, sys, subprocess

try:
  import cubit
except ImportError:
  ########################################################################
  # this section will set the correct environment paths for cubit.       #
  # this is only necessary, if you import cubit into python outside of   #
  # the cubit gui/journal editor. paths should be found automatically.   #
  # if not or paths are missing, change the following paths and set the  #
  # environment variables below.                                         #
  #                                                                      #
  #    PATHS HAVE TO BE SET BEFORE OTHER IMPORTS (especially numpy)!     #

  cubit_search_path = sys.executable
  pillow_search_path = sys.path

  # find the path of a file
  def find_path(file, search_path=None):
    out = subprocess.run(["locate", "-b", "\\" + file], capture_output=True)
    if out.returncode == 0:
      for path in out.stdout.split():
        if not search_path or (sp in path.decode("utf-8") for sp in search_path):
          # returns first occurence
          return os.path.dirname(path.decode("utf-8"))
    else:
      for root, _, files in os.walk(search_path):
        if file in files:
          return root

  CUBIT_PATH = find_path('coreform_cubit', cubit_search_path)
  PILLOW_PATH = find_path('_imaging.cpython-38-x86_64-linux-gnu.so', pillow_search_path)

  assert(CUBIT_PATH)
  assert(PILLOW_PATH)

  # this will evaluate to false, if matviz_cubit is imported inside of cubit
  # gui/journal editor and if paths are set before starting cubit, see wiki
  if not CUBIT_PATH in os.environ['LD_LIBRARY_PATH']:
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.path.join(CUBIT_PATH) + ':' + os.path.join(CUBIT_PATH, 'plugins/mediaservice') + ':' + env['LD_LIBRARY_PATH']
    # env['PYTHONPATH'] = '/home/daniel/code/cfs/share/python:/home/daniel/code/cfs/share/python/matviz_cubit:/home/daniel/code/cfs/share/python/trelis:' + os.environ['PYTHONPATH']
    env['PYTHONPATH'] = os.path.join(PILLOW_PATH) + ':' + os.path.join(CUBIT_PATH) + ':' + os.path.join(CUBIT_PATH, 'python3/lib/python3.8/lib-dynload') + ':' + os.path.join(CUBIT_PATH, 'python3/lib/python3.8/site-packages') + ':' + os.environ['PYTHONPATH']
    env['PATH'] = '/home/daniel/anaconda3/lib/:' + env['PATH']
    env['HDF5_DISABLE_VERSION_CHECK']='2'
    try:
      # restart python interpreter
      os.execve(sys.executable, [sys.executable] + [sys.argv[0]], env)
    except Exception as exc:
      print('Failed re-exec:', exc)
      sys.exit(1)

  ########################################################################

  import cubit
  cubit_version = cubit.get_version()

  # add Cubit libraries to your path
  try:
    sys.path.index('/opt/Coreform-Cubit-{}/bin'.format(cubit_version))
  except ValueError:
    sys.path.append('/opt/Coreform-Cubit-{}/bin'.format(cubit_version))

  # start cubit - this step is key if not working inside cubit gui/journal editor!
  cubit.init(['cubit', '-nographics', '-nojournal', '-batch', '-noecho', '-information', 'off', '-warning', 'off'])

print('running on Cubit version {}'.format(cubit.get_version()))

import time
import numpy as np
#import matviz_2d
from scipy.spatial import ConvexHull


possible_nodesets = {11:'lefty', 12:'bottomy', 13:'righty', 14:'topy',
                     21:'bottom_left', 22:'bottom_right', 23:'top_left', 24:'top_right'}


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

  if design is not None:
    homogeneousDesign = abs(max(design) - min(design)) < 1e-12
  else:
    homogeneousDesign = False

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

    # edge length of outer triangle
    length = (max_bb[0]-min_bb[0])/nx

  else:
    nx = repetitions
    ny = repetitions
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
  else:
    print('param = {:.6f}'.format(param))

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
    # calculate new element centers as gauss integration points
    centers_new = []
    for y in range(ny):
      for x in range(nx):
        # -length/2 will add one column too much. we need this to get
        # the correct orientation of triangles
        #center = np.array([-length/2 + x * length/2, 1/2 * np.sqrt(3)/2 * length + y * np.sqrt(3)/2 * length])
        center = np.array([-length/2 + x * length/2, 1/3 * np.sqrt(3)/2 * length + y * np.sqrt(3)/2 * length])
        centers_new.append(center)
        # kind of integration points -> should be removed and graded triangle bars should be used
        if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
          #    2        4--2           2--5
          #   / \       |   \         /   |
          #  /   \      |    \       /    |
          # 1-----3     1-----3     1-----3
          #centers_new.append(center + [-elem[0]/4.0, -elem[1]/4.0])
          #centers_new.append(center + [+elem[0]/4.0, -elem[1]/4.0])
          #centers_new.append(center + [         0.0, +elem[1]/4.0])
          centers_new.append(center + [      0.0, + 2/5 * 2/3 * np.sqrt(3)/2 * length])
          centers_new.append(center + [-length/5, - 1/5 * np.sqrt(3)/2 * length])
          centers_new.append(center + [+length/5, - 1/5 * np.sqrt(3)/2 * length])
        else:
          # upside down
          # 1-----3     1-----3     1-----3
          #  \   /      |    /       \    |
          #   \ /       |   /         \   |
          #    2        4--2           2--5
          #centers_new.append(center + [-elem[0]/4.0, +elem[1]/4.0])
          #centers_new.append(center + [+elem[0]/4.0, +elem[1]/4.0])
          #centers_new.append(center + [         0.0, -elem[1]/4.0])
          centers_new.append(center + [      0.0, - 2/5 * 2/3 * np.sqrt(3)/2 * length])
          centers_new.append(center + [-length/5, + 1/5 * np.sqrt(3)/2 * length])
          centers_new.append(center + [+length/5, + 1/5 * np.sqrt(3)/2 * length])

    # interpolate original data at new centers
    c, v = matviz_2d.convert_single_data_interpolation_input(centers, design, None)
    ip_data = matviz_2d.ip.griddata(c, v, centers_new, grad, -1.0)
    # any interpolation but nearest neighbor can only interpolate in the convex hull
    ip_near = matviz_2d.ip.griddata(c, v, centers_new, 'nearest') if grad != 'nearest' else None

  holes = np.zeros((nx*ny,2))
  mid_idx = 0
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
        mid = np.array([-length/2 + x * length/2, 1/2 * np.sqrt(3)/2 * length + y * np.sqrt(3)/2 * length])
        
        # gauss integration -> should be removed and graded triangle bars should be used
        _, area = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+0)
        _, area1 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+1)
        _, area2 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+2)
        _, area3 = matviz_2d.get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+3)
        area = -27/48 * area + 25/48 * (area1+area2+area3)
        area = area[0]

      area *= np.sqrt(3)/4*length**2

      if area < np.sqrt(3)/4*length**2 - np.pi*rad**2:
        holes[mid_idx,:] = mid
        mid_idx += 1
        # when we want a parallelogram for homogenization (param not None), all holes are the same size.
        # same is true for homogeneous design.
        # in these cases we only draw one hole and copy + rotate + translate it
        if (param is not None and x == 0 and y == 0) or (param is None and (not homogeneousDesign or (x == 1 and y == 0))):
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
          # than threshold times the area of a circle with radius rad.
          threshold = 1.2
          if np.sqrt(3)/4*length**2 - area < np.pi*rad**2*threshold:
            # midpoint of incircle of triangle
            edgelengths = np.linalg.norm(nodes[[1,2,0],:]-nodes[[2,0,1],:], axis=1)
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
          if param is not None:
            if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
              rotation_angle = 0
              translation_x = mid[0] - length/2
              translation_y = mid[1] - np.sqrt(3)/4 * length
            else:
              rotation_angle = 180
              translation_x = mid[0] + length/2
              translation_y = mid[1] + np.sqrt(3)/4 * length
          else:
            translation_x = mid[0]
            if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
              rotation_angle = 180
              translation_y = mid[1] + np.sqrt(3)/4 * length
            else:
              rotation_angle = 0
              translation_y = mid[1] - np.sqrt(3)/4 * length
          cubit.silent_cmd('surface {} copy rotate {} about z nomesh'.format(triangles[0], rotation_angle))
          id = cubit.get_last_id("surface")
          cubit.silent_cmd('move Surface {} x {} y {} z 0 include_merged '.format(id, translation_x, translation_y))
          triangles[triangleidx] = cubit.get_last_id("surface")
          triangleidx += 1
      else:
        pass

  holes = holes[~np.all(holes == 0, axis=1)]

  mystring = ''
  for idx in triangles:
    if idx > 0:
      mystring += '{:.0f} '.format(idx)

  # cut inner triangles from outerShape
  if mystring:
    cubit.silent_cmd('subtract surface {} from surface {} imprint '.format(mystring, outerShape))

  if savefile is not None:
    cubit.cmd('save cub5 "{}.cub5" overwrite journal'.format(savefile))
    with open(savefile + '.hole', 'w+') as fid:
      np.savetxt(fid, holes, fmt='%f %f')

  shape = cubit.get_last_id("surface")

  #name_regions_and_nodes(shape, param)

  end = time.time()
  print('Time for creating geometry: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))))

  return shape, holes


def show_cross_3D_grad(coords, design, grad, samples, thres, radius=None, savefile=None, param=None, repetitions=1):
  radius = param
  height = 10
  
  cubit.silent_cmd('reset')
  
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height, radius))
  cyl1 = cubit.get_last_id("volume")
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height, radius))
  cyl2 = cubit.get_last_id("volume")
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height, radius))
  cyl3 = cubit.get_last_id("volume")
  
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height*np.sqrt(3), radius))
  cyl4 = cubit.get_last_id("volume")
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height*np.sqrt(3), radius))
  cyl5 = cubit.get_last_id("volume")
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height*np.sqrt(3), radius))
  cyl6 = cubit.get_last_id("volume")
  cubit.silent_cmd('create Cylinder height {} radius {} '.format(height*np.sqrt(3), radius))
  cyl7 = cubit.get_last_id("volume")
  
  cubit.silent_cmd('rotate Volume {} angle 90  about X include_merged '.format(cyl2))
  cubit.silent_cmd('rotate Volume {} angle 90  about Y include_merged '.format(cyl3))
  angle = np.arccos(1/np.sqrt(3))*180/np.pi
  cubit.silent_cmd('rotate Volume {} angle {}  about origin 0 0 0 direction 1 1 0 include_merged '.format(cyl4, angle))
  cubit.silent_cmd('rotate Volume {} angle {}  about origin 0 0 0 direction -1 1 0 include_merged '.format(cyl5, angle))
  cubit.silent_cmd('rotate Volume {} angle {}  about origin 0 0 0 direction 1 -1 0 include_merged '.format(cyl6, angle))
  cubit.silent_cmd('rotate Volume {} angle {}  about origin 0 0 0 direction -1 -1 0 include_merged '.format(cyl7, angle))
  
  cubit.silent_cmd('brick x {} '.format(height))
  cube = cubit.get_last_id("volume")
  cubit.silent_cmd('unite body {:d} {:d} {:d} {:d} {:d} {:d} {:d} '.format(cyl1, cyl2, cyl3, cyl4, cyl5, cyl6, cyl7))
  cubit.silent_cmd('intersect body {:d} {:d} '.format(cube, cyl1))
  shape = cyl1
  
  cubit.silent_cmd('list body {} geometry'.format(shape))
  
  for xx in range(repetitions):
    for yy in range(repetitions):
      for zz in range(repetitions):
        cubit.silent_cmd('Volume {}  copy move x {} y {} z {}'.format(shape, xx*height, yy*height, zz*height))
  
  cubit.silent_cmd('unite body all')
  shape = cubit.get_entities("volume")[0]
  
  return shape, []


# names regions and nodes of a geometry
# @param shape - id of geometry
# @param param - if None, macroscopic naming is applied, if not None, unit cell will be named
def name_regions_and_nodes(shape, param):

  is2d = cubit.get_volume_count() == cubit.get_surface_count()
  entity = "surface" if is2d else "volume"
  subentity = "curve" if is2d else "surface"

  # name region
  cubit.silent_cmd('block 1 add ' + entity + ' {} '.format(shape))
  cubit.silent_cmd('block 1 name "mech"')

  # name boundary nodes
  if param is not None:
    cubit.silent_cmd('nodeset 1 add curve 1 ')
    cubit.silent_cmd('nodeset 2 add curve 2 ')
    cubit.silent_cmd('nodeset 3 add curve 3 ')
    cubit.silent_cmd('nodeset 4 add curve 4 ')
    cubit.silent_cmd('nodeset 5 add vertex 1')
    cubit.silent_cmd('nodeset 6 add vertex 2')
    cubit.silent_cmd('nodeset 7 add vertex 3')
    cubit.silent_cmd('nodeset 8 add vertex 4')
  else:
    bb = cubit.get_bounding_box(entity, shape)

    all_subentities = cubit.get_entities(subentity)

    # find curves belonging to boundaries
    left = []
    right = []
    bottom = []
    top = []
    back = []
    front = []

    for subent in all_subentities:
      if is2d:
        start_point = cubit.curve(subent).position_from_fraction(0)
        end_point = cubit.curve(subent).position_from_fraction(1)
        if abs(start_point[0]-bb[0]) < 1e-12 and abs(end_point[0]-bb[0]) < 1e-12:
          left.append(subent)
        if abs(start_point[0]-bb[1]) < 1e-12 and abs(end_point[0]-bb[1]) < 1e-12:
          right.append(subent)
        if abs(start_point[1]-bb[3]) < 1e-12 and abs(end_point[1]-bb[3]) < 1e-12:
          bottom.append(subent)
        if abs(start_point[1]-bb[4]) < 1e-12 and abs(end_point[1]-bb[4]) < 1e-12:
          top.append(subent)

        bottom_left_id, top_left_id, bottom_right_id, top_right_id = None, None, None, None
        for curve in left:
          vertices = cubit.get_relatives("curve", curve, "vertex")
          for vertex in vertices:
            if bottom_left_id is None and abs(cubit.vertex(vertex).coordinates()[1]-bb[3]) < 1e-12:
              bottom_left_id = cubit.vertex(vertex).id()
            if top_left_id is None and abs(cubit.vertex(vertex).coordinates()[1]-bb[4]) < 1e-12:
              top_left_id = cubit.vertex(vertex).id()
        for curve in right:
          vertices = cubit.get_relatives("curve", curve, "vertex")
          for vertex in vertices:
            if bottom_right_id is None and abs(cubit.vertex(vertex).coordinates()[1]-bb[3]) < 1e-12:
              bottom_right_id = cubit.vertex(vertex).id()
            if top_right_id is None and abs(cubit.vertex(vertex).coordinates()[1]-bb[4]) < 1e-12:
              top_right_id = cubit.vertex(vertex).id()
    
        print('top_right_id: {}'.format(top_right_id))
        print('top_left_id: {}'.format(top_left_id))
        print('bottom_left_id: {}'.format(bottom_left_id))
        print('bottom_right_id: {}'.format(bottom_right_id))

      else:
        a = cubit.surface(subent).get_param_range_U()
        b = cubit.surface(subent).get_param_range_V()
        p1 = cubit.surface(subent).position_from_u_v(a[0],b[0])
        p2 = cubit.surface(subent).position_from_u_v(a[1],b[0])
        p3 = cubit.surface(subent).position_from_u_v(a[0],b[1])
        if abs(p1[0]-bb[0]) < 1e-12 and abs(p2[0]-bb[0]) < 1e-12 and abs(p3[0]-bb[0]) < 1e-12:
          left.append(subent)
        if abs(p1[0]-bb[1]) < 1e-12 and abs(p2[0]-bb[1]) < 1e-12 and abs(p3[0]-bb[1]) < 1e-12:
          right.append(subent)
        if abs(p1[1]-bb[3]) < 1e-12 and abs(p2[1]-bb[3]) < 1e-12 and abs(p3[1]-bb[3]) < 1e-12:
          bottom.append(subent)
        if abs(p1[1]-bb[4]) < 1e-12 and abs(p2[1]-bb[4]) < 1e-12 and abs(p3[1]-bb[4]) < 1e-12:
          top.append(subent)
        if abs(p1[2]-bb[6]) < 1e-12 and abs(p2[2]-bb[6]) < 1e-12 and abs(p3[2]-bb[6]) < 1e-12:
          back.append(subent)
        if abs(p1[2]-bb[7]) < 1e-12 and abs(p2[2]-bb[7]) < 1e-12 and abs(p3[2]-bb[7]) < 1e-12:
          front.append(subent)

    # generate string from list
    left_ids = ' '.join('{:d}'.format(x) for x in left)
    right_ids = ' '.join('{:d}'.format(x) for x in right)
    bottom_ids = ' '.join('{:d}'.format(x) for x in bottom)
    top_ids = ' '.join('{:d}'.format(x) for x in top)
    back_ids = ' '.join('{:d}'.format(x) for x in back)
    front_ids = ' '.join('{:d}'.format(x) for x in front)

    print('Adding {} '.format(len(left)) + subentity + '(s) to lefty.')
    print('Adding {} '.format(len(bottom)) + subentity + '(s) to bottomy.')
    print('Adding {} '.format(len(right)) + subentity + '(s) to righty.')
    print('Adding {} '.format(len(top)) + subentity + '(s) to topy.')
    if not is2d:
      print('Adding {} '.format(len(back)) + subentity + '(s) to backy.')
      print('Adding {} '.format(len(front)) + subentity + '(s) to fronty.')
  
    cubit.silent_cmd('nodeset 1 add ' + subentity + ' {} '.format(left_ids))
    cubit.silent_cmd('nodeset 2 add ' + subentity + ' {} '.format(bottom_ids))
    cubit.silent_cmd('nodeset 3 add ' + subentity + ' {} '.format(right_ids))
    cubit.silent_cmd('nodeset 4 add ' + subentity + ' {} '.format(top_ids))
    if not is2d:
      cubit.silent_cmd('nodeset 5 add ' + subentity + ' {} '.format(back_ids))
      cubit.silent_cmd('nodeset 6 add ' + subentity + ' {} '.format(front_ids))

    if is2d:
      cubit.silent_cmd('nodeset 5 add vertex {}'.format(bottom_left_id))
      cubit.silent_cmd('nodeset 6 add vertex {}'.format(bottom_right_id))
      cubit.silent_cmd('nodeset 7 add vertex {}'.format(top_left_id))
      cubit.silent_cmd('nodeset 8 add vertex {}'.format(top_right_id))

  # west, south, east and north are cubit identifiers, thus we use -y
  cubit.silent_cmd('nodeset 1 name "lefty"')
  cubit.silent_cmd('nodeset 2 name "bottomy"')
  cubit.silent_cmd('nodeset 3 name "righty"')
  cubit.silent_cmd('nodeset 4 name "topy"')
  if not is2d:
    cubit.silent_cmd('nodeset 5 name "backy"')
    cubit.silent_cmd('nodeset 6 name "fronty"')

  if is2d:
    cubit.silent_cmd('nodeset 5 name "bottom_left"')
    cubit.silent_cmd('nodeset 6 name "bottom_right"')
    cubit.silent_cmd('nodeset 7 name "top_left"')
    cubit.silent_cmd('nodeset 8 name "top_right"')


# find matching surfaces, curves and vertices for periodic mesh
def get_pairs(shape):
  left = cubit.get_nodeset_surfaces(1)
  bottom = cubit.get_nodeset_surfaces(2)
  right = cubit.get_nodeset_surfaces(3)
  top = cubit.get_nodeset_surfaces(4)
  back = cubit.get_nodeset_surfaces(5)
  front = cubit.get_nodeset_surfaces(6)

  bb = cubit.get_bounding_box("volume", shape)

  periodic_sides = [[0, left, right, bb[2]], [1, bottom, top, bb[5]], [2, back, front, bb[8]]]

  pairs = []

  for axis_idx, master_list, slave_list, distance in periodic_sides:
    for master in master_list:
      # get a point on the master surface
      msurf = cubit.surface(master)
      a = msurf.get_param_range_U()
      b = msurf.get_param_range_V()
      p = list(msurf.position_from_u_v( (a[0]+a[1])/2, (b[0]+b[1])/2 ))
      # project the point to the other side
      p[axis_idx] += distance
      for slave in slave_list:
        ssurf = cubit.surface(slave)
        # check if projected point is on slave surface
        if ssurf.point_containment(p):
          mcurves = cubit.get_relatives("surface", master, "curve")
          scurves = cubit.get_relatives("surface", slave, "curve")
          assert(len(mcurves) == len(scurves))
          if len(mcurves) == 1:
            # for circles we have only one curve and one vertex
            mvertex = cubit.get_relatives("curve", mcurves[0], "vertex")
            svertex = cubit.get_relatives("curve", scurves[0], "vertex")
            #p1 = cubit.vertex(mvertex[0]).coordinates()
            #p2 = cubit.vertex(svertex[0]).coordinates()
            #diff = [y - x for x, y in zip(p1,p2)]
            #assert(abs(diff[0]-distance) < 1e-12 and abs(diff[1]) < 1e-12 and abs(diff[2]) < 1e-12)
            pairs.append([master, slave, mcurves[0], scurves[0], mvertex[0], svertex[0]])
            break #for slave in slave_list
          for mcurve in mcurves:
            # for more curves, we have to find matching curve and vertex
            mvertices = cubit.get_relatives("curve", mcurve, "vertex")
            mstart_point = np.array(cubit.vertex(mvertices[0]).coordinates())
            mend_point = np.array(cubit.vertex(mvertices[1]).coordinates())
            mstart_point[axis_idx] += distance
            mend_point[axis_idx] += distance
            for scurve in scurves:
              svertices = cubit.get_relatives("curve", scurve, "vertex")
              sstart_point = np.array(cubit.vertex(svertices[0]).coordinates())
              send_point = np.array(cubit.vertex(svertices[1]).coordinates())
              n1 = np.linalg.norm(mstart_point-sstart_point)
              n2 = np.linalg.norm(mstart_point-send_point)
              n3 = np.linalg.norm(mend_point-sstart_point)
              n4 = np.linalg.norm(mend_point-send_point)
              # mcurve and scurve are either congruent or mirrored
              if (abs(n1) < 1e-12 and abs(n4) < 1e-12) or (abs(n2) < 1e-12 and abs(n3) < 1e-12):
                if abs(n1) < 1e-12:
                  pairs.append([master, slave, mcurve, scurve, mvertices[0], svertices[0]])
                else:
                  pairs.append([master, slave, mcurve, scurve, mvertices[0], svertices[1]])
                break #for scurve in scurves
            # this is a trick to break both for loops. 
            # else case is only executed, if the inner loop does not break
            else:
              continue
            break #for mcurve in mcurves
        else:
          continue
        break #for slave in slave_list

  return pairs


# mesh a given geometry
# @param shape - id of geometry
# @param meshsize - approximate size of FE
# @param filename - filename of meshfile
def mesh_shape(shape, meshsize, filename):
  
  is2d = cubit.get_volume_count() == cubit.get_surface_count()
  entity = "surface" if is2d else "volume"
  subentity = "curve" if is2d else "surface"

  start = time.time()

  # calculate convex hull of shape
  vertices_ids = cubit.get_relatives(entity, shape, "vertex")
  vertices = np.zeros((len(vertices_ids), 2 if is2d else 3))
  for ii, vertex_id in enumerate(vertices_ids):
    vertex = cubit.vertex(vertex_id)
    vertices[ii,:] = vertex.coordinates()[:2 if is2d else 3]
  hull = ConvexHull(vertices)

  if is2d:
    block_id = cubit.get_block_id("surface", shape)
    if block_id > 0:
      cubit.silent_cmd('block {} element type TRI3'.format(block_id))

    # mesh shape
    cubit.silent_cmd('surface {} size {} '.format(shape, meshsize))
    cubit.silent_cmd('surface {} scheme trimesh'.format(shape))
    cubit.silent_cmd('Trimesher geometry sizing on')
    cubit.silent_cmd('mesh surface {}'.format(shape))

    nelem = cubit.get_tri_count()
    #bb = cubit.get_bounding_box("surface", shape)
    abs_volume = cubit.get_surface_area(shape)
    #rel_volume = cubit.get_surface_area(shape) / bb[2] / bb[5]
    #rel_meshed_volume = cubit.get_meshed_volume_or_area("surface",[shape]) / bb[2] / bb[5]

  else:
    height = 10
    # mesh shape
    cubit.silent_cmd('volume {} size {} '.format(shape, meshsize))
  
    pairs = get_pairs(shape)
  
    for source_surf, target_surf, source_curf, target_curf, source_vert, target_vert in pairs:
      cubit.silent_cmd('surface {:d} scheme trimesh'.format(source_surf))
      cubit.silent_cmd('mesh surface {:d}'.format(source_surf))
      cubit.silent_cmd('copy mesh surface {:d}  onto  surface {:d}  source curve {:d}  source vertex {:d}  target curve {:d}  target vertex {:d}   smooth '.format(
        source_surf, target_surf, source_curf, source_vert, target_curf, target_vert))
    cubit.silent_cmd('volume all scheme tetmesh ')
    cubit.silent_cmd('mesh volume {}'.format(shape))
  
    nelem = cubit.get_tet_count()
    abs_volume = cubit.get_total_volume([shape])
  
  rel_volume = abs_volume / hull.volume
  rel_meshed_volume = cubit.get_meshed_volume_or_area(entity, [shape]) / hull.volume
  
  cubit.cmd('export ansys cdb "{}.cdb"  overwrite '.format(filename))

  end = time.time()
  print('Time for meshing: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))))

  print('Absolute Surface Area/Volume: {}'.format(abs_volume))
  print('Relative Surface Area/Volume: {}'.format(rel_volume))
  print('Relative Meshed Area/Volume: {} ({} elements)\n'.format(rel_meshed_volume, nelem))

  with open('{}.log'.format(filename), 'a+') as fid:
    print(filename, file=fid)
    print('Absolute Surface Area/Volume: {}'.format(abs_volume), file=fid)
    print('Relative Surface Area/Volume: {}'.format(rel_volume), file=fid)
    print('Relative Meshed Area/Volume: {} ({} elements)\n'.format(rel_meshed_volume, nelem), file=fid)
    print('Time for meshing: {}'.format(time.strftime('%H h %M m %S s', time.gmtime(end-start))), file=fid)

    print('Number of nodes in lefty:  {}'.format(cubit.get_nodeset_node_count(1)), file=fid)
    print('Number of nodes in bottomy: {}'.format(cubit.get_nodeset_node_count(2)), file=fid)
    print('Number of nodes in righty:  {}'.format(cubit.get_nodeset_node_count(3)), file=fid)
    print('Number of nodes in topy: {}'.format(cubit.get_nodeset_node_count(4)), file=fid)
    print('Number of nodes in backy:  {}'.format(cubit.get_nodeset_node_count(5)), file=fid)
    print('Number of nodes in fronty: {}'.format(cubit.get_nodeset_node_count(6)), file=fid)


def mesh_shape_with_triangle(shape, meshsize, filename, holes=None, periodic=False):
#  cubit.silent_cmd('curve all size {}'.format(meshsize))
#  cubit.silent_cmd('curve all scheme equal')
#  cubit.silent_cmd('mesh curve all')

  def split_curve(curve, nsegs):
    cc = cubit.curve(curve)
    has_curvature = any( np.abs( cc.curvature(cc.position_from_fraction(0.5))) > 1e-10)
    if not has_curvature:
      for i in range(1, nsegs+1):
        cubit.silent_cmd("create vertex on curve {} fraction {} from start".format(curve, i/(nsegs+1)))
    else:
      old_vertices = cubit.get_entities("vertex")
      for i in range(1, nsegs):
        cubit.silent_cmd("create vertex on curve {} fraction {} from start".format(curve, i/nsegs))
      new_vertices = cubit.get_entities("vertex")
      curves_to_delete.append(curve)

      vertices = cubit.get_relatives("curve", curve, "vertex")
      new_vertices = [vertex for vertex in new_vertices if vertex not in old_vertices]
      new_vertices.insert(0, vertices[0])
      new_vertices.append(vertices[1])
      for vert1, vert2 in zip(new_vertices, new_vertices[1:]):
        cubit.silent_cmd("create curve vertex {} {}".format(vert1, vert2))

  all_vertices = cubit.get_entities("vertex")
  all_curves = cubit.get_entities("curve")

  nodeset_id_list = cubit.get_nodeset_id_list()
  boundary_curve = {}
  for ns in nodeset_id_list:
    boundary_curve[ns] = cubit.get_nodeset_curves(ns)

  curves_to_delete = []
  for curve in all_curves:
    cc = cubit.curve(curve)
    nsegs = int(cc.length() / meshsize)
    split_curve(curve, nsegs)
    # has_curvature = any( np.abs( cc.curvature(cc.position_from_fraction(0.5))) > 1e-10)
    # if has_curvature:
    #   nsegs = int(cc.length() / meshsize)
    #   if nsegs > 1:
    #     #cubit.silent_cmd("split curve {} segment {} merge".format(curve, nsegs))
    #     split_curve(curve, nsegs)
    #     curves_to_delete.append(curve)
    # else:
    #   for bc in boundary_curve.values():
    #     if curve in bc:
    #       nsegs = int(cc.length() / meshsize)
    #       if nsegs > 1:
    #         #pass
    #         split_curve(curve, nsegs)

  all_vertices = cubit.get_entities("vertex")
  all_curves = cubit.get_entities("curve")
  all_curves = [curve for curve in all_curves if curve not in curves_to_delete]
  print(all_curves)

  # vertex_data contains the original vertex index, the vertex coordinates and a marker
  vertex_data = np.zeros((len(all_vertices), 4))
  for i, vertex in enumerate(all_vertices):
    vertex_data[i, 0] = vertex
    vertex_data[i, 1:3] = cubit.vertex(vertex).coordinates()[0:2]

  # curve_data contains the original curve index, the row indices of the associated vertices of a curve and a marker
  curve_data = np.zeros((len(all_curves), 4))
  for i, curve in enumerate(all_curves):
    curve_data[i, 0] = curve
    vertices = cubit.get_relatives("curve", curve, "vertex")
    for j, vertex in enumerate(vertices):
      curve_data[i, 1+j] = np.where(vertex_data[:,0] == vertex)[0][0]

  # set markers
  nodeset_id_list = cubit.get_nodeset_id_list()
  for ns in nodeset_id_list:
    vertex_list = cubit.get_nodeset_vertices(ns)
    idx = [x in vertex_list for x in vertex_data[:,0]]
    if vertex_list:
      key = list(possible_nodesets.keys())[list(possible_nodesets.values()).index(cubit.get_exodus_entity_name("nodeset", ns))]
      vertex_data[idx, 3] = key

    curve_list = cubit.get_nodeset_curves(ns)
    idx = [x in curve_list for x in curve_data[:,0]]
    if curve_list:
      key = list(possible_nodesets.keys())[list(possible_nodesets.values()).index(cubit.get_exodus_entity_name("nodeset", ns))]
      curve_data[idx, 3] = key

  # create poly file
  with open('{}.poly'.format(filename), 'w+') as fid:
    print('# generated by matviz_cubit.py', file = fid)

  # write poly file. append is necessary for numpy.savetxt
  with open('{}.poly'.format(filename), 'a+') as fid:
    # vertices
    vertices = np.c_[np.arange(0, len(all_vertices)), vertex_data[:,1:]]
    print('{:d} 2 0 1'.format(len(all_vertices)), file=fid)
    np.savetxt(fid, vertices, fmt='%d %f %f %d')
    
    # segments
    segments = np.c_[np.arange(0, len(all_curves)), curve_data[:,1:]]
    print('{:d} 1'.format(len(all_curves)), file=fid)
    np.savetxt(fid, segments, fmt='%d %f %f %d')

    # holes
    if holes is not None:
      print('{:d}'.format(len(holes)), file=fid)
      holes = np.c_[np.arange(0, len(holes)), holes]
      np.savetxt(fid, holes, fmt='%d %f %f')

  print('\nMeshing with triangle...')

  cmd = ["/home/daniel/code/triangle/triangle", "-q30 -zja{:f}".format(pow(meshsize,2)), "-p", "{}.poly".format(filename)]
  if periodic:
    # -Y  No new vertices on the boundary.
    cmd.append("-s")
  
  print(cmd)
  out = subprocess.run(cmd, capture_output=True)
  print(out.stdout.decode("utf-8"))


# this is around 50-100 times faster than the cubit mesher
def mesh_shape_with_gmsh(shape, meshsize, filename):
  import gmsh

  cubit.silent_cmd('export step "{}.stp" overwrite'.format(filename))

  # Cubit 2022.4 crashes for gmsh.initialize()
  # thus we write the code to a file and call the system python
  with open('gmsh_script.py', 'w+') as fid:
    print("""
#!python3
import gmsh
import os
import numpy as np

if not gmsh.is_initialized():
  gmsh.initialize()
gmsh.clear()

# import shape
gmsh.merge('{}' + '.stp')
gmsh.model.geo.synchronize()

xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1,-1)

eps = 1e-6
bottom = gmsh.model.getEntitiesInBoundingBox(xmin-eps, ymin-eps, zmin-eps, xmax+eps, ymin+eps, zmin+eps, 1)
top = gmsh.model.getEntitiesInBoundingBox(xmin-eps, ymax-eps, zmin-eps, xmax+eps, ymax+eps, zmin+eps, 1)
left = gmsh.model.getEntitiesInBoundingBox(xmin-eps, ymin-eps, zmin-eps, xmin+eps, ymax+eps, zmin+eps, 1)
right = gmsh.model.getEntitiesInBoundingBox(xmax-eps, ymin-eps, zmin-eps, xmax+eps, ymax+eps, zmin+eps, 1)

if not left:
  for curve in gmsh.model.getEntities(1):
    lxmin, lymin, lzmin, lxmax, lymax, lzmax = gmsh.model.getBoundingBox(*curve)
    if (abs(lymin - ymin) < eps and abs(lymax - ymax) < eps) and abs(lxmin - xmin) < eps:
      left.append(curve)
    if (abs(lymin - ymin) < eps and abs(lymax - ymax) < eps) and abs(lxmax - xmax) < eps:
      right.append(curve)

# name regions and nodes
gmsh.model.addPhysicalGroup(1, [curve[1] for curve in left], name="lefty")
gmsh.model.addPhysicalGroup(1, [curve[1] for curve in bottom], name="bottomy")
gmsh.model.addPhysicalGroup(1, [curve[1] for curve in right], name="righty")
gmsh.model.addPhysicalGroup(1, [curve[1] for curve in top], name="topy")
gmsh.model.addPhysicalGroup(2, [1], name="mech")

# set periodicity of mesh
translation = [1, 0, 0, 1,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]
for i in left:
  # get the bounding box of each left surface
  lxmin, lymin, lzmin, lxmax, lymax, lzmax = gmsh.model.getBoundingBox(i[0], i[1])
  # For all the matches, we compare the corresponding bounding boxes...
  for j in right:
    lxmin2, lymin2, lzmin2, lxmax2, lymax2, lzmax2 = gmsh.model.getBoundingBox(j[0], j[1])
    lxmin2 -= 1
    lxmax2 -= 1
    # ...and if they match, we apply the periodicity constraint
    if (abs(lxmin2 - lxmin) < eps and abs(lxmax2 - lxmax) < eps
    and abs(lymin2 - lymin) < eps and abs(lymax2 - lymax) < eps
    and abs(lzmin2 - lzmin) < eps and abs(lzmax2 - lzmax) < eps):
      gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], translation)

translation = [1, 0, 0, 0.5,  0, 1, 0, np.sqrt(3)/2,  0, 0, 1, 0,  0, 0, 0, 1]
for i in bottom:
  # get the bounding box of each left surface
  lxmin, lymin, lzmin, lxmax, lymax, lzmax = gmsh.model.getBoundingBox(i[0], i[1])
  # For all the matches, we compare the corresponding bounding boxes...
  for j in top:
    lxmin2, lymin2, lzmin2, lxmax2, lymax2, lzmax2 = gmsh.model.getBoundingBox(j[0], j[1])
    lxmin2 -= 0.5
    lxmax2 -= 0.5
    lymin2 -= ymax
    lymax2 -= ymax
    # ...and if they match, we apply the periodicity constraint
    if (abs(lxmin2 - lxmin) < eps and abs(lxmax2 - lxmax) < eps
    and abs(lymin2 - lymin) < eps and abs(lymax2 - lymax) < eps
    and abs(lzmin2 - lzmin) < eps and abs(lzmax2 - lzmax) < eps):
      gmsh.model.mesh.setPeriodic(1, [j[1]], [i[1]], translation)

# mesh shape
gmsh.option.setNumber("Mesh.MeshSizeFactor", {});
gmsh.model.mesh.generate(2)

gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write('{}' + ".msh")
  """.format(filename, meshsize, filename), file=fid)
    
  with open('gmsh_script.py', 'a') as fid:
    print("""
print('Number of nodes: {:d}'.format(len(gmsh.model.mesh.getNodes(-1, -1)[0])))

#gmsh.fltk.run()
gmsh.finalize()
    """, file=fid)

  # run the script
  out = subprocess.run([sys.executable, "gmsh_script.py"], capture_output=True)
  assert(out.returncode == 0)
  print(out.stdout.decode("utf-8"))

  # cleanup
  os.remove('{}.stp'.format(filename))
  os.remove('step_export.log')
  os.remove('gmsh_script.py')


def convert_poly_to_cdb(meshfilename):

  node_data = np.loadtxt(meshfilename + '.1.node', skiprows=1)
  elem_data = np.loadtxt(meshfilename + '.1.ele', skiprows=1)
  with open(meshfilename + '.1.poly') as f:
    f.readline()
    second_line = f.readline()
  poly_data = np.loadtxt(meshfilename + '.1.poly', skiprows=2, max_rows=int(second_line.split()[0]))

  node_data[:,0] += 1
  elem_data += 1
  poly_data[:,:3] += 1

  nnode = node_data.shape[0]
  nelem = elem_data.shape[0]

  node_lists = {}
  #for id in seg_list_ids:
  for id in possible_nodesets.keys():
    node_lists[possible_nodesets[id]] = np.unique(poly_data[poly_data[:,-1] == int(id), 1:-1])
    if len(node_lists[possible_nodesets[id]]) == 0:
      node_lists[possible_nodesets[id]] = node_data[node_data[:,-1] == id, 0]

  header = """/COM,ANSYS RELEASE 15.0
! Generated by matviz_cubit
/PREP7
/TITLE,
*IF,_CDRDOFF,EQ,1,THEN
_CDRDOFF= 
*ELSE
NUMOFF,NODE, {:>8d}
NUMOFF,ELEM, {:>8d}
NUMOFF,TYPE,        1
*ENDIF
DOF,DELETE
ET,       1,182"""

  outfile = meshfilename + '.cdb'

  with open(outfile, 'w+') as fid:
    print(header.format(nnode, nelem), file=fid)

  node_data_ = np.c_[node_data[:,0], np.zeros((nnode,2)), node_data[:,1:-1], np.zeros((nnode,1))]
  with open(outfile, 'a+') as fid:
    print('NBLOCK,6,SOLID, {:>9d}, {:>9d}'.format(nnode, nnode), file=fid)
    print('(3i9,6e20.13)', file=fid)
    np.savetxt(fid, node_data_, fmt='%9d%9d%9d%20.13e%20.13e%20.13e')

  elem_data_ = np.c_[np.ones((nelem,4)), np.zeros((nelem,4)), 4*np.ones((nelem,1)), np.zeros((nelem,1)), elem_data, elem_data[:,-1]]
  with open(outfile, 'a+') as fid:
    print('N,R5.3,LOC, -1\nEBLOCK, 19, SOLID', file=fid)
    print('(19i9)', file=fid)
    np.savetxt(fid, elem_data_, fmt=''.join('%9d' for _ in range(15)))
    print('       -1', file=fid)

    print('CMBLOCK,mech    ,ELEM,{:>8d}'.format(nelem), file=fid)
    print('(8i10)', file=fid)
    nlines = int(nelem/8)
    extra = nelem % 8
    extra_data = elem_data[nlines*8:,0].astype(int)
    np.savetxt(fid, np.reshape(elem_data[:nlines*8,0], (nlines,8)), fmt=''.join('%10d' for _ in range(8)))
    print(''.join('{:>10d}' for _ in range(extra)).format(*extra_data), file=fid)

    for node_list_name in node_lists:
      nl = node_lists[node_list_name]
      print('CMBLOCK,{:<8s},NODE,{:>8d}'.format(node_list_name, len(nl)), file=fid)
      print('(8i10)', file=fid)
      nlines = int(len(nl)/8)
      extra = len(nl) % 8
      np.savetxt(fid, np.reshape(nl[:nlines*8], (nlines,8)), fmt=''.join('%10d' for _ in range(8)))
      print(''.join('{:>10d}' for _ in range(extra)).format(*nl[nlines*8:].astype(int)), file=fid)


if __name__ == "__main__":
  print('Hello!')
