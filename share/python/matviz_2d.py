from matviz_rot import *
import platform
from PIL import Image, ImageDraw, ImageFont, ImageOps
import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
# matplotlib.use('TkAgg')
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
from scipy import ndimage
import numpy as np
import scipy.interpolate as ip
import math


# create and prepare a matplot figure where patched might be "added" to
def create_figure(min, max, res, for_save):
  # we set the aspect ratio and also the resolution such that we can export as png
  # the problem is that we set the size of the figure but export the subplot w/o axes which is smaller than the figure

  # the dirty solution is to create the figure twice scaled by the error
  dpi_x = (res / 100) * (max[0] - min[0])
  dpi_y = dpi_x * (max[1] - min[1]) / (max[0] - min[0])

  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_y))
  ax = fig.add_subplot(111)
  if for_save:
    # we need to correct the ratio
    wrong = ax.get_window_extent().size
    ratio = dpi_x / dpi_y
    dpi_x *= res / wrong[0]
    dpi_y *= (dpi_y * 100 / ratio) / wrong[1]
    fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_y))
    matplotlib.pyplot.axis('off')
    ax = fig.add_subplot(111)
    # the second figure would make problems with matplotlib.pyplot.show()

  ax.set_xlim(min[0], max[0])
  ax.set_ylim(min[1], max[1])
  return fig, ax

# @param min/max minimal/maximal real node (not barycenter)
# @return image, draw, dim of image, dx to scale from node coord to image coords
def create_image(min, max, res, color="white"):
  if np.size(res) == 1:
    dim = (res, int(res * (max[1] - min[1]) / (max[0] - min[0])))
  else:
    dim = res

  dx = dim[0] / (max[0] - min[0])
  dy = dim[1] / (max[1] - min[1])

  # print "dx=" + str(dx) + " dy=" + str(dy) + " res=" + str(res) + " dim=" + str(dim) + " max=" + str(max)
  im = Image.new("RGB", dim, color)
  draw = ImageDraw.Draw(im)

  return im, draw, dim, dx, dy


# @return phi, r
def to_polar(x, y):
  return numpy.sqrt(x ** 2 + y ** 2), numpy.arctan2(y, x)

# polar coordiantes to cartesian
def to_cart(phi, r):
  return r * numpy.cos(phi), r * numpy.sin(-phi)

# calculate volume for s1, s2. Assume regular grid!
def calc_volume(s1, s2):
  vol = 0.0
  for i in range(len(s1)):
    vol += s1[i] + s2[i] - s1[i] * s2[i]
  if vol.ndim > 1:
    vol = vol[0]
  return (vol / len(s1))

def color_code(color_map, value):
  c = color_map.to_rgba(value)
  return "rgb(" + str(int(255 * c[0])) + ", " + str(int(255 * c[1])) + "," + str(int(255 * c[2])) + ")"


# generate polygon vertices out of rotation data
# to be applied as draw.polygon(result, fill="green", outline="black")
#
# from paraview_fmo import *
# c0 = [12.6, 11.7, 2.3, 0, 0, 8.41]
# from PIL import Image, ImageDraw
# im = Image.new("RGB",(200,200), "white")
# draw = ImageDraw.Draw(im)
# t = to_polygons(data, 100, 100, 2)
# draw.polygon(t, fill="green", outline="black")
# im.show()
# @param only_pos if True only positive values are drawn, with False only the negative (to be used with another color)
def to_polygons(angle, data, x_offset, y_offset, scale, only_pos):
  tupl = []
  for i in range(len(data)):
    r = data[i]
    if only_pos:
      r = r if r > 0 else 0
    else:
      r = numpy.abs(r) if r < 0 else 0
    r = numpy.abs(r)
    x = r * numpy.cos(angle[i])
    y = r * numpy.sin(angle[i])
    tupl.append((x_offset + scale * x, y_offset + scale * y))
  return tupl

# give the corners to draw a rotated rectangle as polygon
def to_rectangle_center(height, width, center, angle = 0.0):

  # print "h=" + str(height) + " w=" + str(width) + " a=" + str(angle) + " x=" + str(x_offset) + " y=" + str(y_offset)

  # height = numpy.max((height,1))

  tupl = []

  for x in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]:
    p = (x[0] * width / 2, x[1] * height / 2)
    r = (cos(angle) * p[0] - sin(angle) * p[1], sin(angle) * p[0] + cos(angle) * p[1])
    # r = (cos(angle) * p[0] + sin(angle)*p[1], - sin(angle) * p[0] + cos(angle) * p[1])
    tupl.append((center[0] + r[0], center[0] + r[1]))

  return tupl

# draws a rotated frustum
# @param direction either vertical or horizontal, not all! or 'sagittal'
def to_frustum_center(start, end, center, elem, scale, direction):
  # 4 ------- 3
  # |         |
  # |         |
  # 1---------2


  angle = 0.5 * (start[2] + end[2])
    # print "horizontal line
  if direction == 'vertical':
    angle += 0.5 * numpy.pi

  idx = 0 if direction == 'vertical' else 0
  alt_idx = 1 if direction == 'vertical' else 1

  # print "to_frustum_centet: start=" + str(start) + " end=" + str(end) + " center=" + str(center) + " d=" + str(direction) + ' idx=' + str(idx)

  val_1 = start[idx]
  val_2 = end[idx]

  tupl = []
  points = []

  # print start
  # print end
  # print elem

  #WARNING CHECK the scaling!!!

  # we scale the element scale, such it overlaps. therefore we need to scale reciproc with the vale, which cancels!
  points.append((-1.0 * scale * elem[alt_idx] / 2, 1. / max((1, scale)) * -val_1 * scale * elem[idx] / 2))
  points.append((1.0 * scale * elem[alt_idx] / 2, 1. / max((1, scale)) * -val_2 * scale * elem[idx] / 2))
  points.append((1.0 * scale * elem[alt_idx] / 2, 1. / max((1, scale)) * val_2 * scale * elem[idx] / 2))
  points.append((-1.0 * scale * elem[alt_idx] / 2, 1. / max((1, scale)) * val_1 * scale * elem[idx] / 2))

  for i in range(4):
    # print "i=" + str(i + 1) + " -> " + str(points[i])
    r = (cos(angle) * points[i][0] - sin(angle) * points[i][1], sin(angle) * points[i][0] + cos(angle) * points[i][1])
    tupl.append(((center[0] + r[0]), (center[1] + r[1])))

  return tupl

# give the corners to draw a rotated rectangle as polygon
def to_rectangle_corner(lower, upper):

  # print "to_rectangle_corner " + str(lower) + " -> " + str(upper)

  tupl = []

  tupl.append(lower)
  tupl.append((upper[0], lower[1]))
  tupl.append(upper)
  tupl.append((lower[0], upper[1]))

  return tupl

# helper
def get_interpol_data(coords, data, fallback, x, eval=True):
  if not eval:
    return None, None
  v = data[x]
  if v[0] == -1.0:
    v = fallback[x]
  return coords[x], v

# helper for get_interpolation and Fields
# @return 2d locations and 2d data
def convert_two_data_interpolation_input(centers, s1, s2, angle):
  # convert to 2D
  c = numpy.zeros((len(centers), 2))
  c[:, 0] = [x[0] for x in centers]
  c[:, 1] = [x[1] for x in centers]

  v = numpy.zeros((len(s1),  2 if angle is None else 3))
  for i in range(len(s1)):
    v[i][0] = s1[i][0]
    v[i][1] = s2[i][0]
    if angle is not None:
      v[i][2] = angle[i]

  return c, v



# helper for get_interpolation and Fields
# @return 2d locations and 2d data
def convert_single_data_interpolation_input(centers, s1, angle):
  # convert to 2D
  c = numpy.zeros((len(centers), 2))
  c[:, 0] = [x[0] for x in centers]
  c[:, 1] = [x[1] for x in centers]

  v = numpy.zeros((len(s1),  1 if angle is None else 2))
  for i in range(len(s1)):
    v[i][0] = s1[i][0]
    if angle is not None:
      v[i][1] = angle[i]

  return c, v


# helper which returns an interpolated grid and one nearest neighbor interpolated grid
def get_interpolation(coords, grad, sample, s1, s2, angle=None):
  assert(sample == 'elem_nodes' or sample == 'edge_centers' or sample == 'elem_centers')

  centers, min, max, elem = coords
  # where we want nodes
  nx = int((max[0] - min[0]) / elem[0])
  ny = int((max[1] - min[1]) / elem[1])
  out = None
  if sample == 'elem_nodes':
    out = numpy.zeros(((nx + 1) * (ny + 1), 2))
    for y in range(ny + 1):
      for x in range(nx + 1):
        out[y * (nx + 1) + x][0] = float(x) / nx * max[0]
        out[y * (nx + 1) + x][1] = float(y) / ny * max[1]
  elif sample == 'edge_centers':  # 'edge_centers' this is much more complicated, draw a coarse grid and mark all element edge centers!
    # ---5-------6----
    # |      |       |
    # 2      3       4
    # |      |       |
    # ---0-------1----

    out = numpy.zeros(((2 * nx + 1) * (2 * ny + 1), 2))
    # 0,1,5,6
    for y in range(ny + 1):
      for x in range(nx):
        out[y * (2 * nx + 1) + x][0] = float(x) / nx * max[0] + 0.5 * elem[0]
        out[y * (2 * nx + 1) + x][1] = float(y) / ny * max[1]
        # print "out[" + str(y * (2*nx+1) + x) + "] = " + str(out[y * (2*nx+1) + x])
    # 2.3.4
    for y in range(ny):
      for x in range(nx + 1):
        out[nx + y * (2 * nx + 1) + x][0] = float(x) / nx * max[0]
        out[nx + y * (2 * nx + 1) + x][1] = float(y) / ny * max[1] + 0.5 * elem[1]
        # print "out[" + str(nx + y * (2*nx+1) + x) + "] = " + str(out[nx + y * (2*nx+1) + x])
  else: # sample == 'elem_centers'
    out = numpy.zeros((nx * ny, 2))
    for y in range(ny):
      for x in range(nx):
        out[y * nx + x][0] = (float(x) + 0.5) / nx * max[0]
        out[y * nx + x][1] = (float(y) + 0.5) / ny * max[1]


  if s2 is None:
    c, v = convert_single_data_interpolation_input(centers, s1, angle)
  else:
    c, v = convert_two_data_interpolation_input(centers, s1, s2, angle)

  ip_data = ip.griddata(c, v, out, grad, -1.0)
  # any interpolation but nearest neighbor can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation
  ip_near = ip.griddata(c, v, out, 'nearest') if grad != 'nearest' else None

  return ip_data, ip_near, out, nx, ny

# visualize the orientational stiffness
# @param grad is 'none' or 'linear'
# @return the image
def show_frame_grad(coords, design, grad, direction, nx):

  s1 = design['s1']
  s2 = design['s2']

  centers, min, max, elem = coords
  im, draw, dim, dx, dy = create_image(min, max, nx, "white")

  height = elem[1] * dy
  length = elem[0] * dx

  # print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " height=" + str(height) + " length=" + str(length) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'elem_nodes', s2, s1) # for some strange reason we need to switch?!

  for y in range(ny + 1):
    for x in range(nx + 1):
      start, v_start = get_interpol_data(out, ip_data, ip_near, y * (nx + 1) + x)
      right, v_right = get_interpol_data(out, ip_data, ip_near, y * (nx + 1) + x + 1, eval=True if x < nx else False)
      upper, v_upper = get_interpol_data(out, ip_data, ip_near, (y + 1) * (nx + 1) + x, eval=True if y < ny else False)

      # print "start=" + str(start) + " v_start=" + str(v_start) + " right=" + str(right) + " v_right=" + str(v_right)

      # 4 ------- 3
      # |         |
      # |         |
      # 1---------2

      # print horizontal line
      if not direction == 'vertical':
        if x < nx:  # the right most vertical line has no right for the horizontal line
          n1x = start[0] * dx
          n2x = right[0] * dx
          n1y = dim[1] - start[1] * dy - height * 0.5 * v_start[1]
          n4y = dim[1] - start[1] * dy + height * 0.5 * v_start[1]
          n2y = dim[1] - right[1] * dy - height * 0.5 * v_right[1]
          n3y = dim[1] - right[1] * dy + height * 0.5 * v_right[1]
          if int(n1y) == int(n4y):
            n4y = n1y + 1.01

          if int(n2y) == int(n3y):
            n3y = n2y + 1.01

          tupels = []
          tupels.append((n1x, n1y))
          tupels.append((n2x, n2y))
          tupels.append((n2x, n3y))
          tupels.append((n1x, n4y))

          draw.polygon(tupels, fill="black")

      # print vertical line
      if not direction == 'horizontal':
        if y < ny:  # the evaluation of the most upper horizontal line has no upper for vertical data
          n1y = dim[1] - start[1] * dy
          n4y = dim[1] - upper[1] * dy
          n1x = start[0] * dy - length * 0.5 * v_start[0]
          n2x = start[0] * dy + length * 0.5 * v_start[0]
          n4x = upper[0] * dy - length * 0.5 * v_upper[0]
          n3x = upper[0] * dy + length * 0.5 * v_upper[0]

          if int(n1x) == int(n2x):
            n2x = n1x + 1.01

          if int(n3x) == int(n4x):
            n4x = n3x + 1.01

          tupels = []
          tupels.append((n1x, n1y))
          tupels.append((n2x, n1y))
          tupels.append((n3x, n4y))
          tupels.append((n4x, n4y))

          draw.polygon(tupels, fill="black")

  return im

# visualize the orientational stiffness
# @return the image
def show_rot_cross_grad(coords, design, grad, direction, nx, scale, color, do_save):

  s1 = design['s1']
  s2 = design['s1']
  angle = design['angle']

  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)
  delta_angle = numpy.max(angle[:]) - numpy.max(angle[:])

  if scale == -1:
    scale = -1.0 if delta_angle == 0.0 else -0.8

  # print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'edge_centers', s1, s2, angle)

  # out is rather complex set, see get_interpolation for details!
  if not direction == 'vertical':
    for y in range(ny):
      for x in range(nx):
        start, v_start = get_interpol_data(out, ip_data, ip_near, nx + y * (2 * nx + 1) + x)
        right, v_right = get_interpol_data(out, ip_data, ip_near, nx + y * (2 * nx + 1) + x + 1)

        center = ((0.5 * (start[0] + right[0]), start[1]))
        # print 'start=' + str(start) + ' v_start=' + str(v_start)

        pol = to_frustum_center(v_start, v_right, center, elem, scale, 'horizontal')
        draw_verts(pol, sub, 'black')

  if not direction == 'horizontal':
    for y in range(ny):
      for x in range(nx):
        start, v_start = get_interpol_data(out, ip_data, ip_near, (y + 1) * (2 * nx + 1) + x)
        upper, v_upper = get_interpol_data(out, ip_data, ip_near, y * (2 * nx + 1) + x)
        center = ((start[0], 0.5 * (start[1] + upper[1])))
        pol = to_frustum_center(v_upper, v_start, center, elem, scale, 'vertical')
        draw_verts(pol, sub, 'black')

  return (fig, sub)

def show_modified_frame(coords, design, nx):
  #TODO: implement scale, direction, color, do_save
  s1 = design['s1']
  s2 = design['s1']
  angle = design['angle']

  centers, min, max, elem = coords
  im, draw, dim, dx, dy = create_image(min, max, nx, "black")
  height = elem[1] * dy
  length = elem[0] * dx
  for i in range(len(s1)):
    coord = centers[i]
    x_off = (coord[0] + min[0] - 0.5 * elem[0]) * dim[0]
    y_off = (coord[1] + min[1] - 0.5 * elem[1]) * dim[1]
    ver = s2[i, 0]
    hor = s1[i, 0]

    pix = im.load()
    eps = 1e-8
    offx = int((length / 2.) * (ver) + 0.5+eps)
    offy = int((height / 2.) * (hor) + 0.5+eps)
    for i in range(offx, int(length+eps) - offx):
       for j in range(offy, int(height+eps) - offy):
          pix[int(x_off+i+eps),int(y_off+j+eps)] = (255,255,255)
      # 2D frame structure
    # modify frame for stress minimization
    for i in range(offx, int(length - offx+eps)):
      for j in range(offy, int(height - offy+eps)):
        if math.ceil((int(length+eps) - 2.*offx) / 3.) <= math.ceil((int(height+eps) - 2.*offy) / 3.):
          r = math.ceil((int(length+eps) - 2.*offx) / 3.)
        else:
          r = math.ceil((int(height+eps) - 2.*offy) / 3.)
        m = [offx + r, offy + r]
        if i - offx < r and j - offy < r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
          pix[int(x_off+i+eps),int(y_off+j+eps)] = (0,0,0)
        m = [int(length+eps) - offx - r - 1, offy + r]
        if i >= int(length+eps) - offx - r and j - offy < r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
          pix[int(x_off+i+eps),int(y_off+j+eps)] = (0,0,0)
        m = [int(length+eps) - offx - r - 1, int(height+eps) - offy - r - 1]
        if i >= int(length+eps) - offx - r and j >= int(height+eps) - offy - r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
          pix[int(x_off+i+eps),int(y_off+j+eps)] = (0,0,0)
        m = [offx + r, int(height+eps) - offy - r - 1]
        if i - offx < r and j >= int(height+eps) - offy - r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
          pix[int(x_off+i+eps),int(y_off+j+eps)] = (0,0,0)
  return im

# visualizes the oriental stiffness as frame with smooth inner corners; creates a vector image
def show_modified_frame_old(coords, design, direction, nx, scale, color, do_save):
  print('image is only correct if eps<= s1,s2 <= 0.5, otherwise scaling is necessary; rotation is not implemented currently')

  s1 = design['s1']
  s2 = design['s1']
  angle = design['angle']

  s1 /= 2.
  s2 /= 2.
  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)
  delta_angle = numpy.max(angle) - numpy.max(angle)

  if scale == -1.0:
    scale = 1.02 if delta_angle == 0.0 else 0.8
  length = scale * (elem[0])

  max_val = numpy.max([numpy.max(s1), numpy.max(s2)])
  min_val = numpy.min([numpy.min(s1), numpy.min(s2)])

  for i in range(len(s1)):

    coord = centers[i]

    x_off = (coord[0] + min[0])
    y_off = (coord[1] + min[1])

    # we need downscale the values when we overscale due to overlapping
    v = [0, 0]
    v[0] = s1[i, 0] / numpy.max((scale, 1.))
    v[1] = s2[i, 0] / numpy.max((scale, 1.))
    theta = angle[i]
    c = [0, 0]
    c[0] = str(1.0 - v[0] / max_val) if color == "grayscale" else 'black'
    c[1] = str(1.0 - v[1] / max_val) if color == "grayscale" else 'black'

    # vmax, vmin are indices which decide whether s1-rectangle or s2-rectangle are drawn first
    vmax,vmin = (0,1) if v[0] > v[1] else (1,0)
    vmin = (vmax + 1) % 2
    pol = to_rectangle_center(length, length, (x_off, y_off), theta + vmin * numpy.pi / 2)
    draw_verts(pol, sub, 'black')

    pol = to_rectangle_center(length * (1. - 2.*v[0]), length * (1. - 2.*v[1]), (x_off, y_off), theta + vmax * numpy.pi / 2)
    draw_verts(pol, sub, 'white')

    r = (1. - 2.*v[0]) * length / 3. if (1. - 2.*v[0]) * length / 3. >= (1. - 2.*v[1]) * length / 3. else (1. - 2.*v[1]) * length / 3.
    l1 = length * (1. - 2.*v[0])
    l2 = length * (1. - 2.*v[1])
    t = [abs(l1 / 2. - r / 2.), abs(l2 / 2. - r / 2.)]
    c1 = abs(l1 / 2. - r)
    c2 = abs(l2 / 2. - r)
    eps = 1./nx

    # TODO: rotation not implemented yet
    theta = 0.
    vmax = 0.
    vmin = 0.
    # lower left corner
    #t = (cos(theta + vmax * numpy.pi / 2) * -abs(l1 / 2. - r / 2.) - sin(theta + vmax * numpy.pi / 2) *-abs(l2 / 2. - r / 2.), sin(theta + vmax * numpy.pi / 2) * -abs(l1 / 2. - r / 2.) + cos(theta + vmax * numpy.pi / 2) * -abs(l2 / 2. - r / 2.))
    pol = to_rectangle_center(r, r, (x_off-t[0], y_off-t[1]), theta + vmax * numpy.pi/2)
    #pol[:][0] -= t[0]
    #pol[:][1] -= t[1]
    draw_verts(pol, sub, 'black')
    c = (x_off - c1, y_off - c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')

    #t = (cos(theta + vmax * numpy.pi / 2) * abs(l1 / 2. - r / 2.) - sin(theta + vmax * numpy.pi / 2) *abs(l2 / 2. - r / 2.), sin(theta + vmax * numpy.pi / 2) * abs(l1 / 2. - r / 2.) + cos(theta + vmax * numpy.pi / 2) * abs(l2 / 2. - r / 2.))
    # lower right corner
    pol = to_rectangle_center(r, r, (x_off + t[0], y_off - t[1]), 0.0)
    draw_verts(pol, sub, 'black')
    c = (x_off + c1, y_off - c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')

    # upper right corner
    pol = to_rectangle_center(r, r, (x_off + t[0], y_off + t[1]), theta + vmax * numpy.pi / 2)
    draw_verts(pol, sub, 'black')
    c = (x_off + c1, y_off + c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')

    # upper left corner
    pol = to_rectangle_center(r, r, (x_off - t[0], y_off + t[1]), theta + vmax * numpy.pi / 2)
    draw_verts(pol, sub, 'black')
    c = (x_off - c1, y_off + c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')
  return (fig, sub)

# visualize the orientational stiffness
# @param grad is 'none' or 'linear'
# @return the image
def show_frame(coords, design, nx, scale):

  s1 = design['s1']
  s2 = design['s2']

  if scale == -1.:
    scale = 1.
  assert(scale <= 1.)
  centers, min, max, elem = coords
  im, draw, dim, dx, dy = create_image(min, max, nx, "white")
  height = scale*elem[1] * dy
  length = scale*elem[0] * dx
  #print 'elem0 = ' +str(elem[0]) + ' dx = ' +str(dx)
  #print 'height = ' +str(height)+', length= '+str(length)
  #print 'dim = '+str(dim)
  eps = 1e-8
  warn = False
  for i in range(len(s1)):
    coord = centers[i]
    #print 'coord = '+str(coord)
    x_off = int(coord[0]*dx+min[0]-height/(2.)+eps)
    y_off = int(dy-coord[1]*dy+min[1]-height/(2.)+eps)

    #x_off = (coord[0] + min[0] - 0.5 * elem[0]) * dim[0]
    #y_off = nx-1-(length/scale)/2.-(coord[1] + min[1] - 0.5 * elem[1]) * dim[1]
    ver = s2[i, 0]/numpy.max((scale, 1.))
    hor = s1[i, 0]/numpy.max((scale, 1.))
    #print 'hor = '+str(hor)+' ver = '+str(ver)
    #print 'x_off = ' +str(x_off)+', y_off= '+str(y_off)
    pix = im.load()
    for i in range(x_off, int(length+eps) + x_off):
       for j in range(y_off, int(height+eps) + y_off):
          #print('index '+str(i)+ ', index '+str(j))
          pix[i,j] = (0,0,0)
    offx = int((length / 2.) * (ver) + 0.5+eps)
    offy = int((height / 2.) * (hor) + 0.5+eps)
    if offx == 0:
      if offy > 0:
        offx = 1
        warn = True
    elif offy == 0:
      if offx > 0:
        offy = 1
        warn = True
    #print 'offx = ' +str(offx)+', offy= '+str(offy)
    for i in range(offx, int(length+eps) - offx):
       for j in range(offy, int(height+eps) - offy):
          pix[x_off+i,y_off+j] = (255,255,255)
  if warn:
    print('Warning: minimal thickness scaled to 1 Pixel.')
  return im


# @return the image
def show_rot_cross(coords, design, direction, nx, scale, color, do_save):

  s1 = design['s1']
  s2 = design['s1']
  angle = design['angle']

  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)
  delta_angle = numpy.max(angle) - numpy.min(angle)

  if scale == -1.0:
    scale = -1.02 if delta_angle == 0.0 else -0.8

  length = scale * (elem[0])

  max_val = numpy.max([numpy.max(s1), numpy.max(s2)])
  min_val = numpy.min([numpy.min(s1), numpy.min(s2)])
  sm = cmx.ScalarMappable(colors.Normalize(min_val, max_val), cmap=plt.get_cmap('gray' if color == 'grayscale' else color)) if not color == 'black' else 'black'

  for i in range(len(s1)):

    coord = centers[i]

    center = (coord[0] + min[0], coord[1] + min[1])

    # we need to downscale the values when we overscale due to overlapping
    v = [0, 0]
    v[0] = s1[i, 0] / numpy.max((scale, 1.))
    v[1] = s2[i, 0] / numpy.max((scale, 1.))
    theta = angle[i]
    c = [0, 0]
    c[0] = sm.to_rgba(max_val-v[0]) if color != 'black' else 'black'
    c[1] = sm.to_rgba(max_val-v[1]) if color != 'black' else 'black'


  # c[0] = str(1.0 - v[0] / max_val) if color == "grayscale" else 'black'
  # c[1] = str(1.0 - v[1] / max_val) if color == "grayscale" else 'black'

    # print 'S=' + str(s1[i,0]) + '/' + str(s2[i,0])  + ' v=' + str(v) + ' c=' + str(c)
    theta = angle[i]
    # a
    if direction == 'horizontal':
      pol = to_rectangle_center(length * v[0], length, (center[0], dim[1] - center[1]), theta)
      draw_verts(pol, sub, c[0])
    # b
    elif direction == 'vertical':
      pol = to_rectangle_center(length * v[1], length, center, theta + numpy.pi/2.)
      draw_verts(pol, sub, c[1])
    else:
      # vmax, vmin are indices which decide whether s1-rectangle or s2-rectangle are drawn first
      vmax,vmin = (0,1) if v[0] > v[1] else (1,0)
      pol = to_rectangle_center(length * v[vmin], length, center, theta + vmin * numpy.pi/2.)
      # draw_verts(pol, sub, str(1.0 - c[vmin]))
      draw_verts(pol, sub, c[vmin])
      pol = to_rectangle_center(length * v[vmax], length, center, theta + vmax * numpy.pi/2.)
      draw_verts(pol, sub, c[vmax])
  return (fig, sub)

# @return the image
def show_sheared_rot_cross(coords, design, direction, nx, scale, color, do_save):

  s1 = design['s1']
  s2 = design['s2']
  sh1 = design['sh1']
  angle = design['angle']

  centers, min, max, elem = coords

  fig, sub = create_figure(min, max, nx, do_save)

  if scale == -1.0:
    scale = 0.8

  length =  scale * (elem[0])

  max_val = numpy.max([numpy.max(s1), numpy.max(s2)])
  min_val = numpy.min([numpy.min(s1), numpy.min(s2)])

  for i in range(len(s1)):

    coord = centers[i]

    center = (coord[0] + min[0], coord[1] + min[1])

    # we need downscale the values when we overscale due to overlapping
    v = [0,0]
    v[0] = s1[i,0] / numpy.max((scale, 1.))
    v[1] = s2[i,0] / numpy.max((scale, 1.))
    phi = (sh1[i] - .5) * numpy.pi
    c = [0,0]
    c[0] = str(1.0 - v[0] / max_val) if color == "grayscale" else 'black'
    c[1] = str(1.0 - v[1] / max_val) if color == "grayscale" else 'black'

    #print 'S=' + str(s1[i,0]) + '/' + str(s2[i,0])  + ' v=' + str(v) + ' c=' + str(c)
    theta = angle[i]

    if direction == 'horizontal':
      pol = to_rectangle_center(length * v[0], length, (center[0], dim[1] - center[1]), theta)
      draw_verts(pol, sub, c[0])
    elif direction == 'vertical':
      pol = to_rectangle_center(length * v[1], length, center, theta + phi + numpy.pi/2.)
      draw_verts(pol, sub, c[1])
    else:
      vmax = 0 if v[0] > v[1] else 1
      vmin = (vmax + 1) % 2
      shearingangle = [0,0]
      shearingangle[vmax] = phi
      angle = theta + shearingangle[0] + vmin*numpy.pi/2
      pol = to_rectangle_center(length * v[vmin], length, center, angle)
      draw_verts(pol, sub, c[vmin])
      angle = theta + shearingangle[1] + vmax*numpy.pi/2
      pol = to_rectangle_center(length * v[vmax], length, center, angle)
      draw_verts(pol, sub, c[vmax])

  return (fig, sub)

# @return the image
def show_cross_bar(coords, design, direction, nx, scale, color, do_save):

  s1 = design['s1']
  s2 = design['s1']
  s3 = design['s3']
  angle = design['angle']

  centers, min, max, elem = coords

  fig, sub = create_figure(min, max, nx, do_save)

  if scale == -1.0:
    scale = 0.8

  length =  scale * (elem[0])

  max_val = numpy.max([numpy.max(s1), numpy.max(s2), numpy.max(s3)])
  min_val = numpy.min([numpy.min(s1), numpy.min(s2), numpy.min(s3)])

  for i in range(len(s1)):

    coord = centers[i]

    center = (coord[0] + min[0], coord[1] + min[1])

    # we need downscale the values when we overscale due to overlapping
    v = [0,0,0]
    v[0] = s1[i,0] / numpy.max((scale, 1.))
    v[1] = s2[i,0] / numpy.max((scale, 1.))
    v[2] = s3[i,0] / numpy.max((scale, 1.))
    c = [0,0,0]
    c[0] = str(1.0 - v[0] / max_val) if color == "grayscale" else 'black'
    c[1] = str(1.0 - v[1] / max_val) if color == "grayscale" else 'black'
    c[2] = str(1.0 - v[2] / max_val) if color == "grayscale" else 'black'

    if direction == 'vertical':
      pol = to_rectangle_center(length * v[0], length, center, numpy.pi/2.)
      draw_verts(pol, sub, c[0])
    else:
      vidx = numpy.argsort(v)
      vmin = vidx[0]
      vmid = vidx[1]
      vmax = vidx[2]
      angle = [numpy.pi/2, -numpy.pi/4, numpy.pi/4]
      #pol = to_rectangle_center(length * v[vmin], length * sqrt(1+pow(sin(2*angle[vmin]),2)), center, angle[vmin])
      pol = to_rectangle_center(length * v[vmin], length, center, angle[vmin])
      draw_verts(pol, sub, c[vmin])
      #pol = to_rectangle_center(length * v[vmid], length * sqrt(1+pow(sin(2*angle[vmid]),2)), center, angle[vmid])
      pol = to_rectangle_center(length * v[vmid], length, center, angle[vmid])
      draw_verts(pol, sub, c[vmid])
      #pol = to_rectangle_center(length * v[vmax], length * sqrt(1+pow(sin(2*angle[vmax]),2)), center, angle[vmax])
      pol = to_rectangle_center(length * v[vmax], length, center, angle[vmax])
      draw_verts(pol, sub, c[vmax])

  return (fig, sub)

def show_frame2(coords, design, nx, scale, color, do_save):

  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)

  if scale == -1.0:
    scale = -1.02
  length = scale * (elem[0])

  max_val = numpy.max([numpy.max(design[0]), numpy.max(design[1])])
  min_val = numpy.min([numpy.min(design[0]), numpy.min(design[1])])

  if color == 'black':
    sm = 'black'
  else:
    sm = cmx.ScalarMappable(colors.Normalize(min_val, max_val),
        cmap=plt.get_cmap('gray' if color == 'grayscale' else color))

  for i in range(len(design[0])):
    coord = centers[i]

    x_off = (coord[0] + min[0])
    y_off = (coord[1] + min[1])

    # we need downscale the values when we overscale due to overlapping
    v = [0, 0]
    c = [0, 0]
    for t in range(2):
      v[t] = design[t][i, 0] / numpy.max((scale, 1.))
      c[t] = sm.to_rgba(max_val-v[t]) if not color == 'black' else 'black'

    # draw thicker rectangle pair last (on top)
    for t in sorted([0, 1], key=lambda t: v[t]):
      width = v[t] / 2
      corner = [x_off - length / 2, y_off - length / 2]

      for k in range(2):
        if k == 0:
          lower = (0, 0)
          upper = ((1, width) if t == 0 else (width, 1))
        else:
          lower = ((0, 1 - width) if t == 0 else (1 - width, 0))
          upper = (1, 1)

        pol = to_rectangle_corner(lower, upper)
        pol = [(x * length + corner[0],
                y * length + corner[1]) for x, y in pol]
        draw_verts(pol, sub, c[t])

  return (fig, sub)

def show_framed_cross(coords, design, nx, scale, color, do_save):

  design = design['microparams']

  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)

  if scale == -1.0:
    scale = -1.02
  length = scale * (elem[0])

  max_val = numpy.max([numpy.max(design[0]), numpy.max(design[1]),
                       numpy.max(design[2]), numpy.max(design[3])])
  min_val = numpy.min([numpy.min(design[0]), numpy.min(design[1]),
                       numpy.min(design[2]), numpy.min(design[3])])

  if color == 'black':
    sm = 'black'
  else:
    sm = cmx.ScalarMappable(colors.Normalize(min_val, max_val),
        cmap=plt.get_cmap('gray' if color == 'grayscale' else color))

  for i in range(len(design[0])):
    coord = centers[i]

    x_off = (coord[0] + min[0])
    y_off = (coord[1] + min[1])

    # we need downscale the values when we overscale due to overlapping
    v = [0, 0, 0, 0]
    c = [0, 0, 0, 0]
    for t in range(4):
      v[t] = design[t][i, 0] / numpy.max((scale, 1.))
      c[t] = sm.to_rgba(max_val-v[t]) if not color == 'black' else 'black'

    # draw thicker rectangles last (on top)
    for t in sorted([0, 1, 2, 3], key=lambda t: v[t]):
      width = v[t] / 2
      corner = [x_off - length / 2, y_off - length / 2]
      draw_poly = lambda pol: draw_verts([(x * length + corner[0],
                                           y * length + corner[1])
                                          for x, y in pol], sub, c[t])

      if t < 2:
        # horizontal/vertical bars
        for k in range(2):
          if k == 0:
            lower = (0, 0)
            upper = ((1, width) if t == 0 else (width, 1))
          else:
            lower = ((0, 1 - width) if t == 0 else (1 - width, 0))
            upper = (1, 1)

          pol = to_rectangle_corner(lower, upper)
          draw_poly(pol)
      else:
        # diagonal bars
        pol = [(0, 0), (width, 0), (1, 1 - width),
               (1, 1), (1 - width, 1), (0, width)]
        if t == 2: pol = [(x, 1 - y) for x, y in pol]
        draw_poly(pol)

  return (fig, sub)


def show_triangle_grad(coords, design, grad, samples, res, thres, save, access, equilateral=True, radius=None, param=None, repetitions=1):

  # hourglass or diamond
  type = 'diamond'

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
    return int(np.round(start)), int(np.round(end))
  
  def draw_chord(x, y, tupels, idx, point_on_bisec, radius, boundary=False):
    p = tupels[idx]
    p1 = tupels[(idx+1)%3]
    p2 = tupels[(idx+2)%3]

    v1 = np.array(p1) - np.array(p)
    v2 = np.array(p2) - np.array(p)
    assert(np.linalg.norm(v1))
    assert(np.linalg.norm(v2))
    # angle at corner of triangle
    angle = np.arccos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))

    # midpoint of arc
    vec = np.array(point_on_bisec) - np.array(p)
    assert(np.linalg.norm(vec) > 1e-8)
    mid_circ = np.array(p) + vec/np.linalg.norm(vec) * radius / np.sin(angle/2)
    mid_circ = np.round(mid_circ)
    #draw.point([mid_circ[0],mid_circ[1]], fill="yellow") #debug

    # get bounding box of circle
    cp_e = (mid_circ[0] + radius)
    cp_n = (mid_circ[1] + radius)
    cp_w = (mid_circ[0] - radius)
    cp_s = (mid_circ[1] - radius)
    bbox = [cp_w, cp_s, cp_e, cp_n]

    start, end = get_angles_for_chord(x, y, idx, angle, boundary)

    # output for meshing with triangles
    #with open('matviz.out', 'a+') as f:
    #  print('circ {}, {}, {}'.format(bbox, start, end), file=f)
    
    a1 = np.zeros((2,))
    a2 = np.zeros((2,))
    a1[0] = mid_circ[0] + radius * np.cos(start/180*np.pi)
    a1[1] = mid_circ[1] - radius * np.sin(start/180*np.pi)
    a2[0] = mid_circ[0] + radius * np.cos(end/180*np.pi)
    a2[1] = mid_circ[1] - radius * np.sin(end/180*np.pi)
    a1, a2 = np.round(a1), np.round(a2)

    # we will flip the image later upside down, i.e. for plotting we have
    # to go clockwise and switch start and end
    start, end = -end, -start

    # draw section of circle
    draw.chord(bbox, start, end, fill="white")

    return a1, a2

  s1 = design['s1']

  if np.max(s1) < 2.1:
    s1 = np.minimum(s1,1)

  # usually we read physical design, so we have to recalculate the non-penalized design
  # this assumes RAMP with parameter 2.2
  if param is None and access == "smart":
    p = 2.2 # RAMP with p=2.2
    s1 = (p+1)*s1 / (p*s1+1)

  centers, mini, maxi, elem = coords

  # set size dx/dy/dz of one cell
  if samples is not None:
    tmp = samples.split(',')
    samples = [int(tmp[0]), int(tmp[0])] if len(tmp) == 1 else [int(tmp[0]), int(tmp[1])]
      

    dx = (maxi[0] - mini[0]) / samples[0]
    dy = (maxi[1] - mini[1]) / samples[1]
    elem = np.array([dx, dy, 0])

  if equilateral:
    # rescale element for equilateral triangle
    # height = sqrt(3) / 2 * baseline
    # triangles form tesellation, i.e. their surrounding rectangles overlap and thus element width is baseline / 2
    elem = np.array([elem[0], elem[0] * np.sqrt(3), 0])

  nx = int((maxi[0]-mini[0]) / elem[0])
  ny = int(np.round((maxi[1]-mini[1]) / elem[1]))

  if equilateral:
    if param is None:
      # recalculate elem size, such that we always draw complete triangles in y-direction
      elem[1] = (maxi[1]-mini[1]) / ny
    else:
      # if we generate base cells, we draw multiple triangles in x- and y-direction
      # after drawing, we crop the image to obtain a periodic base cell with an "X" and three horizontal bars
      if type == 'hourglass':
        nx = 2*repetitions+1
        ny = repetitions
      elif type == 'diamond':
        nx = 4*repetitions
        ny = repetitions
      else:
        raise Exception('Type not implemented')
      maxi[0] = nx*elem[0]
      maxi[1] = ny*elem[1]

    # calculate new element centers
    centers_new = []
    for y in range(ny):
      for x in range(nx):
        center = np.array([elem[0]/2.0 + x * elem[0], elem[1]/2.0 + y * elem[1]])
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

    # get image resolution
    res = (res, int(res * ny/nx * elem[1] / elem[0])) if param is None else [res, int(res * ny/nx * np.sqrt(3))]

  im, draw, dim, dx, dy = create_image(mini, maxi, res, "white")

  # output for meshing with triangles
  #with open('matviz.out', 'a+') as f:
  #  print('size {}'.format(dim), file=f)

  # if the parameter is zero, we do not have to draw anything and return a white image
  if param == 0.0:
    # remove half of an element on both sides (times dx for pixels)
    # crop argument = (left, upper, right, lower)
    im = im.crop((maxi[0]/nx*1/2*dx, mini[1]*dy, maxi[0]/nx*(nx-1/2)*dx, maxi[1]*dy))
    return im

  # size of element in pixels
  height = elem[1] * dy
  length = elem[0] * dx

  if radius:
    assert(radius >= 0)
    #assert(radius <= 1)
    radius = np.round(radius * length/np.sqrt(3)/4)

  if equilateral:
    # interpolate original data at new centers
    c, v = convert_single_data_interpolation_input(centers, s1, None)
    ip_data = ip.griddata(c, v, centers_new, grad, -1.0)
    # any interpolation but nearest neighbor can only interpolate in the convex hull
    ip_near = ip.griddata(c, v, centers_new, 'nearest') if grad != 'nearest' else None
  else:
    ip_near = s1
    ip_data = s1

  for y in range(ny):
    for x in range(nx):
      rad = radius

      # if the value is -1 we use the nearest interpolation
      # mid IS NOT THE MIDPOINT OF THE TRIANGLE BUT OF elem (=rectangle)
      mid, v = get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+0)
      _, v1 = get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+1)
      _, v2 = get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+2)
      _, v3 = get_interpol_data(centers_new if equilateral else centers, ip_data, ip_near, (y * nx + x)*4+3)

      # kind of integration -> should be removed and graded triangle bars should be used
      v = (3*v + v1+v2+v3)/6

      if v < thres:
        v = 0.0
      #if v > thres:
      #  v = 1.0

      # ensure param is in bounds
      if param is not None:
        v = np.min((1.0, np.max((0.0, param))))

      if isinstance(v,np.ndarray):
        v = v[0]

      if v < -1e-12 or v > 1+1e-12:
        raise ValueError("argument v={:.3f} out of bounds".format(v))

      v = np.min((1.0, np.max((0.0, v))))
      # this will generate a linear param-volume realtionship
      # else it will be volume = -param^2 + 2*param
      v = 1-np.sqrt(1-v)

      # midpoint of equilateral triangle is h/3. Thus v/2 should be in [0,1/3].
      v *= 2.0/3.0

      # black outer triangle
      # This would not be necessary if im was created as black image.
      # But then, param == 0.0, would have to be treated in another way above.
      n1x = mid[0] * dx - length
      n2x = mid[0] * dx
      n3x = mid[0] * dx + length
      n4x = mid[0] * dx - length * 0.5
      n5x = mid[0] * dx + length * 0.5
      if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
        #    2        4--2           2--5
        #   / \       |   \         /   |
        #  /   \      |    \       /    |
        # 1-----3     1-----3     1-----3
        n1y = mid[1] * dy - height * 0.5
        n2y = mid[1] * dy + height * 0.5
      else:
        # upside down
        # 1-----3     1-----3     1-----3
        #  \   /      |    /       \    |
        #   \ /       |   /         \   |
        #    2        4--2           2--5
        n1y = mid[1] * dy + height * 0.5
        n2y = mid[1] * dy - height * 0.5
      n3y = n1y
      n4y = n2y
      n5y = n2y

      tupels = []
      tupels.append((n1x, n1y))
      if x == 0: # left boundary
        tupels.append((n4x, n4y))
      tupels.append((n2x, n2y))
      if x == nx-1: # right boundary
        tupels.append((n5x,n5y))
      tupels.append((n3x, n3y))

      #draw.polygon(tupels, fill=(0,0,int(np.random.rand(1)[0]*256))) #colored for debug
      draw.polygon(tupels, fill="black")

      # v = 2/3 is massive material, everything below is porous
      # due to rounding we substract the size of half a pixel
      if v < 2.0/3.0 - 1/res[0]/2:
        # white inner triangle
        # assumes equilateral triangle with 60° angle
        # v/2*h / offsetx = tan(60°/2) = 1/sqrt(3)
        offsetx = np.sqrt(3) * v * 0.5 * height
        offsety = np.sqrt(pow(offsetx,2) + pow(v * 0.5 * height,2)) # assumes equilateral triangle
        n1x = np.floor(mid[0] * dx) - np.ceil(length) + np.floor(offsetx)
        n2x = np.floor(mid[0] * dx)
        n3x = np.floor(mid[0] * dx) + np.ceil(length) - np.floor(offsetx)

        # raise a ValueError, if the radius is too large for the parameter
        # too large radius would result no longer in a triangle with round vertices
        if rad > (n3x - n1x)/2/sqrt(3):
          if param:
            raise ValueError('Bending parameter to large')
          else:
            continue
            #rad = np.floor((n3x - n1x)/2/sqrt(3))

        if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
          n1y = mid[1] * dy - height * 0.5 + height * 0.5 * v
          n2y = mid[1] * dy + height * 0.5 - offsety
          n3y = n1y
          height_tria = np.abs(n1y - n2y)
          mid_tria = (n2x, n2y - 2/3 * height_tria)
        else:
          n1y = mid[1] * dy + height * 0.5 - height * 0.5 * v
          n2y = mid[1] * dy - height * 0.5 + offsety
          n3y = n1y
          height_tria = np.abs(n1y - n2y)
          mid_tria = (n2x, n2y + 2/3 * height_tria)

        tupels = []
        tupels.append((n1x, n1y))
        tupels.append((n2x, n2y))
        tupels.append((n3x, n3y))

        if v > 0 and rad is not None and rad > 0:
          new_tupels = []
          for idx in range(3):
            # draw rounded corner
            a1, a2 = draw_chord(x, y, tupels, idx, mid_tria, rad)
            
            # add end points of arc to form new white inner "triangle" (actually hexagon)
            if (x % 2 == 0 and y % 2 == 0) or (x % 2 == 1 and y % 2 == 1):
              # tupel goes counterclockwise -> switch a1 and a2
              a1, a2 = a2, a1
            new_tupels.append((a1[0],a1[1]))
            new_tupels.append((a2[0],a2[1]))
          tupels = new_tupels

        # output for meshing with triangles
        #with open('matviz.out', 'a+') as f:
        #  print('tria {}'.format(tupels), file=f)

        #draw.polygon(tupels, fill=(0,int(np.random.rand(1)[0]*256),0)) #colored for debug
        draw.polygon(tupels, fill="white")

        # white triangle at boundary
        if x == 0 or x == nx-1:
          offsetx = np.sqrt(3) * v * 0.5 * height

          fac = 1 if x == 0 else -1
          n1x = mid[0] * dx - fac * length * 0.5
          n2x = np.round(mid[0] * dx) - np.round(fac * offsetx)
          n3x = n1x

          # n1x has to be at the boundary
          assert((x == 0 and np.round(n1x) == 0) or (x == nx-1 and np.round(n1x) == res[0]))

          if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
            # 1--2    2--1
            # | /      \ |
            # |/        \|
            # 3          3 
            n1y = mid[1] * dy + height * 0.5 - height * 0.5 * v
            n2y = n1y
            n3y = np.floor(mid[1] * dy + height * v)
          else:
            # upside down
            n1y = mid[1] * dy - height * 0.5 + height * 0.5 * v
            n2y = n1y
            n3y = np.ceil(mid[1] * dy - height * v)

          # for some v the boundary triangle might actually not exist 
          if (x == 0 and n2x-n1x < 1e-10) or (x == nx-1 and n1x-n2x < 1e-10):
            continue

          tupels = []
          tupels.append((n1x, n1y))
          tupels.append((n2x, n2y))
          tupels.append((n3x, n3y))

          if v > 0 and rad is not None and rad > 0:
            # we have to calculate a point on the angle bisector
            vec1 = np.array(tupels[0]) - np.array(tupels[1])
            vec2 = np.array(tupels[2]) - np.array(tupels[1])
            vec1 = vec1
            vec2 = vec2/np.linalg.norm(vec2) * np.linalg.norm(vec1)
            p = tupels[1] + (vec1 + vec2) / 2

            # draw rounded corner
            a1, a2 = draw_chord(x, y, tupels, 1, p.tolist(), rad, True)

            # tupel goes counterclockwise -> switch a1 and a2
            if x == 0: #left boundary
              if y % 2 == 0:
                a1, a2 = a2, a1
            else: #right boundary
              if (x % 2 == 0 and y % 2 == 1) or (x % 2 == 1 and y % 2 == 0):
                a1, a2 = a2, a1

            new_tupels = []
            new_tupels.append(tuple(np.round(tupels[0]).tolist()))
            new_tupels.append(tuple(a1.tolist()))
            new_tupels.append(tuple(a2.tolist()))
            new_tupels.append(tuple(np.round(tupels[2]).tolist()))
            tupels = new_tupels

          #draw.polygon(tupels, fill=(int(np.random.rand(1)[0]*256),0,0)) #colored for debug
          draw.polygon(tupels, fill="white")

  # the coordinate system of the image is centered at the upper left corner
  # with +x from left to right and +y from top to bottom
  # however, our meshes are centered at the lower left corner
  # with +x from left to right and +y from bottom to top
  # thus we flip the image upside down
  im = ImageOps.flip(im)

  if param is not None:
    # crop argument = (left, upper, right, lower)
    # remove half of an element on both sides (times dx for pixels)
    # +1 for periodicity:
    # will add one row of black pixels in y direction (cut happens outside of image region)
    # but not in x direction (cut happens inside of image region)
    if type == 'hourglass':
      im = im.crop((maxi[0]/nx*1/2*dx, mini[1]*dy, maxi[0]/nx*(nx-1/2)*dx+1, maxi[1]*dy+1))
    elif type == 'diamond':
      im = im.crop((maxi[0]/nx*1/2*dx, mini[1]*dy, maxi[0]/nx*(nx-1/2-(repetitions-1))*dx, maxi[1]*dy))
    
  return im

def draw_verts(verts, axsubplot, color):
  codes = [Path.MOVETO] + ((len(verts) - 1) * [Path.LINETO]) + [Path.CLOSEPOLY]
  verts = list(verts)
  verts.append((0, 0))

  path = Path(verts, codes)
  patch = matplotlib.patches.PathPatch(path, edgecolor='none', facecolor=color, lw=1)
  axsubplot.add_patch(patch)

def draw_circle(center, radius, axsubplot, col):
  circle = matplotlib.pyplot.Circle(center, radius, color=col, fill=True)
  axsubplot.add_artist(circle)

# draws a thick circle, where the thickness is determined automatically by the radius
def draw_thick_circle(draw, center, radius):

  for o in numpy.arange(0.0, 0.01 * radius, 0.05):
    draw.ellipse((center[0] - radius + o, center[1] - radius + o, center[0] + radius - o, center[1] + radius - o), fill=None, outline="black")


# visualize the orientational stiffness
# @return the image
def show_orientational_stiffness(coords, angle, data, nx, scale=-1.0, axes=False):

  centers, min, max, elem = coords

  max_val = numpy.max(data[:])
  min_val = numpy.min(data[:])
  data_offset = -2 * min_val if min_val < 0 else 0

  im, draw, dim, dx, dy = create_image(min, max, nx)

  if scale == -1.0:
    dist = 1.0 if len(centers) == 1 else centers[1][0] - centers[0][0]
    scale = 0.35 * dx * dist / max_val

  sm = cmx.ScalarMappable(colors.Normalize(vmin=min_val, vmax=max_val), cmap=plt.get_cmap('jet'))

  if numpy.ndim(data) == 1:
    data = [data]
    angle = [angle]

  for vals, angles, coord in zip(data, angle, centers):
    angle_max = angles[numpy.argmax(numpy.abs(vals))]
    angle_min = angles[numpy.argmin(numpy.abs(vals))]

    #print(" largest entry: {:>13.6e}".format(max_val) + "  in direction " + str(to_vector(angle_max)))
    #print("smallest entry: {:>13.6e}".format(min_val) + "  in direction " + str(to_vector(angle_min)))

    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy

    pol = to_polygons(angles, vals, x_off, dim[1] - y_off, scale, True)
    m = numpy.max(vals)
    if axes:
      pol.append(pol[0])
      draw.line(pol, fill="red", width=3)
    else:
      draw.polygon(pol, fill=color_code(sm, m), outline="black")
    if min_val < 0:
      pol = to_polygons(angles, vals, x_off, dim[1] - y_off, scale, False)
      draw.polygon(pol, fill="gray", outline="black")

  if axes:
    # draw circles for absolute value
    num_axes = 6
    max = numpy.max(data[:])
    mult = pow(10,-np.floor(np.log10(max/num_axes)))
    axes_step = np.ceil(max*mult/num_axes) / mult
    num_axes = int(np.ceil(max/axes_step))
    x_center = (coord[0] + min[0]) * dx
    y_center = (coord[1] + min[1]) * dy
    for r in range(1, num_axes + 1):
      radius = r * axes_step * scale
      twoPointList = [x_center-radius, y_center-radius, x_center+radius, y_center+radius]
      draw.ellipse(twoPointList, outline="black", width=1)
      label = str(r*axes_step)
      pos = (x_center + radius/sqrt(2), y_center - radius/sqrt(2))
      
      fnt = None
      try:
        fnt = ImageFont.truetype("arial.ttf", int(dx/50))
      except OSError:  
        fnt = ImageFont.truetype("Geneva.ttf", int(dx/50))
      draw.text(pos, label, anchor="lb", fill="black", font=fnt)
    # draw lines for angles
    for seg in range(12):
      phi = seg * -30 / 180 * numpy.pi
      end = (x_center + radius * numpy.cos(phi), y_center + radius * numpy.sin(phi))
      draw.line((x_center, y_center) + end, fill="black")
      pos = (x_center + radius * numpy.cos(phi) * 1.05, y_center + radius * numpy.sin(phi) * 1.05)
      draw.text(pos, str(seg*30), anchor="mm", fill = "black", font=fnt)

  return im


# @param aux see  perform_voigt_tensor_samping()
# @return list of angles and list of data which might be aux
def perform_rotations(tensors, samples, aux_code=None):

  # tensors contains just one tensor
  if numpy.ndim(tensors) == 2:
    angle, data, aux = perform_voigt_tensor_samping(tensors, samples, aux_code)
    return angle, data if aux_code == None else aux

  # tensors is an array of tensors
  res_angle = []
  res_data = []
  for i in range(len(tensors)):
    tensor = tensors[i]
    assert(numpy.ndim(tensor) == 2)
    angle, data, aux = perform_voigt_tensor_samping(tensor, samples, aux_code)
    res_angle.append(angle)
    res_data.append(data if aux_code == None else aux)

  return res_angle, res_data
