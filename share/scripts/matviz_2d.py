from matviz_rot import *
from hdf5_tools import *

import Image, ImageDraw
import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
# matplotlib.use('Agg')
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
from scipy import ndimage
import numpy as np
import scipy.interpolate as ip
import math


# # create and prepare a matplot figure where patched might be "added" to
def create_figure(min, max, res, for_save):
  # we set the aspect ration and also the resolution such that we can export as png
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

# # @param centers barycenters
# @param min/max minimal/maximal real node (not barycenter)
# @retuen image, draw, dim of image, dx to scale from node coord to image coords  
def create_image(min, max, nx, color="white"):
  dim = (nx, int(nx * (max[1] + min[1]) / (max[0] + min[0])))
  
  dx = dim[0] / (max[0] + 2.0 * min[0])
  dy = dim[1] / (max[1] + 2.0 * min[1])
 
  # print "dx=" + str(dx) + " dy=" + str(dy) + " nx=" + str(nx) + " dim=" + str(dim) + " max=" + str(max) 
  im = Image.new("RGB", dim, color)
  draw = ImageDraw.Draw(im)
  
  return im, draw, dim, dx, dy  


# # @return phi, r
def to_polar(x, y):
  return numpy.sqrt(x ** 2 + y ** 2), numpy.arctan2(y, x)

# # polar coordiantes to cartesian
def to_cart(phi, r):
  return r * numpy.cos(phi), r * numpy.sin(-phi)
# calculate volume for s1, s2. Assume regular grid!
def calc_volume(s1, s2):
  vol = 0.0
  for i in range(len(s1)):
    vol += s1[i] + s2[i] - s1[i] * s2[i]
  return (vol / len(s1))[0]  # somehow this is a numpy.ndarry


# # generate polygon vertices out of rotation data
# to be applied as draw.polygon(result, fill="green", outline="black")
#
# from paraview_fmo import *
# c0 = [12.6, 11.7, 2.3, 0, 0, 8.41]
# import Image, ImageDraw
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
    y = r * numpy.sin(-angle[i])
    tupl.append((x_offset + scale * x, y_offset + scale * y))
  return tupl    

# # give the corners to draw a rotated rectangles as polygon
def to_rectangle_center(height, width, angle, x_offset, y_offset):
  
  # print "h=" + str(height) + " w=" + str(width) + " a=" + str(angle) + " x=" + str(x_offset) + " y=" + str(y_offset)
  
  # height = numpy.max((height,1))
  
  tupl = []

  for x in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]:
    p = (x[0] * width / 2, x[1] * height / 2)
    r = (cos(angle) * p[0] - sin(angle) * p[1], sin(angle) * p[0] + cos(angle) * p[1])
    # r = (cos(angle) * p[0] + sin(angle)*p[1], - sin(angle) * p[0] + cos(angle) * p[1])
    tupl.append((x_offset + r[0], y_offset + r[1]))

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

# # give the corners to draw a rotated rectangles as polygon
def to_rectangle_corner(lower, upper):
  
  # print "to_rectangle_corner " + str(lower) + " -> " + str(upper)
  
  tupl = []

  tupl.append(lower)
  tupl.append((upper[0], lower[1]))
  tupl.append(upper)
  tupl.append((lower[0], upper[1]))
  
  return tupl  

# # helper for show_rot_frame_grad
def get_interpol_data(coords, data, fallback, x, eval=True):
  if not eval:
    return None, None   
  v = data[x]
  if v[0] == -1.0:
    v = fallback[x]
  return coords[x], v  

# # helper for get_interpolation and Fields
# @return 2d locations and 2d data
def convert_two_data_interpolation_input(centers, s1, s2, angle):
  # convert to 2D
  c = numpy.zeros((len(centers), 2))
  c[:, 0] = [x[0] for x in centers]
  c[:, 1] = [x[1] for x in centers]

  v = numpy.zeros((len(s1), 2 if angle == None else 3))
  for i in range(len(s1)):
    v[i][0] = s1[i][0]
    v[i][1] = s2[i][0]
    if angle <> None:
      v[i][2] = angle[i]
      
  return c, v    



# # helper for get_interpolation and Fields
# @return 2d locations and 2d data
def convert_single_data_interpolation_input(centers, s1, angle):
  # convert to 2D
  c = numpy.zeros((len(centers), 2))
  c[:, 0] = [x[0] for x in centers]
  c[:, 1] = [x[1] for x in centers]

  v = numpy.zeros((len(s1), 1 if angle == None else 2))
  for i in range(len(s1)):
    v[i][0] = s1[i][0]
    if angle <> None:
      v[i][1] = angle[i]
      
  return c, v    


# # helper which returns an interpolated grid and one nearest neighbor interpolated grid
def get_interpolation(coords, grad, sample, s1, s2, angle=None):
  assert(sample == 'elem_nodes' or sample == 'edge_centers')
  
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
  else:  # 'edge_centers' this is much more complicated, draw a coarse grid and mark all element edge centers!
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

  c, v = convert_two_data_interpolation_input(centers, s1, s2, angle)
    
  ip_data = ip.griddata(c, v, out, grad, -1.0)
  # any interpolation but nearest neighbor can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation
  ip_near = ip.griddata(c, v, out, 'nearest') if grad <> 'nearest' else None
  
  return ip_data, ip_near, out, nx, ny 

# # visualize the orientational stiffness
# @param grad is 'none' or 'linear'
# @return the image
def show_frame_grad(coords, s1, s2, grad, direction, nx):

  centers, min, max, elem = coords
  im, draw, dim, dx, dy = create_image(min, max, nx, "white")

  height = elem[1] * dy 
  length = elem[0] * dx
 
  # print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " height=" + str(height) + " length=" + str(length) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'elem_nodes', s1, s2)
  

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

# # visualize the orientational stiffness
# @return the image
def show_rot_cross_grad(coords, s1, s2, angle, grad, direction, nx, scale, do_save):

  centers, min, max, elem = coords
  fig, sub = create_figure(min, max, nx, do_save)
  delta_angle = numpy.max(angle[:]) - numpy.max(angle[:]) 

  # print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'edge_centers', s1, s2, angle)
  if scale == -1:
    scale = 1.0 if delta_angle == 0.0 else 0.8 

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
  
# visualizes the oriental stiffness as frame with smooth inner corners; creates a vector image
def show_modified_frame(coords, s1, s2, angle, direction, nx, scale, color, do_save):
  print 'image is only correct if eps<= s1,s2 <= 0.5, otherwise scaling is necessary; rotation is not implemented currently'
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
    pol = to_rectangle_center(length, length, theta + vmin * numpy.pi / 2, x_off, y_off)
    draw_verts(pol, sub, 'black')
    
    pol = to_rectangle_center(length * (1. - 2.*v[0]), length * (1. - 2.*v[1]), theta + vmax * numpy.pi / 2, x_off, y_off)
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
    pol = to_rectangle_center(r, r, theta + vmax * numpy.pi/2, x_off-t[0], y_off-t[1])
    #pol[:][0] -= t[0]
    #pol[:][1] -= t[1]
    draw_verts(pol, sub, 'black')
    c = (x_off - c1, y_off - c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')
  
    #t = (cos(theta + vmax * numpy.pi / 2) * abs(l1 / 2. - r / 2.) - sin(theta + vmax * numpy.pi / 2) *abs(l2 / 2. - r / 2.), sin(theta + vmax * numpy.pi / 2) * abs(l1 / 2. - r / 2.) + cos(theta + vmax * numpy.pi / 2) * abs(l2 / 2. - r / 2.))
    # lower right corner
    pol = to_rectangle_center(r, r, 0., x_off + t[0], y_off - t[1])
    draw_verts(pol, sub, 'black')
    c = (x_off + c1, y_off - c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')
    
    # upper right corner
    pol = to_rectangle_center(r, r, theta + vmax * numpy.pi / 2, x_off + t[0], y_off + t[1])
    draw_verts(pol, sub, 'black')
    c = (x_off + c1, y_off + c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')
    
    # upper left corner
    pol = to_rectangle_center(r, r, theta + vmax * numpy.pi / 2, x_off - t[0], y_off + t[1])
    draw_verts(pol, sub, 'black')
    c = (x_off - c1, y_off + c2)
    center = (cos(theta) * c[0] - sin(theta) * c[1], sin(theta) * c[0] + cos(theta) * c[1])
    draw_circle(center, r-eps, sub, 'white')    
  return (fig, sub)  

# # visualize the orientational stiffness
# @param grad is 'none' or 'linear'
# @return the image
def show_frame(coords, s1, s2, directions, nx):

  centers, min, max, elem = coords
  im, draw, dim, dx, dy = create_image(min, max, nx, "white")
  height = elem[1] * dy 
  length = elem[0] * dx
 
  for i in range(len(s1)):
  
    coord = centers[i]
    
         
    x_off = (coord[0] + min[0] - 0.5 * elem[0]) * dx 
    y_off = (coord[1] + min[1] - 0.5 * elem[1]) * dy
    hor = s2[i, 0]  # it seems that stiff1 and stiff2 are mixed up. This tries to correct it
    ver = s1[i, 0]
    # print "hor=" + str(hor) + " ver=" + str(ver) 
    
    if not directions == 'vertical': 
      # lower horizontal line  
      pol = to_rectangle_corner((x_off, dim[1] - y_off), (x_off + length, dim[1] - y_off - height * 0.5 * ver - 0.5))
      draw.polygon(pol, fill="black")
  
      # upper horizontal line
      pol = to_rectangle_corner((x_off, dim[1] - y_off - height + height * 0.5 * ver - 0.5), (x_off + length, dim[1] - y_off - height))
      draw.polygon(pol, fill="black")

    if not directions == 'horizontal':
      # left vertical line
      pol = to_rectangle_corner((x_off, dim[1] - y_off), (x_off + length * 0.5 * hor + 0.5, dim[1] - y_off - height))
      draw.polygon(pol, fill="black")
  
      # right vertical line
      pol = to_rectangle_corner((x_off + length - length * 0.5 * hor - 0.5, dim[1] - y_off), (x_off + length, dim[1] - y_off - height))
      draw.polygon(pol, fill="black")

  return im  


# @return the image
def show_rot_cross(coords, s1, s2, angle, direction, nx, scale, color, do_save):

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

    # print 'S=' + str(s1[i,0]) + '/' + str(s2[i,0])  + ' v=' + str(v) + ' c=' + str(c)

    # a
    if direction == 'horizontal': 
      pol = to_rectangle_center(length * v[0], length, theta, x_off, dim[1] - y_off)
      draw_verts(pol, sub, c[0])    
    # b
    elif direction == 'vertical':
      pol = to_rectangle_center(length * v[1], length, theta + numpy.pi / 2, x_off, y_off)
      draw_verts(pol, sub, c[1])
    else:
      # vmax, vmin are indices which decide whether s1-rectangle or s2-rectangle are drawn first 
      vmax,vmin = (0,1) if v[0] > v[1] else (1,0)
      pol = to_rectangle_center(length * v[vmin], length, theta + vmin * numpy.pi / 2, x_off, y_off)
      # draw_verts(pol, sub, str(1.0 - c[vmin]))
      draw_verts(pol, sub, c[vmin])
      pol = to_rectangle_center(length * v[vmax], length, theta + vmax * numpy.pi / 2, x_off, y_off)
      draw_verts(pol, sub, c[vmax])
  return (fig, sub)


def color_code(color_map, value):
  c = color_map.to_rgba(value)
  return "rgb(" + str(int(255 * c[0])) + ", " + str(int(255 * c[1])) + "," + str(int(255 * c[2])) + ")"
  
  
def draw_verts(verts, axsubplot, color):
  verts.append((0, 0))
  codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]

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

  
# # visualize the orientational stiffness
# @return the image
def orientational_stiffness(coords, angle, data, nx, scale=-1.0):

  centers, min, max, elem = coords

  max_val = numpy.max(data[:])
  min_val = numpy.min(data[:])
  data_offset = -2 * min_val if min_val < 0 else 0
     
  im, draw, dim, dx, dy = create_image(min, max, nx)   
   
  print "max=" + str(max_val) + " min=" + str(min_val)
   
  if scale == -1.0:
    dist = 1.0 if len(centers) == 1 else centers[1][0] - centers[0][0]
    scale = 0.35 * dx * dist / max_val
     
  sm = cmx.ScalarMappable(colors.Normalize(vmin=min_val, vmax=max_val), cmap=plt.get_cmap('jet'))
    
  for i in range(len(data)):
    coord = centers[i]
    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy

    pol = to_polygons(angle[i], data[i], x_off, dim[1] - y_off, scale, True)
    m = numpy.max(data[i])
    draw.polygon(pol, fill=color_code(sm, m), outline="black")
    if min_val < 0:
      pol = to_polygons(angle[i], data[i], x_off, dim[1] - y_off, scale, False)
      draw.polygon(pol, fill="gray", outline="black")

  return im  
  
  
# @param aux see  perform_cfs_rotation()
# @return list of angles and list of data which might be aux   
def perform_rotations(tensors, notation, samples, name="mechTensor", aux_code=None):

  res_angle = []
  res_data = []
  for i in range(len(tensors)):
    t = tensors[i]
    # print "pr: t=" + str(t)
    tensor = 0
    if name == "mechTensor":
      if notation == "voigt":
        tensor = to_mech_tensor(t.transpose())
      else: 
        tensor = HillMandel2Voigt(to_mech_tensor(t.transpose()))
    if name == "elecTensor":
      tensor = to_elec_tensor(t.transpose())
    if name == "piezoTensor":
      tensor = to_piezo_tensor(t.transpose())

    angle, data, aux = perform_cfs_rotation(tensor, samples, aux_code)
    res_angle.append(angle)
    res_data.append(data if aux_code == None else aux)

  return res_angle, res_data  
    
# import pylab
# pylab.plot(data[:,0],data[:,1])
# pylab.show()  
 
# tensor = []
# t = to_mech_tensor(eval("[1.0,0.5,0.0,0.0,0.0,0.0]"))
# tensor.appt.end(to_mech_vector(t, as_array=True)) 
# centers = []
# centers.append([0.0,0.0,0.0])
# centers.append([1.0,1.0,0.0])

# angle, data = perform_rotations(tensor, 10)

# f = h5py.File("/home/fwein/project/simp/hook.h5")
# s1 = get_element(f, "design_stiff1_plain", "mech")
# s2 = get_element(f, "design_stiff2_plain", "mech") 
# angle = get_element(f, "design_rotAngle_plain", "mech")
 
# t = to_mech_tensor(eval("[2.607069, 1.484705, 0.1626158, 0,0, 0.3030707]")) 
 
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")
# f = h5py.File("/home/fwein/tmp/l_sl-m_20-g_al-p_0.1-gamma_1e-07-tau_0.01.h5")
# r  = centered_elements(f)
# tensor = get_element(f, "mechTensor", "piezo")
# tensor = get_element(f, "mechTensor", "mech")
# rot = perform_rotations(tensor, 10)
# rot = perform_rotations(tensor, t21=True)
# orientational_stiffness(r, rot, 1500, 2.5).show()
