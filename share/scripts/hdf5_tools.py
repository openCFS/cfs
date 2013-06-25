import sys
import h5py
import numpy
import operator
import Image, ImageDraw
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import ndimage
import numpy as np
import scipy.interpolate as ip 
from paraview_fmo import *

# open a hdf5 file as
# f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")

## print tensor, pyhton makes for 6x6 tensors too early line breaks
def dump_tensor(tensor):
  for y in range(tensor.shape[0]):
    sys.stdout.write(str(y+1) + ": ")
    for x in range(tensor.shape[1]):
      #sys.stdout.write(str(tensor[y][x]) + "\t")
      print '%(val)10.4g ' % {"val": tensor[y][x]},
    print ""    
  
# checks if the region exists
# exits with an message if not
def validate_region(hdf5_file, region):
  regions = hdf5_file['/Mesh/Regions']
  if not any(k in regions.keys() for k in [region]):
    print "region '" + region + "' not within regions " + str(regions.keys())

## give back elements with barycenters
# assumes rectangles
# @return list barycenter tuple ordered by elements and min and max node coordinates and first element dimensions
def centered_elements(hdf5_file, region):
  all_elements = hdf5_file['/Mesh/Elements/Connectivity'] # for all regions
  reg_elements = hdf5_file['/Mesh/Regions/' + region + '/Elements']
  types = hdf5_file['/Mesh/Elements/Types']
  all_nodes = hdf5_file['/Mesh/Nodes/Coordinates']
  reg_nodes = hdf5_file['/Mesh/Regions/' + region + '/Nodes']
  
  # determine elem_dim from first element dimensions, ignore region
  node_coords = []
  for n in range(len(all_elements[0])):
    node_coords.append(all_nodes[all_elements[0][n]-1]) # numbers are one-based
  ma = np.array([max(node_coords,key=operator.itemgetter(0))[0], max(node_coords,key=operator.itemgetter(1))[1],  max(node_coords,key=operator.itemgetter(2))[2]])
  mi = np.array([min(node_coords,key=operator.itemgetter(0))[0], min(node_coords,key=operator.itemgetter(1))[1],  min(node_coords,key=operator.itemgetter(2))[2]])
  elem_dim = ma - mi
    
  # determine region dimensions, we need to resort for the desired region! Due to 1 to zero based conversion we need to do it manually :(
  nodes = numpy.zeros((len(reg_nodes), 3))
  for e in range(len(reg_nodes)):
    nodes[e] = all_nodes[reg_nodes[e] - 1]  
  min_dim = min(nodes[:,0]), min(nodes[:,1]) 
  max_dim = max(nodes[:,0]), max(nodes[:,1])  
    
  result = []
  for e in range(len(reg_elements)):
    idx = reg_elements[e] - 1 # cfs writes one based
    if types[idx] == 6:
       nod = all_elements[idx]
       center = numpy.array([0.0, 0.0, 0.0])
       for n in range(len(nod)):
         center += all_nodes[nod[n]-1] # numbers are one-based
         # print "el=" + str(e) + " n=" + str(n) + " node=" + str(nod[n]) + "->" + str(nodes[nod[n]-1]) + " center=" + str(center) 
       center *= 1.0/len(nod)
       result.append(center)
       # print "e=" + str(e) + " idx=" + str(idx) + " nod=" + str(nod) + " center=" + str(center) 
    else:
      assert(False) # 3D? not implented yet   

  return result, min_dim, max_dim, elem_dim     
                
## find minimal and maximal coordinate
# @param coordinates as from centered_elements
def find_corners(centers):
  min = [1e30, 1e30, 1e30]
  max = [-1e30, -1e30, -1e30]
  
  if len(centers) == 1:
    min = [0.0,0.0,0.0]
    max = [1.0,1.0,0.0]
    
  
  for e in range(len(centers)):
    test = centers[e]
    
    for c in range(3):
      if test[c] < min[c]:
        min[c] = test[c]
      if test[c] > max[c]:
        max[c] = test[c]
  return min, max      
    
           
def get_element(hdf5_file, name, region, step=99999):
  ms = hdf5_file['/Results/Mesh/MultiStep_1']
  # for optimization we assume steps to be numbered from 0, in simulation they start from 1
  if step >= len(ms):
    step = max((len(ms) - 2,0)) # reset to last, first element is ResultDescription
  key = "/Results/Mesh/MultiStep_1/Step_" + str(step) + "/" + name + "/" + region + "/Elements/Real"
  try:
    data = ms[key]
    return data
  except:
    raise Exception("cannot access '" + key + "' in " + str(hdf5_file.filename))
  

## @return phi, r
def to_polar(x, y):
  return numpy.sqrt(x**2 + y**2), numpy.arctan2(y, x)

## polar coordiantes to cartesian
def to_cart(phi, r):
  return r * numpy.cos(phi), r * numpy.sin(-phi)

## generate polygon vertices out of rotation data
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
def to_polygons(angle, data, x_offset, y_offset, scale):
  tupl = []
  for i in range(len(data)):
    r = numpy.abs(data[i])
    x = r * numpy.cos(angle[i])
    y = r * numpy.sin(-angle[i])
    tupl.append((x_offset + scale * x,y_offset + scale * y))
  return tupl    

## give the corners to draw a rotated rectangles as polygon
def to_rectangle_center(height, width, angle, x_offset, y_offset):
  
  # print "h=" + str(height) + " w=" + str(width) + " a=" + str(angle) + " x=" + str(x_offset) + " y=" + str(y_offset)
  
  height = numpy.max((height,1))
  
  tupl = []

  for x in [(-1.0,-1.0),(1.0,-1.0),(1.0,1.0),(-1.0,1.0)]:
    p = (x[0] * width/2, x[1] * height/2)
    r = (cos(angle) * p[0] -sin(angle)*p[1], sin(angle) * p[0] + cos(angle) * p[1])
    # r = (cos(angle) * p[0] + sin(angle)*p[1], - sin(angle) * p[0] + cos(angle) * p[1])
    tupl.append((x_offset + r[0], y_offset + r[1]))

  return tupl  

# draws a rotated frustum
# @param direction either vertical or horizontal, not both!
def to_frustum_center(start, end, center, elem_scale, elem, scale, direction):
  # 4 ------- 3
  # |         |
  # |         |
  # 1---------2


  angle = 0.5 * (start[2] + end[2])
    # print "horizontal line
  if direction == 'vertical':
    angle += 0.5 * numpy.pi
    
  idx = 0 if direction == 'vertical' else 1
  alt_idx = 1 if direction == 'vertical' else 0  
    
  val_1 = start[idx] 
  val_2 = end[idx]
  
  tupl = []
  points = []
  
  #print start
  #print end
  #print elem
  
  #WARNING CHECK the scaling!!!
  
  # we scale the element scale, such it overlaps. therefore wen need to scale reciproc with the vale, which cancels!
  points.append((-1.0 * scale * elem[alt_idx]/2, 1./max((1, scale)) * -val_1 * scale * elem[idx]/2))
  points.append(( 1.0 * scale * elem[alt_idx]/2, 1./max((1, scale)) * -val_2 * scale * elem[idx]/2))
  points.append(( 1.0 * scale * elem[alt_idx]/2, 1./max((1, scale)) *  val_2 * scale * elem[idx]/2))
  points.append((-1.0 * scale * elem[alt_idx]/2, 1./max((1, scale)) *  val_1 * scale * elem[idx]/2))
  
  for i in range(4):
    #print "i=" + str(i + 1) + " -> " + str(points[i])
    r = (cos(angle) * points[i][0] -sin(angle)*points[i][1], sin(angle) * points[i][0] + cos(angle) * points[i][1])
    tupl.append(((center[0] + r[0]) * elem_scale[0], (center[1] + r[1]) * elem_scale[1]))
  
  return tupl

## give the corners to draw a rotated rectangles as polygon
def to_rectangle_corner(lower, upper):
  
  # print "to_rectangle_corner " + str(lower) + " -> " + str(upper)
  
  tupl = []

  tupl.append(lower)
  tupl.append((upper[0], lower[1]))
  tupl.append(upper)
  tupl.append((lower[0], upper[1]))
  
  return tupl  

  
## @param centers barycenters
# @param min/max minimal/maximal real node (not barycenter)
# @retuen image, draw, dim of image, dx to scale from node coord to image coords  
def create_image_new(centers, min, max, nx, color = "white"):
  dim = (nx, int(nx *  (max[1] + min[1]) / (max[0] + min[0])))
  
  dx = dim[0] / (max[0] + 2.0 * min[0])
  dy = dim[1] / (max[1] + 2.0 * min[1])
 
  # print "dx=" + str(dx) + " dy=" + str(dy) + " nx=" + str(nx) + " dim=" + str(dim) + " max=" + str(max) 
  im = Image.new("RGB", dim, color)
  draw = ImageDraw.Draw(im)
  
  return im, draw, dim, dx, dy  
  
def create_image(centers, nx,color = "white"):
  # zoom gives proper offsets of the elements
  min, max = find_corners(centers) # we assume the real grid to start at 0/0
  
  dim = (nx, int(nx *  (max[1] + min[1]) / (max[0] + min[0]))) 
  
  dx = dim[0] / (max[0] + 2.0 * min[0])
  dy = dim[1] / (max[1] + 2.0 * min[1])
  
  im = Image.new("RGB", dim, color)
  draw = ImageDraw.Draw(im)
  
  return im, draw, dim, dx, dy, min, max
  
## helper for show_rot_frame_grad
def get_interpol_data(coords, data, fallback, x, eval = True):
  if not eval:
    return None, None   
  v = data[x]
  if v[0] == -1.0:
    v = fallback[x]
  return coords[x], v  

## helper which returns an interolated grid and one nearest neighbor interpolated grid

def get_interpolation(coords, grad, sample, s1, s2, angle = None):
  assert(sample == 'elem_nodes' or sample == 'edge_centers')
  
  centers, min, max, elem = coords
  # convert to 2D
  c = numpy.zeros((len(centers), 2))
  c[:,0] = [x[0] for x in centers]
  c[:,1] = [x[1] for x in centers]
 
  # where we want nodes
  nx = int((max[0] - min[0]) / elem[0])
  print str((max[1] - min[1]) / elem[1])
  ny = int((max[1] - min[1]) / elem[1])
  out = None
  if sample == 'elem_nodes':
    out = numpy.zeros(((nx + 1) * (ny + 1), 2))
    for y in range(ny+1):
      for x in range(nx+1):
        out[y * (nx+1) + x][0] = float(x)/nx * max[0]
        out[y * (nx+1) + x][1] = float(y)/ny * max[1]
  else: # 'edge_centers' this is much more complicated, draw a coarse grid and mark all element edge centers!
    # ---5-------6----
    # |      |       | 
    # 2      3       4
    # |      |       |
    # ---0-------1----
    
    out = numpy.zeros(((2 * nx + 1) * (2 * ny + 1), 2))
    # 0,1,5,6
    for y in range(ny+1):
      for x in range(nx):
        out[y * (2*nx+1) + x][0] = float(x)/nx * max[0] + 0.5 * elem[0]
        out[y * (2*nx+1) + x][1] = float(y)/ny * max[1]
        # print "out[" + str(y * (2*nx+1) + x) + "] = " + str(out[y * (2*nx+1) + x])
    # 2.3.4
    for y in range(ny):
      for x in range(nx+1):
        out[nx + y * (2*nx+1) + x][0] = float(x)/nx * max[0]
        out[nx + y * (2*nx+1) + x][1] = float(y)/ny * max[1] + 0.5 * elem[1]
        # print "out[" + str(nx + y * (2*nx+1) + x) + "] = " + str(out[nx + y * (2*nx+1) + x])

  v = numpy.zeros((len(s1),  2 if angle == None else 3))
  for i in range(len(s1)):
    v[i][0] = s1[i][0]
    v[i][1] = s2[i][0]
    if angle <> None:
      v[i][2] = angle[i][0]
    
  ip_data = ip.griddata(c, v, out, grad, -1.0)
  # any interpolatiob but nearest neighbor can only interpolate in the convex hull,
  # if the value is -1 we use the nearest interpolation
  ip_near = ip.griddata(c, v, out, 'nearest') if grad <> 'nearest' else None
  
  return ip_data, ip_near, out, nx, ny 

## visualize the orientational stiffness
# @param grad is 'none' or 'linar'
# @return the image
def show_frame_grad(coords, s1, s2, grad, direction, nx):

  centers, min, max, elem = coords

  im, draw, dim, dx, dy = create_image_new(centers, min, max, nx,"white")

  height = elem[1] * dy 
  length = elem[0] * dx
 
  # print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " height=" + str(height) + " length=" + str(length) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'elem_nodes', s1, s2)
  

  for y in range(ny+1):
    for x in range(nx+1):
      start, v_start = get_interpol_data(out, ip_data, ip_near, y * (nx+1) + x)
      right, v_right = get_interpol_data(out, ip_data, ip_near, y * (nx+1) + x + 1, eval = True if x < nx else False)
      upper, v_upper = get_interpol_data(out, ip_data, ip_near, (y + 1) * (nx+1) + x, eval = True if y < ny else False)
      
      # print "start=" + str(start) + " v_start=" + str(v_start) + " right=" + str(right) + " v_right=" + str(v_right)
      
      # 4 ------- 3
      # |         |
      # |         |
      # 1---------2

      # print horizontal line
      if not direction == 'vertical':
        if x < nx: # the right most vertical line has no right for the horizontal line
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
        if y < ny: # the evaluation of the most upper horizontal line has no upper for vertical data
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

## visualize the orientational stiffness
# @return the image
def show_rot_cross_grad(coords, s1, s2, angle, grad, direction, nx, scale=-1):

  centers, min, max, elem = coords

  im, draw, dim, dx, dy = create_image_new(centers, min, max, nx,"white") 

  delta_angle = numpy.max(angle[:,0]) - numpy.max(angle[:,0]) 

  #print "elem=" + str(elem) + " dx=" + str(dx) + " dy=" + str(dy) + " min=" + str(min) + " max=" + str(max)
  ip_data, ip_near, out, nx, ny = get_interpolation(coords, grad, 'edge_centers', s1, s2, angle)

  #print "nx="  + str(nx) + " ny=" + str(ny)
  
  if scale == -1:
    scale = 1.0 if delta_angle == 0.0 else 0.8 

  # out is rather complex set, see get_interpolation for details!
  if not direction == 'vertical':
    for y in range(ny):
      for x in range(nx):
        start, v_start = get_interpol_data(out, ip_data, ip_near, nx + y * (2*nx+1) + x)
        right, v_right = get_interpol_data(out, ip_data, ip_near, nx + y * (2*nx+1) + x + 1)
        
        center = ((0.5 * (start[0] + right[0]), max[1] - start[1]))
        
        pol = to_frustum_center(v_start, v_right, center, (dx, dy), elem, scale, 'horizontal') 
        draw.polygon(pol, fill="black")
  
  if not direction == 'horizontal':
    for y in range(ny):
      for x in range(nx):
        start, v_start = get_interpol_data(out, ip_data, ip_near, y * (2*nx+1) + x)
        upper, v_upper = get_interpol_data(out, ip_data, ip_near, (y+1) * (2*nx+1) + x)

        center = ((start[0], max[1] - 0.5 * (start[1] + upper[1])))

        #if y == 0:
        #  print "start=" + str(start) + " right=" + str(right) + " v_start=" + str(v_start) + " v_upper=" + str(v_upper)
        
        pol = to_frustum_center(v_upper, v_start, center, (dx, dy), elem, scale, 'vertical') 
        draw.polygon(pol, fill="black")



  return im  



## visualize the orientational stiffness
# @param grad is 'none' or 'linar'
# @return the image
def show_frame(coords, s1, s2, directions, nx):

  centers, min, max, elem = coords
  
  im, draw, dim, dx, dy = create_image_new(centers, min, max, nx,"white")

  height = elem[1] * dy 
  length = elem[0] * dx
 
  for i in range(len(s1)):
  
    coord = centers[i]
    
    # if coord[1] > 0.04:
    #  continue
          
    x_off = (coord[0] + min[0] - 0.5 * elem[0]) * dx 
    y_off = (coord[1] + min[1]  - 0.5 * elem[1]) * dy

    #print str(coord) + " " + str(x_off) + " " + str(y_off)

    hor = s1[i,0]  # it seems that stiff1 and stiff2 are mixed up. This tries to correct it
    ver = s2[i,0]

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


## visualize the orientational stiffness
# @return the image
def show_rot_cross(coords, s1, s2, angle, direction, nx, scale=-1):

  centers, min, max, elem = coords

  im, draw, dim, dx, dy = create_image_new(centers, min, max, nx,"white") 

  delta_angle = numpy.max(angle[:,0]) - numpy.max(angle[:,0]) 

  if scale == -1:
    scale = 1.02 if delta_angle == 0.0 else 0.8 

  length =  scale * (elem[0]) * dx
  
  for i in range(len(s1)):
  
    coord = centers[i]
    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy

    # we need downscale the values when we overscale due to overlapping 
    v1 = s1[i,0] / numpy.max((scale, 1.))
    v2 = s2[i,0] / numpy.max((scale, 1.))
    theta = angle[i,0]
    
    # b
    if not direction == 'horizontal':
      pol = to_rectangle_center(length * v2, length, theta, x_off, dim[1] - y_off) 
      draw.polygon(pol, fill="black")

    # a
    if not direction == 'vertical': 
      pol = to_rectangle_center(length * v1, length, theta + numpy.pi/2, x_off, dim[1] - y_off) 
      draw.polygon(pol, fill="black")

  return im  




def color_code(color_map, value):
  c = color_map.to_rgba(value)
  return "rgb(" + str(int(255 * c[0])) + ", " + str(int(255*c[1])) + "," + str(int(255*c[2])) + ")"
  
## visualize the orientational stiffness
# @return the image
def orientational_stiffness(centers, angle, data, nx, scale=-1):

  im, draw, dim, dx, dy, min, max = create_image(centers, nx)

  max_val = numpy.max(data[:])
  min_val = numpy.min(data[:])
   
  if scale == -1:
    dist = 1.0 if len(centers) == 1 else centers[1][0] - centers[0][0] 
    scale = 0.35 * dx * dist / max_val
   
  sm = cmx.ScalarMappable(colors.Normalize(vmin=min_val, vmax=max_val), cmap=plt.get_cmap('jet'))
    
  for i in range(len(data)):
    coord = centers[i]
    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy

    pol = to_polygons(angle[i], data[i], x_off, dim[1] - y_off, scale)
    m = numpy.max(data[i])
    draw.polygon(pol, fill=color_code(sm, m), outline="black")

  return im  


  
  
# @param aux see  perform_cfs_rotation()
# @return list of angles and list of data which might be aux   
def perform_rotations(tensors, samples, name = "mechTensor", aux_code = "default"):

  res_angle = []
  res_data  = []
  for i in range(len(tensors)):
    t = tensors[i]
    tensor = 0
    if name == "mechTensor":
      tensor = HillMandel2Voigt(to_mech_tensor(t.transpose()))
    if name == "elecTensor":
      tensor = to_elec_tensor(t.transpose())
    if name == "piezoTensor":
      tensor = to_piezo_tensor(t.transpose())

    angle, data, aux = perform_cfs_rotation(tensor, samples, aux_code)

    res_angle.append(angle)
    res_data.append(data if aux_code == "default" else aux)

  return res_angle, res_data  
    
# import pylab
# pylab.plot(data[:,0],data[:,1])
# pylab.show()  
 
#tensor = []
#t = to_mech_tensor(eval("[1.0,0.5,0.0,0.0,0.0,0.0]"))
#tensor.appt.end(to_mech_vector(t, as_array=True)) 
#centers = []
#centers.append([0.0,0.0,0.0])
#centers.append([1.0,1.0,0.0])

#angle, data = perform_rotations(tensor, 10)

#f = h5py.File("/home/fwein/project/simp/hook.h5")
#s1 = get_element(f, "design_stiff1_plain", "mech")
#s2 = get_element(f, "design_stiff2_plain", "mech") 
#angle = get_element(f, "design_rotAngle_plain", "mech")
 
#t = to_mech_tensor(eval("[2.607069, 1.484705, 0.1626158, 0,0, 0.3030707]")) 
 
#f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")
# f = h5py.File("/home/fwein/tmp/l_sl-m_20-g_al-p_0.1-gamma_1e-07-tau_0.01.h5")
#r  = centered_elements(f)
#tensor = get_element(f, "mechTensor", "piezo")
# tensor = get_element(f, "mechTensor", "mech")
#rot = perform_rotations(tensor, 10)
# rot = perform_rotations(tensor, t21=True)
# orientational_stiffness(r, rot, 1500, 2.5).show()
