import sys
import h5py
import numpy
import Image, ImageDraw
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy import ndimage
import numpy as np

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
  

## give back elements with barycenters
# assumes rectangles
# @return list barycenter tuple ordered by elements and min and max node coordinates and first element dimensions
def centered_elements(hdf5_file):
  elements = hdf5_file['/Mesh/Elements/Connectivity']
  types = hdf5_file['/Mesh/Elements/Types']
  nodes = hdf5_file['/Mesh/Nodes/Coordinates']
  
  # first element dimensions
  node_coords = []
  for n in range(len(elements[0])):
    node_coords.append(nodes[elements[0][n]-1]) # numbers are one-based
 
  elem_dim = (max(node_coords[:][0]) - min(node_coords[:][0]), max(node_coords[:][1]) - min(node_coords[:][1]))
    
  result = []
  
  for e in range(len(elements)):
     if types[e] == 6:
       nod = elements[e]
       center = numpy.array([0.0, 0.0, 0.0])
       for n in range(len(nod)):
         center += nodes[nod[n]-1] # numbers are one-based
         # print "el=" + str(e) + " n=" + str(n) + " node=" + str(nod[n]) + "->" + str(nodes[nod[n]-1]) + " center=" + str(center) 
       center *= 1.0/len(nod)
       result.append(center)

  return result, (min(nodes[:,0]), min(nodes[:,1])), (max(nodes[:,0]), max(nodes[:,1])), elem_dim     
                
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
def to_rectangle(height, width, angle, x_offset, y_offset):
  
  # print "h=" + str(height) + " w=" + str(width) + " a=" + str(angle) + " x=" + str(x_offset) + " y=" + str(y_offset)
  
  height = numpy.max((height,1))
  
  tupl = []

  for x in [(-1.0,-1.0),(1.0,-1.0),(1.0,1.0),(-1.0,1.0)]:
    p = (x[0] * width/2, x[1] * height/2)
    r = (cos(angle) * p[0] -sin(angle)*p[1], sin(angle) * p[0] + cos(angle) * p[1])
    # r = (cos(angle) * p[0] + sin(angle)*p[1], - sin(angle) * p[0] + cos(angle) * p[1])
    tupl.append((x_offset + r[0], y_offset + r[1]))

  return tupl  

## give the corners to draw a rotated rectangles as polygon
def to_frustum(lower, upper):
  
  # print "to_frustum " + str(lower) + " -> " + str(upper)
  
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
  

## visualize the orientational stiffness
# @param grad is 'none' or 'linar'
# @return the image
def show_rot_frame(coords, s1, s2, angle, grad, nx, scale=-1):

  centers, min, max, elem = coords
  
  delta_angle = numpy.max(angle[:,0]) - numpy.max(angle[:,0]) 
   
  im, draw, dim, dx, dy = create_image_new(centers, min, max, nx,"white")

  height = elem[1] * dy 
  length = elem[0] * dx
 
  for i in range(len(s1)):
  
    coord = centers[i]
    
    #if coord[1] > 0.08:
    #  continue
          
    #x_off = (coord[0] + min[0]) * dx 
    #y_off = (coord[1] + min[1]) * dy

    x_off = (coord[0] + min[0] - 0.5 * elem[0]) * dx 
    y_off = (coord[1] + min[1]  - 0.5 * elem[1]) * dy

    #print str(coord) + " " + str(x_off) + " " + str(y_off)

    hor = s1[i,0]  # it seems that stiff1 and stiff2 are mixed up. This tries to correct it
    ver = s2[i,0]
    theta = angle[i,0]
    
    # lower horizontal line  
    pol = to_frustum((x_off, dim[1] - y_off), (x_off + length, dim[1] - y_off - height * ver))
    draw.polygon(pol, fill="black")

    # upper horizontal line
    pol = to_frustum((x_off, dim[1] - y_off - height + height * ver), (x_off + length, dim[1] - y_off - height))
    draw.polygon(pol, fill="black")

    # left vertical line
    pol = to_frustum((x_off, dim[1] - y_off), (x_off + length *hor, dim[1] - y_off - height))
    draw.polygon(pol, fill="black")

    # right vertical line
    pol = to_frustum((x_off + length - length * hor, dim[1] - y_off), (x_off + length, dim[1] - y_off - height))
    draw.polygon(pol, fill="black")


    # pol = to_rectangle(height * v1, length, theta, x_off, dim[1] - y_off + height - 2 * v1)
    # draw.polygon(pol, fill="black")
    
    
    #pol = to_rectangle(length, 0.25 * v2 * dy, theta, x_off, dim[1] - y_off)
    #draw.polygon(pol, fill="black")
        
    #pol = to_rectangle(length, length, theta + numpy.pi/2, x_off, dim[1] - y_off) 
    #draw.polygon(pol, fill="black")

    #pol = to_rectangle(0.9* length, 0.9 *length, theta, x_off, dim[1] - y_off) 

    #a+b
    #pol = to_rectangle(length * (1.-v2), (1.-v1)*length, theta + numpy.pi/2, x_off, dim[1] - y_off) 
    #draw.polygon(pol, fill="white")

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
