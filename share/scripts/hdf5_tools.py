import sys
import h5py
import numpy
import Image, ImageDraw
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx



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
# @return list barycenter tuple ordered by elements
def centered_elements(hdf5_file):
  elements = hdf5_file['/Mesh/Elements/Connectivity']
  types = hdf5_file['/Mesh/Elements/Types']
  nodes = hdf5_file['/Mesh/Nodes/Coordinates']
  
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

  return result     
                
## find mininmal and maximal coordiante
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
    tupl.append((x_offset + scale * x,y_offset + scale * y))
  return tupl    

## give the corners to draw a rotated rectangles as polygon
def to_rectangle(height, width, angle, x_offset, y_offset):
  
  # print "h=" + str(height) + " w=" + str(width) + " a=" + str(angle)
  
  height = numpy.max((height,1))
  
  tupl = []

  for x in [(-1.0,-1.0),(1.0,-1.0),(1.0,1.0),(-1.0,1.0)]:
    p = (x[0] * width/2, x[1] * height/2)
    r = (cos(angle) * p[0] -sin(angle)*p[1], sin(angle) * p[0] + cos(angle) * p[1])
    # r = (cos(angle) * p[0] + sin(angle)*p[1], - sin(angle) * p[0] + cos(angle) * p[1])
    tupl.append((x_offset + r[0], y_offset + r[1]))

  return tupl  
  
  
def create_image(centers, nx):
  # zoom gives proper offsets of the elements
  min, max = find_corners(centers) # we assume the real grid to start at 0/0
  
  dim = (nx, int(nx *  (max[1] + min[1]) / (max[0] + min[0]))) 
  
  dx = dim[0] / (max[0] + 2.0 * min[0])
  dy = dim[1] / (max[1] + 2.0 * min[1])
  
  im = Image.new("RGB", dim, "white")
  draw = ImageDraw.Draw(im)
  
  return im, draw, dim, dx, dy, min, max
  
## visualize the orientational stiffness
# @return the image
def show_rot_rect(centers, s1, s2, angle, show, nx, scale=-1):

  im, draw, dim, dx, dy, min, max = create_image(centers, nx)

  length = 0.8 * (centers[1][0] - centers[0][0]) * dx
  
  sm = cmx.ScalarMappable(colors.Normalize(vmin=0.0, vmax=0.5), cmap=plt.get_cmap('jet'))
  
  for i in range(len(s1)):
  
    coord = centers[i]
    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy

    v1 = s1[i,0]
    v2 = s2[i,0]
    theta = angle[i,0]
    
    #v1 = 0.5
    #v2 = 0.01
    #theta = 2.5
    #v1 = 0.15
    #v2 = 0.5
    #theta = 0.2

    

    # b
    size = length * v2 if show == "thickness" else 2
    pol = to_rectangle(size, length, theta + numpy.pi/2, x_off, dim[1] - y_off) 
    draw.polygon(pol, fill="black" if show == "thickness" else color_code(sm, v2))

    # a
    size = length * v1 if show == "thickness" else 2
    pol = to_rectangle(size, length, theta, x_off, dim[1] - y_off) 
    draw.polygon(pol, fill="black" if show == "thickness" else color_code(sm, v1))

  return im  

def color_code(color_map, value):
  c = color_map.to_rgba(value)
  return "rgb(" + str(int(255 * c[0])) + ", " + str(int(255*c[1])) + "," + str(int(255*c[2])) + ")"


# draws a thick circle, where the thickness is determined automatically by the radius 
def draw_thick_circle(draw, center, radius):
  
  for o in numpy.arange(0.0, 0.01 * radius, 0.05):
    draw.ellipse((center[0] - radius + o, center[1] - radius + o, center[0] + radius - o, center[1] + radius - o), fill=None, outline="black")

  
## visualize the orientational stiffness
# @return the image
def orientational_stiffness(centers, angle, data, nx, scale=-1):

  max_val = numpy.max(data[:])
  min_val = numpy.min(data[:])
  data_offset = -2 * min_val if min_val < 0 else 0

  im, draw, dim, dx, dy, min, max = create_image(centers, nx)
   
  print "max=" + str(max_val) + " min=" + str(min_val)
   
  if scale == -1:
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
def perform_rotations(tensors, notation, samples, name = "mechTensor", aux_code = "default"):

  res_angle = []
  res_data  = []
  for i in range(len(tensors)):
    t = tensors[i]
    #print "pr: t=" + str(t)
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
