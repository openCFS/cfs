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
    step = len(ms) - 2 # reset to last, first element is ResultDescription
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
    
  
## visualize the orientational stiffness
# @return the image
def orientational_stiffness(centers, angle, data, nx, scale=-1):

  # zoom gives proper offsets of the elements
  min, max = find_corners(centers) # we assume the real grid to start at 0/0
  
  dim = (nx, int(nx *  (max[1] + min[1]) / (max[0] + min[0]))) 
  
  dx = dim[0] / (max[0] + 2.0 * min[0])
  dy = dim[1] / (max[1] + 2.0 * min[1])
  
  im = Image.new("RGB", dim, "white")
  draw = ImageDraw.Draw(im)

  max_val = numpy.max(data[:])
  min_val = numpy.min(data[:])
   
  if scale == -1:
    dist = 1.0 if len(centers) == 1 else centers[1][0] - centers[0][0] 
    scale = 0.48 * dx * dist / max_val
   
  sm = cmx.ScalarMappable(colors.Normalize(vmin=min_val, vmax=max_val), cmap=plt.get_cmap('jet'))
    
  for i in range(len(data)):
    coord = centers[i]
    x_off = (coord[0] + min[0]) * dx
    y_off = (coord[1] + min[1]) * dy
    m = numpy.max(data[i])
    c = sm.to_rgba(m)
    rgb = "rgb(" + str(int(255 * c[0])) + ", " + str(int(255*c[1])) + "," + str(int(255*c[2])) + ")"
    pol = to_polygons(angle[i], data[i], x_off, dim[1] - y_off, scale)
    draw.polygon(pol, fill=rgb, outline="black")

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
 
tensor = []
t = to_mech_tensor(eval("[1.0,0.5,0.0,0.0,0.0,0.0]"))
tensor.append(to_mech_vector(t, as_array=True)) 
centers = []
centers.append([0.0,0.0,0.0])
centers.append([1.0,1.0,0.0])

#angle, data = perform_rotations(tensor, 10)

#f = h5py.File("/home/fwein/project/simp/piezo_fmo.h5")
# f = h5py.File("/home/fwein/tmp/l_sl-m_20-g_al-p_0.1-gamma_1e-07-tau_0.01.h5")
#r  = centered_elements(f)
#tensor = get_element(f, "mechTensor", "piezo")
# tensor = get_element(f, "mechTensor", "mech")
#rot = perform_rotations(tensor, 10)
# rot = perform_rotations(tensor, t21=True)
# orientational_stiffness(r, rot, 1500, 2.5).show()
