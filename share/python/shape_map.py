#!/usr/bin/env python
# this ia a tool for shape map optimization 

import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
#matplotlib.use('Agg')
#matplotlib.use('tkagg')

import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
from scipy import ndimage
import numpy as np
from optimization_tools import *
import argparse
import copy
from scipy.interpolate import interp1d

# this is for TW to suppress annying 
#value_data = vtk.util.numpy_support.numpy_to_vtk(tmp_element_data_arr, deep=True)
#/usr/lib64/python3.6/site-packages/vtk/util/numpy_support.py:137: FutureWarning: Conversion of the second argument of issubdtype from `complex` to `np.complexfloating` is deprecated. In future, it will be treated as `np.complex128 == np.dtype(complex).type`.
#  assert not numpy.issubdtype(z.dtype, complex), \
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

try:
  from matviz_vtk import *
except:
  print('could not import matviz_vtk, hope you do not need it: ' + str(sys.exc_info()[0]))

# convert string to integer dof: x= 0 or int to string
def dof(val):
  if val == 0:
    return 'x'
  if val == 1:
    return 'y'
  if val == 2:
    return 'z'
  if val == 'x':
    return 0
  if val == 'y':
    return 1
  if val == 'z':
    return 2
  print(val)
  assert(False)

# 'x' and 'y' -> 'z' 
def flip_dof(dof1, dof2):
  
  if dof1 == 0:
    assert(dof2 == 1 or dof2 == 2)
    if dof2 == 1:
      return 2
    else:
      return 1

  if dof1 == 1:
    assert(dof2 == 0 or dof2 == 2)
    if dof2 == 2:
      return 0
    else:
      return 2

  if dof1 == 2:
    assert(dof2 == 0 or dof2 == 1)
    if dof2 == 0:
      return 1
    else:
      return 0
  
  if dof1 == 'x':
    assert(dof2 == 'y' or dof2 == 'z')
    if dof2 == 'y':
      return 'z'
    else:
      return 'y'

  if dof1 == 'y':
    assert(dof2 == 'x' or dof2 == 'z')
    if dof2 == 'z':
      return 'x'
    else:
      return 'z'

  if dof1 == 'z':
    assert(dof2 == 'x' or dof2 == 'y')
    if dof2 == 'x':
      return 'y'
    else:
      return 'x'

  assert(False)

# for 3D gives the third dof
def missing_dof(dof1, dof2):
  assert(dof1 == 0 or dof1 == 1 or dof1 == 2)
  assert(dof2 == 0 or dof2 == 1 or dof2 == 2)
  assert(dof1 != dof2)
  if dof1 == 0:
    return 1 if dof2 == 2 else 2
  if dof1 == 1:
    return 0 if dof2 == 2 else 2
  if dof1 == 2:
    return 1 if dof2 == 0 else 0
  assert(False)

def create_figure(res, minimal, maximal):

  dpi_x = res / 100.0 

  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_x))
  ax = fig.add_subplot(111)
  ax.set_xlim(min(0,minimal[0]), max(1,maximal[0]))
  ax.set_ylim(min(0,minimal[1]), max(1,maximal[1]))
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

def find_shape(shapes, el):
  for s in shapes:
    if el >= s.el[0] and el <= s.el[-1]:
      return s
  print('none of the ' + str(len(shapes)) + ' shapes has an element with nr ' + str(el)) 
  return None     

def find_shape_by_dof(shapes, dof):
  res = []  
  for s in shapes:
    if s.dof == dof:
      res.append(s)  
  return res     

# transforms ids from 0 to 6 to color codes 'b' to 'k'. Only for matplotlib, not for vtk!
def matplotlib_color_coder(id):
  id = id % 6
  if id == 0:
    return 'b'  
  if id == 1:
    return 'g'  
  if id == 2:
    return 'r'  
  if id == 3:
    return 'c'  
  if id == 4:
    return 'm'  
  if id == 5:
    return 'y'  
  if id == 6:
    return 'k'
  assert(id <= 6)  
  
  
class Shape: 
  def __init__(self, id, dof, dof_b = None, scale = [1.0,1.0,1.0], ref = None):
    self.id = id
    # legacy shape have no ref yet
    self.ref = ref
    self.dof = dof
    self.dof_b = dof_b

    self.scale = scale
    # the element number
    self.el = []
    # node variable a. For 3D there is also b
    self.a = []
    # 3D center nodes case
    self.b = []
    self.profile = []
    self.valid = [] # only necessary for image import
    # the color is in matplotlib code ('b', ...) transform for vtk
    self.color = matplotlib_color_coder(id)  
  
  # print the shape info
  def __str__(self):
    return "id=" + str(self.id) + " dof=" + str(self.dof) + " dof_b=" + str(self.dof_b) \
         + " color=" + self.color + " #el=" + str(len(self.el)) + " #b=" + str(len(self.b)) + " el=" + str(self.el)

  # shape info for a given index
  def to_string(self, idx):
    return "shape=" + str(self.id) + " dof=" + str(self.dof) + " color=" + str(self.color) + " idx=" + str(idx) + " val=" + str(self.a[idx]) + " profile=" + str(self.profile[idx]) + " valid=" + str(self.valid[idx]) + ""
   
  #@return x, y for center point a
  def get_center(self, idx):   
    free = idx * float(1./(len(self.el)-1))
    # fails for 3D!!!
    free *= self.scale[1-self.dof]
    val  = self.a[idx]
    # this is the a middle line
    x = free if self.dof == 1 else val
    y = val  if self.dof == 1 else free
    return x, y
    
  # averages the current profiles (only valid)
  def average_valid_profile(self):
    if len(self.profile) != len(self.valid):
      print(len(self.profile))
      print(len(self.valid))    
    assert(len(self.profile) == len(self.valid))   
    c = 0
    s = 0
    for i in range(len(self.profile)):
      if self.valid[i]:
         c += 1
         s += self.profile[i]
    
    assert(c > 0)
    return s / c

  # return the coordinates for both profile nodes
  #@return x1,y1,x2,y2
  def get_profiles(self, idx):
    x,y = self.get_center(idx)
    if self.dof == 0:
       return x - self.profile[idx], y, x + self.profile[idx], y
    else:    
      return x, y - self.profile[idx], x, y + self.profile[idx]   
      
  # return the line coordinates for the profile
  #@param left (True) or right (False)
  def get_profile(self, idx_1, idx_2, left):
    x_val = []
    y_val = []

    x,y = self.get_center(idx_1)
    if self.dof == 0:
      x += (-1.0 if left else +1.0) * self.profile[idx_1]
    else:    
      y += (-1.0 if left else +1.0) * self.profile[idx_1]  
    
    x_val.append(x)
    y_val.append(y)
    
    x,y = self.get_center(idx_2)
    if self.dof == 0:
      x += (-1.0 if left else +1.0) * self.profile[idx_2]
    else:    
      y += (-1.0 if left else +1.0) * self.profile[idx_2]  
    
    x_val.append(x)
    y_val.append(y)

    return x_val, y_val    
    
## symmetrize assues symmetry on x-axis (mirror)
def symmetrize(shapes):
  x = find_shape_by_dof(shapes, 0) 
  y = find_shape_by_dof(shapes, 1)
  assert(len(x) == 2)

  center = .5* (np.asarray(x[0].a) + np.asarray(x[1].a))
  offset = np.mean(center) - 0.5 # positive is shift to right
  shift = int(offset * len(y[0].a) + .5)

  print("center on x-axis is " + str(np.mean(center)) + " (shift " + str(shift) + " nodes) with var " + str(np.var(center))) 

  assert(np.var(center) < 0.001) 
  
  # avg and shift the left and right x-structures (from bottom to top)
  for i in range(len(x[0].a)):
    lv = x[0].a[i] - offset
    rv = x[1].a[i] - offset
    m = .5*(rv-lv)
    x[0].a[i] = .5 - m
    x[1].a[i] = .5 + m 
    
    lp = x[0].profile[i]
    rp = x[1].profile[i]
    mp = .5*(lp+rp)
    x[0].profile[i] = mp
    x[1].profile[i] = mp 

  # shift the horizontal lines and average left/right counterparts    
  for s in y:
    # we need the original data, otherwise rv might be replaced by shifted lv
    nodes = copy.deepcopy(s.a)
    profiles = copy.deepcopy(s.profile)

    for i in range(int(.5 * len(s.a))):  
      lv = nodes[i+shift]      
      rv = nodes[-(i+1-shift)]
      s.a[i] = s.a[-(i+1)] = .5*(lv+rv) 
      
      lp = profiles[i+shift]      
      rp = profiles[-(i+1-shift)]
      s.profile[i] = s.profile[-(i+1)] = .5*(lp+rp) 
     

# @param file withe grad.dat or density.xml. 2D or 3D!
# @param profile use if not in file
def read_file(filename, set, profile):
  shapes = None

  if filename.endswith('.plot'):
    shapes = read_plot_file(filename, profile)  
  else:
    xml = open_xml(filename)
    if has(xml, '//set[last()]/shapeParamElement/@ref'):
      shapes = read_xml(xml, set, profile)
    else:
      shapes = read_legacy_xml(filename, set, profile)

  # for both input formats
  if len(shapes[-1].a) != len(shapes[-1].profile):
    print("error: no profiles in '" + filename + "' and not given via command line")
    sys.exit(-3)    
  return shapes  

def read_plot_file(filename, profile):
  shapes = []
  all = open(filename).readlines()
  assert(len(all) > 2)
  comment = all[0].split('\t')
  # add first shape
  first = all[1].split('\t')
  shapes.append(Shape(id = int(first[2]), dof=int(first[3])))
  curr = shapes[0]
  for i in range(1,len(all)): # skip first comment
    line = all[i].split()
    id = int(line[2])
    if curr.id != id:
      shapes.append(Shape(id = id, dof=int(line[3])))
      curr = shapes[-1]  
    curr.el.append(int(line[0]))
    curr.a.append(float(line[4]))
    curr.profile.append(profile)
    curr.valid.append(True)    

  return shapes  

# this is the old xml parser before ref was added as attribute to .density.xml. Only 2D!
def read_legacy_xml(filename, set, profile):
  shapes = []
  # might be too much when not the whole domain is design (e.g. lbm with boundary)
  nx, ny, nz = read_mesh_info(filename, silent=True)
  assert(nx == ny)  
    
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  sq = 'last()' if not set else '@id="' + str(set) + '"'
  query = '//set[' + sq + ']/shapeParamElement' 
  sett = root.xpath(query)

  # check for full mesh and be tolerant against lbm meshes. 
  # we don't known if profiles are written
  nshapes = int(len(sett)/nx + .5)
  shape_elems = nx+1 # default case
  if nshapes > len(sett)/nx: 
    # apparently we have the lbm case
    print(nshapes, len(sett), nx, len(sett)/nx, len(sett)/nshapes, int(len(sett)/nshapes))
    # assert(len(sett)/nshapes == int(len(sett)/nshapes))
    shape_elems = int(len(sett)/nshapes)
    print('assume not the full domain is design: ' + str(nshapes) + ' shapes with ' + str(shape_elems) + ' (' + str(nx+1) + ') elements')   
    
  curr = None
  for el in sett:
    nr = int(el.get('nr'))
    v =  float(el.get('design'))
    if nr == 0 or nr % (shape_elems) == 0:
      d = dof(el.get('dof'))
      node = el.get('type') == 'node' or el.get('type') == None
      if node: 
        shapes.append(Shape(id = len(shapes), dof=d, ref = len(shapes)))
        curr = shapes[-1]
      else:
        # we have a profile, hence search the corresponding shape
        if profile != None:
          print('error: profile data given in ' + filename + ' and concurrently via command line')
          sys.exit(-2)  
        curr = find_shape(shapes, nr - shapes[-1].el[-1])   
   
    if node:  
      curr.el.append(nr)
      curr.a.append(v)
    else:  
      curr.profile.append(v)  
      curr.valid.append(True)
 
  # check if there was no profile in the xml
  if len(shapes[-1].a) > 0 and len(shapes[-1].profile) == 0:
     if profile != None: # error message comes below
       for shape in shapes:
         shape.profile = [profile] * len(shape.el)
         shape.valid   = [True] * len(shape.el)
  
  return shapes              

# reads 2D and 3D s
def read_xml(xml, set, profile):
  # find scaling assuming equal mesh sizing
  nx, ny, nz = read_mesh_info_xml(xml)
  scale = [1.0, float(ny)/nx, float(nz)/nx] # nz=1 ignored in 2D 
  
  shapes = []
  sq = 'last()' if set == None else '@id="' + str(set) + '"'
  ref = 0 # the current ref id to read the shape from, incremented at end of loop
  query = '//set[' + sq + ']/shapeParamElement[@ref="' + str(ref) + '"]'
  list = xml.xpath(query) 
  while list:
    # we do not know yet if we are 2D or 3D. For 3D center nodes, the there are two nodes dof the the shape dof is the third by definition
    first_dof = dof(list[0].get('dof')) # might change
    first_shape = int(list[0].get('shape'))
    shape = Shape(id = len(shapes), dof = first_dof, scale = scale, ref = ref)
    for idx, el in enumerate(list):
      nr = int(el.get('nr'))
      v  = float(el.get('design'))
      s  = int(el.get('shape'))
      t  = el.get('type')
      
      if s == first_shape:
        shape.a.append(v)
        shape.el.append(nr)
        # works currently only for 2D. This is the running coordinate
        #print(idx/len(list))
        
      else:
        if t == 'node':
          d  = dof(el.get('dof'))
          assert(d != first_dof)
          shape.dof_b = d # done everytime but who cares ...
          shape.b.append(v)
        else:
          shape.profile.append(v if not profile else profile) 
          shape.valid.append(True) # not really necessary    

    # print(len(shape.a), len(shape.b), len(shape.profile))
    assert(len(shape.a) == len(shape.profile))
    assert(len(shape.b) == 0 or len(shape.b) == len(shape.a)) # 3D case
    assert((len(shape.b) == 0 and shape.dof_b is None) or (len(shape.b) > 0 and shape.dof_b is not None))
    shapes.append(shape)
    
    ref += 1 if len(shape.b) == 0 else 2 
    list = xml.xpath('//set[' + sq + ']/shapeParamElement[@ref="' + str(ref) + '"]')

  return shapes   


# resamples given shapes. ignores valid
def resample(shapes, resample):
  res = []  
  org_space = np.linspace(0, 1.0, num=len(shapes[0].a), endpoint=True)
  new_space = np.linspace(0, 1.0, num=resample+1, endpoint=True)

  for o in shapes:
     s = Shape(o.id, o.dof, ref = o.ref)
     s.el = list(range(len(res) * (resample+1), (len(res)+1) * (resample+1)))
     v = interp1d(org_space, o.a, kind='cubic')
     s.a = v(new_space)
     p = interp1d(org_space, o.profile, kind='cubic')
     s.profile = p(new_space)
     s.valid = (resample+1) * [True]
     res.append(s)

  return res   


## create 2D data from 3D data. Removes all shapes with z-orientation.
# assumes values for dof_a and dof_b=z to be the same 
# in the sym case the forst shape is drwan four times in 2D
def downscale(shapes, mode):
  assert mode == 'straight' or mode == 'sym'
  assert len(shapes) > 0
  assert len(shapes[0].b) > 0
  
  res = []
  base_el = 0
  if mode == 'straight':
    for org in shapes:
      assert org.dof is not None and org.dof_b is not None
      if org.dof == 2 or org.dof_b == 2: # skip case with orientation z
        s = Shape(len(res), org.dof if org.dof != 2 else org.dof_b, ref=len(res))
        s.el = list(range(base_el, base_el + len(org.el)))
        base_el += len(org.el)
        s.a  = org.a if org.dof != 2 else org.b
        s.profile = org.profile
        res.append(s)
    if len(res) == 0:
      print('downscaling removed all ', len(shapes), ' shapes. Apparently all were with z-orientation.')
      sys.exit()    
  else:
    assert mode == 'sym'    
    assert len(shapes) == 12 # not really a necessary 
    org = shapes[0]
    n = len(org.el)
    
    s = Shape(0, dof=org.dof, ref=0)
    s.el = list(range(0, n))
    s.a  = org.a  
    s.profile = org.profile
    res.append(s)

    s = Shape(1, dof=org.dof, ref=1)
    s.el = list(range(n, 2*n))
    s.a  = 1-numpy.array(org.a)  
    s.profile = org.profile
    res.append(s)

    s = Shape(2, dof=1 if org.dof == 0 else 0, ref=2)
    s.el = list(range(2*n, 3*n))
    s.a  = org.a  
    s.profile = org.profile
    res.append(s)

    s = Shape(3, dof=1 if org.dof == 0 else 0, ref=3)
    s.el = list(range(3*n, 4*n))
    s.a  = 1-numpy.array(org.a)  
    s.profile = org.profile
    res.append(s)

  return res   

## create 3D data from 2D data. Meant for square symmetry
# this makes 3d center node shapes with second shape of dof=z and same values
def upscale(shapes, mode):
  assert mode == 'straight' or mode == 'sym'
  assert len(shapes) > 0
  assert len(shapes[0].b) == 0
  assert shapes[0].dof == 0 or shapes[0].dof == 1
  
  res = []
  
  if mode == 'straight':
    for org in shapes:
      s = Shape(len(res), org.dof, dof_b=2, ref=len(res))
      s.el = org.el
      s.a  = org.a
      s.b  = org.a # we copy the same information
      s.profile = org.profile
      #s.valid = org.valid
      res.append(s)
  else: # sym    
    assert(len(shapes) == 4)
    org = shapes[0]
    n = len(org.el)
    assert org.el == list(range(0,n))
    upscale_3d_sym_block(org, res, org.dof, 2)
    upscale_3d_sym_block(org, res, missing_dof(org.dof, 2), 2)
    upscale_3d_sym_block(org, res, org.dof, missing_dof(org.dof, 2))
  return res

def upscale_3d_sym_block(org, res, dof_a, dof_b):
  n = len(org.el)
  s = Shape(len(res), dof_a, dof_b=dof_b, ref=len(res))
  el_base = 0 if len(res) == 0 else res[-1].el[-1] + n + 1

  s.el = list(range(el_base,el_base+n))
  s.a  = org.a
  s.b  = org.a
  s.profile = org.profile
  res.append(s)
  
  s = Shape(len(res), dof_a, dof_b=dof_b, ref=len(res))
  s.el = list(range(el_base+2*n,el_base+3*n)) # first shape has elemes for a and b
  s.a  = 1-numpy.array(org.a)  
  s.b  = org.a
  s.profile = org.profile
  res.append(s)
    
  s = Shape(len(res), dof_a, dof_b=dof_b, ref=len(res))
  s.el = list(range(el_base+4*n,el_base+5*n)) # first shape has elemes for a and b
  s.a  = org.a  
  s.b  = 1-numpy.array(org.a)
  s.profile = org.profile
  res.append(s)
     
  s = Shape(len(res), dof_a, dof_b=dof_b, ref=len(res))
  s.el = list(range(el_base+6*n,el_base+7*n)) # first shape has elemes for a and b
  s.a  = 1-numpy.array(org.a)
  s.b  = 1-numpy.array(org.a)
  s.profile = org.profile
  res.append(s)

# searches for data (< 128) within a 1D line
# returns tupes with center and profile in meters (assume unit cube!)
# @param line vector with len == nx of content int 0 ... 255
def find_data_in_line(line):

  data = []

  # the token
  in_data = line[0] < 128
  start = 0 if in_data else -1
  end   = -1 # not found yet
  for i in range(1,len(line)):
    v = line[i]
    if v < 128 and not in_data: # data starts
      assert(start == -1)
      assert(end == -1)
      in_data = True
      start = i
      end   = -1
    elif v >= 128 and in_data: # data ends
      assert(start != -1)
      assert(end == -1)
      in_data = False
      end = i
      center = 1.0/len(line) * (start+end)/2.0
      profile = 1.0/len(line) * float(end-start)
      data.append((center, profile))   
      start = -1
      end = -1           
  # check if we ended with data
  if end != -1:
    center = 1.0/len(line) * (start+end)/2.0
    profile = 1.0/len(line) * float(end-start)
    data.append((center, profile))

  return data        
  
# fixes the shape data which is not valid by interpolating between the last valid
def repair_shapes(shapes):
  cnt = 0  
  for shape in shapes:
    start = -1
    end   = -1
    assert(len(shape.valid) == len(shape.a))
    # find invalid ranges
    in_valid = shape.a[0]
    assert(in_valid) # we shall start good
    for i in range(len(shape.valid)):
      if not shape.valid[i] and in_valid: # invalid starts
        assert(start == -1)        
        assert(end == -1)
        in_valid = False
        start = i
        end = -1
      elif shape.valid[i] and not in_valid: #invalid ends
        assert(start != -1)
        assert(end == -1)
        in_valid = True
        end = i-1
        #print 'last valid region: ' + shape.to_string(start-1)
        for j in range(end+1-start):
          #print start
          #print end
          #print end+2-start  
          #print 1.0/(end+2-start)
          bal = (j+1)*1.0/(end+2-start)
          val = shape.a[start-1] + bal * (shape.a[end+1] - shape.a[start-1])  
          #print 'ivalid: bal = ' + str(bal) + " val=" + str(val) + " -> " + shape.to_string(start+j)
          shape.a[start+j] = val
          shape.profile[start+j] = shape.average_valid_profile()
          cnt += 1
        #print 'invalid region: start ' + str(start) + " -> " + shape.to_string(start)
        #print 'invalid region: end ' + str(end) + " -> " + shape.to_string(end)
        #print 'first valid region: ' + shape.to_string(end+1)
        # reset
        start = -1
        end = -1      
  print("repaired " + str(cnt) + " times")        

# fills shape content for import_form_image based on the previous shapes data
def smart_append_data(n, data, shapes):
  assert(len(data) >= len(shapes))
  # easy case where there is no crossing from y to x or vice versa
  if len(data) == len(shapes):
    for i in range(len(data)):
       shape = shapes[i] 
       a, p = data[i] 
       # el is missing!
       shape.a.append(a)
       # validate profile, it might be too large
       if len(shape.profile) > 0 and p >= 1.3 * shape.average_valid_profile():
         shape.profile.append(shape.average_valid_profile())
         shape.valid.append(False)
       else:
         shape.profile.append(p)
         shape.valid.append(True)
  # we copy the last if we have not enough data
  if len(data) < len(shapes):
    for  s in shapes:
      s.a.append(s.a[-1])  
      s.profile.append(s.average_valid_profile())
      s.valid.append(False)
  # we search the best data if we have too much and throw it away if best is not good enough
  if len(data) > len(shapes):
    for shape in shapes:
      ref_a = shape.a[-1]  
      best = None  
      for dat in data:
        if best == None or abs(dat[0] - ref_a) < abs(best[0] - ref_a):
           best = dat  
      if abs(best[0] - ref_a) < 1.5/n:
        # use best found  
        # print 'best ' + str(best) + ' for ref ' + str(ref_a) + ' for line data ' + str(data)
        shape.a.append(best[0])  
        shape.profile.append(best[1])
        shape.valid.append(True)
      else:    
        #  print 'best ' + str(best) + ' is not good enough for ref ' + str(ref_a) + ' with err ' + str(abs(best[0] - ref_a)) + ' and criteria ' + str(1.5/n)
        # copy old stuff
        shape.a.append(shape.a[-1])  
        shape.profile.append(shape.average_valid_profile())
        shape.valid.append(False)

  
def import_from_image(filename, resample, repair):
  img = Image.open(args.input).transpose(Image.FLIP_TOP_BOTTOM)
  img = img.convert("L")
  dat = numpy.asarray(img, dtype="uint8")
  assert(len(dat.shape) == 2)
  nx, ny = dat.shape
  assert(nx == ny) # check what is all needed!
  if resample == None:
    resample = nx 
  shapes = []  
  # we assume that the first line determines the number of shapes
  # start with columns
  num = len(find_data_in_line(dat[:,0]))
  print("we assume " + str(num) + " horizontal shapes")
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=1))
  for i in numpy.linspace(0,nx-1, resample+1): #resample +1 because for 10 elements we have 11 variables
     smart_append_data(resample, find_data_in_line(dat[:,int(i)]), shapes[-num:]) # the last num shapes just appended

  num = len(find_data_in_line(dat[0,:]))
  print("we assume " + str(num) + " vertical shapes")
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=0))
  for i in numpy.linspace(0,nx-1, resample+1):
     smart_append_data(resample, find_data_in_line(dat[int(i),:]), shapes[-num:]) # the last num shapes just appended
     
  cnt = 0   

  # add the element numbers   
  for shape in shapes:
    shape.el = list(range(cnt, cnt + len(shape.a)))
    cnt += len(shape.a)
        
  if repair:       
    repair_shapes(shapes)       
  return shapes
    
# creates a matplotlib figure     
def plot_data(res, shapes, unit):
  # find extreme bounds to also visualize negative node positions
  minimal = [0.0]*2
  maximal = [1.0]*2
  if not unit:
    minimal = [1e9]*2
    maximal = [-1e9]*2
    for shape in shapes:
      minimal[shape.dof] = min(minimal[shape.dof], min(shape.a))
      maximal[shape.dof] = max(maximal[shape.dof], max(shape.a))
  
  fig, sub = create_figure(res, minimal, maximal)
  
  for shape in shapes:
    n = len(shape.el)

    # circles
    for i in range(n):
      x,y = shape.get_center(i)
      c=plt.Circle((x,y), 0.005, color = shape.color)
      fig.gca().add_artist(c)
      
    for i in range(0,n-1):
      x1, y1 = shape.get_profile(i, i+1, True) # left      
      l = plt.Line2D(x1,y1, marker='.', color=shape.color)        
      sub.add_line(l)
      x2, y2 = shape.get_profile(i, i+1, False) # right          
      l = plt.Line2D(x2,y2, marker='.', color=shape.color)                                    
      sub.add_line(l)
      
      if x1 == x2:
        sub.fill_between(x1,y1,y2,color=shape.color, alpha=.1, linewidth=0)
      else:
        sub.fill_betweenx(y1,x1,x2,color=shape.color, alpha=.1, linewidth=0)
  return fig, sub

# create vtk polydata tesselation
def create_2d_vtk(shapes):
  # create vtk cells and points
  points = vtk.vtkPoints()
  cells = vtk.vtkCellArray()
  
  for shape in shapes:
    # we initialize the loop with the -1 values
    cx,cy = shape.get_center(0)
    last_center = points.InsertNextPoint(cx, cy, 0.0)
  
    # last (-1) 'left' and 'right' profile nodes
    xl, yl, xr, yr = shape.get_profiles(0)
    last_left = points.InsertNextPoint(xl, yl, 0.0)
    last_right = points.InsertNextPoint(xr, yr, 0.0)
  
  
    for i in range(1,len(shape.el)):
      # this center node
      cx,cy = shape.get_center(i)
      this_center = points.InsertNextPoint(cx, cy, 0.0)
  
      # this 'left' and 'right' profile nodes
      xl, yl, xr, yr = shape.get_profiles(i)
      this_left = points.InsertNextPoint(xl, yl, 0.0)
      this_right = points.InsertNextPoint(xr, yr, 0.0)
    
      # two quadrialterals:
      # last_left--------this_left
      #    |                 |
      # last_center------this_center
      #    |                 |
      # last_right-------this_right
      # we divide the quadrilaterals by triangles
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_left)
      tri.GetPointIds().SetId(1, last_center)
      tri.GetPointIds().SetId(2, this_left)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_center)
      tri.GetPointIds().SetId(1, this_center)
      tri.GetPointIds().SetId(2, this_left)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_center)
      tri.GetPointIds().SetId(1, last_right)
      tri.GetPointIds().SetId(2, this_center)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_right)
      tri.GetPointIds().SetId(1, this_right)
      tri.GetPointIds().SetId(2, this_center)
      cells.InsertNextCell(tri)
  
      last_center = this_center
      last_left = this_left
      last_right = this_right
  
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

  return polydata


# helper for export
def ShapeParamElementXML(out, type, dof_val, id, ref, elnr, design):
   assert(type == 'node' or type == 'profile')
   assert(len(elnr) > 0)
   assert(len(design) == len(elnr))
   for i in range(len(elnr)):  
     out.write('    <shapeParamElement nr="' + str(elnr[i]) + '" type="' + type + '"') 
     if dof_val is not None:
       out.write(' dof="' + dof(dof_val) + '"')
     out.write(' shape="' + str(id) + '" ref="' + str(ref) + '"')
     out.write(' design="' + str(design[i]) + '"/>\n')
   
def export(shapes, filename, suppress_profile):
  d2 = len(shapes[0].b) == 0
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  nx = len(shapes[0].a) - 1 # todo! this can be done exactly!
  ny = nx # todo! this can be done exactly!
  nz = 1 if d2 else nx # todo! this can be done exactly!
  out.write('    <mesh x="' + str(nx) + '" y="' + str(ny) + '" z="' + str(nz) + '"/>\n')  
  out.write('  </header>\n')
  out.write('  <set id="shape_map.py">\n')
  # <shapeParamElement nr="0" type="node" dof="x" design="0.3"/>
  for shape in shapes:
    if d2:
      ShapeParamElementXML(out, 'node', shape.dof, shape.id, shape.ref, shape.el, shape.a)
    else:
      ShapeParamElementXML(out, 'node', shape.dof, 2*shape.id, 2*shape.ref, shape.el, shape.a)
      b_el = list(numpy.array(shape.el) + len(shape.el))
      ShapeParamElementXML(out, 'node', shape.dof_b, 2*shape.id+1, 2*shape.ref, b_el, shape.b)
  if not suppress_profile:
    el_base = shapes[-1].el[-1]+1 + ( 0 if d2 else len(shape.el)) # successive elements are already taken by the last b-shape   
    for shape in shapes:
      # also 3D has only one profile
      ref = shape.ref if d2 else 2*shape.ref
      id  = (1 if d2 else 2) * shapes[-1].id + shape.id + (1 if d2 else 2)
      ShapeParamElementXML(out, 'profile', None, id, ref, list(el_base + numpy.array(range(0, len(shape.el)))), shape.profile)
      el_base += len(shape.el)
  out.write('  </set>\n')
  out.write(' </cfsErsatzMaterial>\n')
  print("saved '" + filename + "'")


##see also sketch_tool.py  
def abaqus(shapes, filename):
  for shape in shapes:
    n = len(shape.a)
    assert len(shape.profile) == n
    
    for i in range(n):
      
    
    
      ShapeParamElementXML(out, 'node', shape.dof, shape.id, shape.ref, shape.el, shape.a)


# small permutation helper
#@param pm the indices define from where to take the value
def permutate(pm, v0, v1, v2):  
  assert(len(pm) == 3)
  v = [v0, v1, v2]
  #print('permutate: pm=',pm, '(',v0,v1,v2,') -> ',v[pm[0]], v[pm[1]], v[pm[2]]) 
  return (v[pm[0]], v[pm[1]], v[pm[2]])
  
  
## creates a 3d vtk polydate from a shape
#@param res how fine to resolve a circle, shall be even!
def create_3d_vtk_shape(shape, res):
  assert(len(shape.a) == len(shape.b) == len(shape.profile))
  assert(res % 2 == 0)

  # here we assume a cylinder standing on the xy plane, growing in the z-direction
  # this needs to permutated when we store the actual nodal coordinates
  pm = None
  if shape.dof == 0 and shape.dof_b == 1: # this is the assumed case xy
    pm = [0, 1, 2] 
  if shape.dof == 1 and shape.dof_b == 0: # yx plane, could be probably be same as above
    pm = [1, 0, 2] 
  if shape.dof == 0 and shape.dof_b == 2: # xz plane
    pm = [0, 2, 1]
  if shape.dof == 2 and shape.dof_b == 0: # zx plane
    pm = [1, 2, 0]
  if shape.dof == 1 and shape.dof_b == 2: # yz plane
    pm = [2, 0, 1]
  if shape.dof == 2 and shape.dof_b == 1: # zy plane
    pm = [2, 1, 0]
  
  # centers points
  ns = len(shape.a)
  z = np.linspace(0,1,ns)
  l = np.linspace(0,2*np.pi, res)

  points = []
  
  for s in range(ns-1):
    b = s*2*res # base
    c1x = shape.a[s] + shape.profile[s] * np.cos(l)
    c1y = shape.b[s] + shape.profile[s] * np.sin(l)
    c2x = shape.a[s+1] + shape.profile[s+1] * np.cos(l)
    c2y = shape.b[s+1] + shape.profile[s+1] * np.sin(l)
    
    for r in range(res):
      points.append(permutate(pm, c1x[r], c1y[r], z[s]))
      points.append(permutate(pm, c2x[r], c2y[r], z[s+1]))
      # points.append((c1x[r], c1y[r], z[s]))
      # points.append((c2x[r], c2y[r], z[s+1]))
  
  cells = []
  for s in range(ns-1):
    b = s*2*res # base
    for r in range(res):
      id1 = b + 2*r   # base line start
      id2 = b + 2*r+1 # upper line start
      id3 = (s+1)*2*res + 2*r # base line end
      id4 = (s+1)*2*res + 2*r+1 # upper line end
      cells.append((id1, id3, id2))
      cells.append((id3, id4, id2))
    
  cells = []
  for d in range(int(len(points)/4)):
    cells.append((4*d+0, 4*d+2, 4*d+1))  
    cells.append((4*d+2, 4*d+1, 4*d+3))
    if (d < int(len(points)/4)-1):
      cells.append((4*d+2, 4*d+4, 4*d+3))
      cells.append((4*d+4, 4*d+5, 4*d+3))
    
  poly = fill_vtk_polydata(points, cells)

  # setup the color scalar
  color = vtk.vtkUnsignedCharArray()
  color.SetNumberOfComponents(3)
  color.SetName("color")
  for c in range(len(cells)):
    color.InsertNextTuple3(0, 0, shape.id)  
  poly.GetCellData().SetScalars(color)
  
  return poly  

  
# __name__ is 'shape_map' if imported or '__main__' if run as commandline tool
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml or .grad.dat file or a image")
  parser.add_argument("--set", help="set nr within a .density.file", type=int)
  parser.add_argument("--profile", help="give the profile if it is not in input or to overwrite", type=float)
  parser.add_argument("--scale_profile", help="legacay files have a doubled profile value. Scale them by .5", type=float)
  parser.add_argument('--resample', help="resample to this resolution", type=int)
  parser.add_argument('--repair', help="interpolate unsure data (when parsing from image)", action='store_true')
  parser.add_argument('--symmetrize', help="mirror on x-axis", action='store_true')
  parser.add_argument('--transform', help="transforms 2D <-> 3D. 'sym' assumes square/cubic symmetry", choices=['straight','sym'])
  parser.add_argument('--export', help="write a density.xml file with shapeParam variables")
  parser.add_argument('--suppress_profile', help="do not export profile", action='store_true')
  parser.add_argument('--unbounded', help="do not restrict the visualization on a unit square", action='store_true')
  parser.add_argument('--save', help="save the image to the given name with the given format. Might be png, pdf, eps, vtp")
  parser.add_argument('--abaqus', help="write 2D data as abaqus sketch python script fragment to given filename")
  parser.add_argument('--noshow', help="don't show the image", action='store_true')
  parser.add_argument('--nooutline', help="don't show outline of the design domain for 3D", action='store_true')
  args = parser.parse_args()
  
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()
  
  shapes = []
  if args.input.endswith('.xml') or args.input.endswith('.xml.gz') or args.input.endswith('.plot'):
    shapes = read_file(args.input, args.set, args.profile)
    if args.resample:
      shapes = resample(shapes, args.resample)  
  else:
    shapes = import_from_image(args.input, args.resample, args.repair)
  
  print('average profile is ' + str(1.0/len(shapes) * sum([ s.average_valid_profile() for s in shapes])))
    
  if args.scale_profile:
    for s in shapes:
      s.profile = list(args.scale_profile * numpy.array(s.profile))  
    
  if args.symmetrize:
    symmetrize(shapes)  
  
  if args.transform:
    if len(shapes[0].b) == 0:
      shapes = upscale(shapes, args.transform)
    else:
      shapes = downscale(shapes, args.transform)
  
  if args.export:  
    export(shapes, args.export, args.suppress_profile)

  # do we do 3d? 
  d3 = len(shapes[0].b) > 0
  
  if args.abaqus:
    if d3:
      print('abaqus export only for 2D')
      os.sys.exit(1)
    abaqus(args.abaqus)   
  
  # vtp generation exclusivly triggerd by saving an vtp file
  if d3 or (args.save and args.save.endswith('.vtp')):
    poly = None
    if d3:
      ap = vtk.vtkAppendPolyData()
      for shape in shapes:
        poly = create_3d_vtk_shape(shape, 26)
        ap.AddInputData(poly)
      if not args.nooutline:
        ap.AddInputData(generate_outline_box([1,1,1], [0,0,0])) # make smarter if you need it  
      ap.Update()
      poly = ap.GetOutput()
    else:
      poly = create_2d_vtk(shapes)
    if args.save:
      show_write_vtk(poly,800,args.save)
    if not args.noshow:
      show_vtk(poly,800,show_edges=True if d3 else False) 
  else:
    fig, sub = plot_data(800, shapes, not args.unbounded)
    if args.save:
      print("write '" + args.save + "'")
      fig.savefig(args.save)
    if not args.noshow:
      fig.show()
      input("Press Enter to terminate.")
      
      
else:
  #f = 'shape_map_3d.density.xmp'
  #print(f)
  #shapes = read_file(f, None)
  #dump_shapes(shapes)
  #fig, sub = plot_data(800, shapes)
  #fig.show()
  import vtk
  from matviz_vtk import *
  poly = create_3d_vtk_shape()
  #show_write_vtk(poly, 200, 'polydata.vtp')
  show_vtk(poly, 800)
