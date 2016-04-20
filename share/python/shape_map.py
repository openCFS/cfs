#!/usr/bin/env python
# this ia a tool for shape map optimization 

import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
#matplotlib.use('Agg')
matplotlib.use('tkagg')

import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
from scipy import ndimage

from optimization_tools import *
import argparse


def create_figure(res):

  dpi_x = res / 100.0 
  
  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_x))
  ax = fig.add_subplot(111)
  ax.set_xlim(0.0, 1.0)
  ax.set_ylim(0.0, 1.0)
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print s   

def find_shape(shapes, el):
  for s in shapes:
    if el >= s.el[0] and el <= s.el[-1]:
      return s
  print 'none of the ' + str(len(shapes)) + ' shapes has an element with nr ' + str(el) 
  return None     
class Shape: 
  def __init__(self, id, dof):
    self.id = id
    self.dof = dof
    self.el = []
    self.val = []
    self.profile = []
    self.valid = [] # only necessary for image import
    if id == 0:
      self.color = 'b'  
    if id == 1:
      self.color = 'g'  
    if id == 2:
      self.color = 'r'  
    if id == 3:
      self.color = 'c'  
    if id == 4:
      self.color = 'm'  
    if id == 5:
      self.color = 'y'  
    if id == 6:
      self.color = 'k'
    assert(id <= 6)  
  
  def __str__(self):
    return "id=" + str(self.id) + " dof=" + str(self.dof) + " color=" + self.color + " #el=" + str(len(self.el)) + " el=" + str(self.el) 
  #@return x, y for center point a
  def get_center(self, idx):   
    free = idx * float(1./(len(self.el)-1))
    val  = self.val[idx]
    # this is the a middle line
    x = free if self.dof == 1 else val
    y = val  if self.dof == 1 else free
    return x, y
    
  # averages the current profiles (only valid)
  def average_valid_profile(self):
    assert(len(self.profile) == len(self.valid))   
    c = 0
    s = 0
    for i in range(len(self.profile)):
      if self.valid[i]:
         c += 1
         s += self.profile[i]
    
    assert(c > 0)
    return s / c
      
  # return the line coordinates for the profile
  #@param left (True) or right (False)
  def get_profile(self, idx_1, idx_2, left):
    x_val = []
    y_val = []

    x,y = self.get_center(idx_1)
    if self.dof == 0:
      x += (-1 if left else +1) * self.profile[idx_1]
    else:    
      y += (-1 if left else +1) * self.profile[idx_1]  
    
    x_val.append(x)
    y_val.append(y)
    
    x,y = self.get_center(idx_2)
    if self.dof == 0:
      x += (-1 if left else +1) * self.profile[idx_2]
    else:    
      y += (-1 if left else +1) * self.profile[idx_2]  
    
    x_val.append(x)
    y_val.append(y)

    return x_val, y_val    
    

# @param file withe grad.plot or density.xml
# @param profile use if not in file
def read_file(filename, profile):
  shapes = []

  if filename.endswith('.plot'):  
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
      if curr.id <> id:
        shapes.append(Shape(id = id, dof=int(line[3])))
        curr = shapes[-1]  
      curr.el.append(int(line[0]))
      curr.val.append(float(line[4]))
      curr.profile.append(profile)    
      
  else:
   nx, ny, nz = read_mesh_info(filename, silent=True)
   assert(nx == ny)  
      
   tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
   root = tree.getroot()
   query = '//set[last()]/shapeParamElement' 
   sett = root.xpath(query)

   curr = None
   for el in sett:
     nr = int(el.get('nr'))
     v =  float(el.get('design'))
     if nr == 0 or nr % (nx + 1) == 0:
       dof = 0 if el.get('dof') == 'x' else 1
       node = el.get('type') == 'node' or el.get('type') == None
       if node: 
         shapes.append(Shape(id = len(shapes), dof=dof))
         curr = shapes[-1]
       else:
         # we have a profile, hence search the corresponding shape
         if profile <> None:
           print 'error: profile data given in ' + filename + ' and concurrently via command line'
           sys.exit(-2)  
         curr = find_shape(shapes, nr - shapes[-1].el[-1])   
     
     if node:  
       curr.el.append(nr)
       curr.val.append(v)
     else:  
       curr.profile.append(v)  
  
  if len(shapes[-1].val) <> len(shapes[-1].profile):
    print "error: no profiles in '" + filename + "' and not given via command line"
    sys.exit(-3)    
  return shapes  

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
  if end <> -1:
    center = 1.0/len(line) * (start+end)/2.0
    profile = 1.0/len(line) * float(end-start)
    data.append((center, profile))

  return data        

# fills shape content for import_form_image based on the previous shapes data
def smart_append_data(n, data, shapes):
  assert(len(data) >= len(shapes))
  # easy case where there is no crossing from y to x or vice versa
  if len(data) == len(shapes):
    for i in range(len(data)):
       shape = shapes[i] 
       a, p = data[i] 
       # el is missing!
       shape.val.append(a)
       # validate profile, it might be too large
       if len(shape.profile) > 0 and p >= 1.2 * shape.average_valid_profile():
         shape.profile.append(shape.average_valid_profile())
         shape.valid.append(False)
       else:
         shape.profile.append(p)
         shape.valid.append(True)
  # we copy the last if we have not enough data
  if len(data) < len(shapes):
    for  s in shapes:
      s.val.append(s.val[-1])  
      s.profile.append(s.average_valid_profile())
      s.valid.append(False)
  # we search the best data if we have too much and throw it away if best is not good enough
  if len(data) > len(shapes):
    for shape in shapes:
      ref_a = shape.val[-1]  
      best = None  
      for dat in data:
        if best == None or abs(dat[0] - ref_a) < abs(best[0] - ref_a):
           best = dat  
      if abs(best[0] - ref_a) < 1.0/n:
        # use best found  
        #print 'best ' + str(best) + ' for ref ' + str(ref_a) + ' for line data ' + str(data)
        shape.val.append(best[0])  
        shape.profile.append(best[1])
        shape.valid.append(False)
      else:    
        # print 'best ' + str(best) + ' is not good enough for ref ' + str(ref_a) + ' for line data ' + str(data)
        # copy old stuff
        shape.val.append(shape.val[-1])  
        shape.profile.append(shape.average_valid_profile())
        shape.valid.append(False)

  
def import_from_image(filename):
  img = Image.open(args.input).transpose(Image.FLIP_TOP_BOTTOM)
  img = img.convert("L")
  dat = numpy.asarray(img, dtype="uint8")
  assert(len(dat.shape) == 2)
  nx, ny = dat.shape
  assert(nx == ny) # check what is all needed!

  shapes = []  
  # we assume that the first line determines the number of shapes
  # start with columns
  num = len(find_data_in_line(dat[:,0]))
  print "we assume " + str(num) + " horizontal shapes"
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=1))
  for i in range(nx):
     smart_append_data(nx, find_data_in_line(dat[:,i]), shapes[-num:]) # the last num shapes just appended

  num = len(find_data_in_line(dat[0,:]))
  print "we assume " + str(num) + " vertical shapes"
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=0))
  for i in range(nx):
     smart_append_data(nx, find_data_in_line(dat[i,:]), shapes[-num:]) # the last num shapes just appended
     
  cnt = 0   

  # add the element numbers   
  for shape in shapes:
    shape.el = range(cnt, cnt + len(shape.val))
    cnt += len(shape.val)
         
  return shapes
     
def plot_data(res, shapes):
  fig, sub = create_figure(res)
  
  for shape in shapes:
    n = len(shape.el)

    # circles
    for i in range(n):
      x,y = shape.get_center(i)
      c=plt.Circle((x,y), 0.005, color = shape.color)
      fig.gca().add_artist(c)
      
    for i in range(0,n-1):
      x, y = shape.get_profile(i, i+1, True) # left          
      l = plt.Line2D(x,y, color= shape.color)                                    
      sub.add_line(l)
      x, y = shape.get_profile(i, i+1, False) # right          
      l = plt.Line2D(x,y, color= shape.color)                                    
      sub.add_line(l)
      
  return fig, sub

# __name__ is 'shape_map' if imported or '__main__' if run as commandline tool
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml or .grad.plot file or a image")
  parser.add_argument("--profile", help="give the profile if it is not in input", type=float)
  parser.add_argument('--noshow', help="don't show the image", action='store_true')
  args = parser.parse_args()
  shapes = []
  if args.input.endswith('.xml') or args.input.endswith('.plot'):
    shapes = read_file(args.input, args.profile)
  else:
    shapes = import_from_image(args.input)
  fig, sub = plot_data(800, shapes)
  if not args.noshow:
    fig.show()
    raw_input("Press Enter to terminate.")
else:
  f = 'shape_map_mech_39.grad.plot'
  print f
  #shapes = read_file(f, profile=0.1)
  #dump_shapes(shapes)
  #fig, sub = plot_data(800, shapes)
  #fig.show()
