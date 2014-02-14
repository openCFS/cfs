from matviz_2d import *
import matplotlib.pyplot
from matplotlib.path import Path
from scipy import ndimage
import numpy 
import scipy.interpolate


## create and prepare a matplot figure where patched might be "added" to
def create_figure(min, max):
  fig = matplotlib.pyplot.figure()
  ax = fig.add_subplot(111)
  ax.set_xlim(min[0],max[0])
  ax.set_ylim(min[1],max[1])
  
  return ax


## Subclass for Fields
# this performs an interpolation, stores the data and allows access 
class Data:
  def __init__(self, min, max, nx, input_space, input_data, method):
    
    self.NO_VAL = -12345.6
    
    self.min = min
    self.max = max
    self.nx = nx      
    dx = (max[0] - min[0]) / nx 
    self.dx = dx
    ny = int((max[1] - min[1] + 1e-3) / dx)
    self.ny = ny
    # dy might be slightly different from dx but is close
    dy = (max[1] - min[1]) / ny
    self.dy = dy

    out = numpy.zeros((nx*ny, 2))
    for y in range(ny):
      for x in range(nx):
        out[y * nx + x][0] = min[0] + dx * x + 0.5 * dx 
        out[y * nx + x][1] = min[0] + dy * y + 0.5 * dy

    self.data = scipy.interpolate.griddata(input_space, input_data, out, method, self.NO_VAL)
  
    if method <> 'nearest':
      self.nearest = scipy.interpolate.griddata(input_space, input_data, out, 'nearest')



  ## tests if the coordinates are within min, max  
  def within(self, x, y): 
    if x < self.min[0] or x > self.max[0]:
      return False
    if y < self.min[1] or y > self.max[1]:
      return False
    return True
    
  ## gives for the coordinates the value line. Does not further interpolate but assumes piecewise constant data.
  # has a fallback to nearest    
  def getData(self, x, y):  
    i = int((x - self.min[0]) / self.dx)
    #j = int((self.max[1] - y - self.min[1]) / self.dx)
    j = int((y - self.min[1]) / self.dy)
    #print 'getData x=' + str(x) + ' y=' + str(y) + ' -> i=' + str(i) + ' j=' + str(j) + ' idx=' + str(self.nx * j + i),
    v = self.data[self.nx * j + i]
    if v[0] == self.NO_VAL:
      v = self.nearest[self.nx * j + i]
    
    #print ' -> ' + str(v) 
    return v
    
  def dump_data(self):
    for j in range(self.ny):
      for i in range(self.nx):
        x = self.min[0] + i * self.dx + 0.5 * self.dx
        y = self.min[1] + j * self.dy + 0.5 * self.dy
        o = self.getData(x, y)
        print 'i=' + str(i) + ' j=' + str(j) + ' -> ' + str((x, y)) + ' -> ' + str(o) 
        
  ## dumps the data as image for debug purpose
  def dump_image(self, data_idx):
    out = numpy.zeros((self.ny, self.nx))

    minv = min(self.data[:,data_idx])
    maxv = max(self.data[:,data_idx])
    
    print minv
    print maxv
    
    for j in range(self.ny):
      for i in range(self.nx):
        #v = self.data[self.nx*j + i][data_idx]
        o = self.getData(i * self.dx, j * self.dy)[data_idx]
        s = (o + -1.0 * minv) / (maxv - minv) * 255
        #print ' scale from ' + str(o) + ' -> ' + str(s) + ' base=' + str(o + -1.0 * minv) + ' scale=' + str((maxv - minv) * 255)
        out[j,i] = s
        
    img = Image.fromarray(out)
    #img.resize((600,600), Image.ANTIALIAS)
    img.show()      

## this contains the physical fields stiff and angle in the macroscopic and fine discretization.
# The macroscopic is where we consider cells. This might be the original discretization
# The fine discretization is used for doing the streamlines. This might be the same as the macroscopic 
class Fields:

  ## constructs the data sets macro and fine
  def __init__(self, coords, s1, s2, angle, macro_samples):
    
    centers, min, max, elem = coords
    # convert to 2D
    c, v = convert_interpolation_input(centers, s1, s2, angle)
    
    nx = int((max[0] - min[0]) / elem[0])

    self.macro = Data(min, max, nx if macro_samples == None else macro_samples, c, v, 'linear')
    self.fine  = Data(min, max, 4 * nx, c, v, 'linear')
 
 
  ## calculates the streamline from a given point up to the end
  # @param idx 0 for s1 and 1 for s2 which changes the angle!
  # return a list of coordinates. The first entry is the starting point
  def streamline(self, x, y, idx, steplength):
    
    assert(idx == 0 or idx == 1)
    
    trace = []
    trace.append(((x, y), self.fine.getData(x, y)[idx]))
    
    while True:
      here = self.fine.getData(x, y)
      angle = here[2] + idx * numpy.pi/2 # for idx==0 nothing changes
      
      #print 'x=' + str(x) + ' y=' + str(y) + ' data=' + str(here),
      x += steplength * numpy.cos(angle)
      y += steplength * numpy.sin(angle)
      #print ' next_x=' + str(x) + ' next_y=' + str(y)
      
      if self.fine.within(x, y) and len(trace) < 2000: # prevent circular streams
        trace.append(((x, y), here[idx]))
      else:
        break

    return trace  
      

## draws a trace  
# @param scale scale stiffness to color. 1.0 for normal or 2.0 if max stiff = 0.5
def draw_trace(fig, trace, scale):
  # we search for portions in the trace with similar value which need to be at least two segments wide.
  # Small value errors are ok. value equals color
  
  vertices = [i[0] for i in trace]
  values   = [i[1] for i in trace]

  start = 0
  pos = 1
  val = values[0]
  while pos < len(trace):
    # up to the value changes or we reach the end
    if abs(values[pos] - val) > 0.1 or pos == len(trace)-1:
      assert(pos > start)
      #print 's=' + str(start) + ' e=' + str(pos)
      path = Path(vertices[start:pos])
      patch = matplotlib.patches.PathPatch(path, edgecolor=str(1.0 - (val * scale)), facecolor='none', lw=1)
      fig.add_patch(patch)

      start = pos
      val = values[pos]      
    pos += 1
      
## sort data such that s1 > s2. Will change angle by +/ pi/2. Sorts in-place!
def sort_data(s1, s2, angle):
  
  for i in range(len(s1)):
    if s1[i] < s2[i]:
      #print ' before: s1=' + str(s1[i]) + ' s2=' + str(s2[i]) + ' a=' + str(angle[i]),  
      t = float(s1[i]) # without float() the reference is copied and t changes value!
      s1[i] = float(s2[i])
      s2[i] = t
      angle[i] = float(angle[i]) - numpy.pi/2 if angle[i] > numpy.pi/2 else float(angle[i]) + numpy.pi/2 
      #print '  after: s1=' + str(s1[i]) + ' s2=' + str(s2[i]) + ' a=' + str(angle[i]) 
 
def show_streamline(coords, s1, s2, angle, dir, samples):            
 
  centers, min, max, elem = coords
  
  #print len(centers)
  #print min
  #print max
  #print elem
  
  sort_data(s1, s2, angle)

  fields = Fields(coords, s1, s2, angle, samples)
  #fields.macro.dump(s1,0)
  #fields.macro.dump_data()
  
  fig = create_figure(min, max)

  #trace = fields.streamline(0.251, 0.251, 0.1)
  #draw_trace(fig, trace)

  #interpret horizontal as s1. Note that we sort!!
  dirs = [1,0] # 1 is weaker, draw first and second the darker s1 lines
  if dir == 'horizontal':
    dirs = [0]
  if dir == 'vertical':
    dirs = [1]

  if True:
    for idx in dirs:
      for j in range(fields.macro.ny):
        for i in range(fields.macro.nx):
         trace = fields.streamline(i * fields.macro.dx, j * fields.macro.dy, idx, 0.02)
         draw_trace(fig, trace, 2.0)
  
  
  matplotlib.pyplot.show()
