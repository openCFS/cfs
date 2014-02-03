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


## this contains the physical fields stiff and angle in the macroscopic and fine discretization.
# The macroscopic is where we consider cells. This might be the original discretization
# The fine discretization is used for doing the streamlines. This might be the same as the macroscopic 
class Fields:
  ## this performs an interpolation, stores the data and allows access 
  class Data:
    def __init__(self, min, max, nx, input_space, input_data):
      
      self.nx = nx      
      dx = (max[0] - min[0]) / nx # == dy
      self.ny = int((max[1] - min[1] + 1e-3) / dx)
      
      out = numpy.zeros(nx*ny, 2)
      for y in range(ny):
        for x in range(nx):
          out[y * nx + x][0] = min[0] + dx * x + 0.5 * dx 
          out[y * nx + x][1] = min[0] + dx * y + 0.5 * dx
  
      self.data = scipy.interpolate.griddata(input_space, input_data, out, 'linear', -1.0)
  

  def __init__(self, coords, s1, s2, angle):
    
    centers, min, max, elem = coords
    # convert to 2D
    c, v = convert_interpolation_input(centers, s1, s2, angle)
    
    self.macro = Data(min, max, 6, c, v)
 
def show_streamline(coords, s1, s2, angle):            
  
  centers, min, max, elem = coords
  
  print len(centers)
  
  print min
  print max
  
  print elem
  
  fileds = Fields(coords, s1, s2, angle)
  
  
  
  fig = create_figure(min, max)
  
  v = [(0.2,0.2),(0.4,0.25),(0.5,0.5),(0.4,0.7)]
  
  path = Path(v)
  
  patch = matplotlib.patches.PathPatch(path, facecolor='none', lw=2)
  fig.add_patch(patch)
  #matplotlib.pyplot.show()
