#!/usr/bin/env python
# this ia a tool for visualization of density.xml files with density and rotation design variables
import numpy as np
from numpy.linalg import norm
from numpy import sin
from numpy import cos
import os
import optimization_tools as ot
import cfs_utils as ut
from PIL import Image, ImageDraw
import sys
import argparse
import glob
import os

import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
#matplotlib.use('Agg')
#  matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
import matplotlib.patches as patches
import imageio

# minimal and maximal are vectors.
def create_figure(res, minimal, maximal):

  dpi_x = res / 100.0 * (maximal[0]-minimal[0])
  dpi_y = res / 100.0 * (maximal[1]-minimal[1])

  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_y))
  ax = fig.add_subplot(111)
  ax.set_xlim(min(0,minimal[0]), max(1,maximal[0]))
  ax.set_ylim(min(0,minimal[1]), max(1,maximal[1]))
  return fig, ax

#@ return image
def density_to_rotation_image(filename, set=None, color='gray', angles_only=False):
  if not ot.is_valid_density_file(filename):
    print("not a valid density file given!")
    sys.exit(1)

  # read density.xml which needs to contain design values for density and rotAngle
  xml = ut.open_xml(filename)
  if ut.has(xml, '//design[@name="density"]') and not angles_only:
    des, nr_vec = ot.read_multi_design(filename, 'density', 'rotAngle', matrix=True, set=set, attribute='physical')
    # extract density field and angle from matrix
    density = des[:,:,:,0]
    angle = des[:,:,:,1]
  # fallback if density does not exist, assume constant full material
  else:
    des, nr_vec = ot.read_multi_design(filename, 'rotAngle', matrix=True, set=set)
    # extract angle field from matrix
    angle = des[:,:,:,0]
    density = np.ones(angle.shape)
  x, y, z = ot.getDim(des)
  
  if z > 1:
    sys.exit("can handle only 2D data")

  # could respect non-unit regions and out of bounds movement
  minimal = [0,0]
  min_dim = min((x,y))
  maximal = [x/min_dim,y/min_dim] # normalize smaller dimension to 1, as there is no other element information in .density.xml
  
  fig, sub = create_figure(800, minimal, maximal)
  
  # element length, up to now assumed as uniform in x and y direction
  dx = max((1/x,1/y))

  max_val = 1
  min_val = 0
  sm = cmx.ScalarMappable(colors.Normalize(min_val, max_val), cmap=plt.get_cmap('gray' if color == 'grayscale' else color))

  for i in range(y):
    for j in range(x):
      # center of finite element
      center = [(j+.5)*dx, (i+.5)*dx]
      # we cut off density values larger than 1
      # these can occur in feature mapping when using certain functions for the smooth maximum
      rho = np.min((1.,density[j,i,0])) if not angles_only else 1
      theta = angle[j,i,0]
      # we plot this vector centered at the point "center" above
      vec = [.5*dx*np.cos(theta), .5*dx*np.sin(theta)]
      #c = sm.to_rgba(1-rho)
      X = [center[0]-vec[0], center[0]+vec[0]]
      Y = [center[1]-vec[1], center[1]+vec[1]]
      sub.add_line(plt.Line2D(X,Y,color=str(1-rho)))

  return (fig, sub)





parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="the density.xml file to visualize")
parser.add_argument('--save', help="optional filename to write image")
parser.add_argument('--angles_only', help="shows angles without density scaling", action='store_true')
parser.add_argument('--set', help="optional label of set, default is the last one")
parser.add_argument('--animgif', help="write single animated gif for all sets to filename")
parser.add_argument('--savesets', help="write <savesets>_XXX.png files for all sets")

args = parser.parse_args()

inp = args.input if len(args.input) != 1 else glob.glob(args.input[0]) # for Windows potentially globalize 

for file in args.input:
  img, den = density_to_rotation_image(file, args.set, color='gray', angles_only=args.angles_only)
  if args.save:
    print("saving image to file " + args.save)
    img.savefig(args.save)
  elif (args.animgif or args.savesets):
    # read sets: TODO - read XML only once!
    xml = ut.open_xml(inp[0])
    sets = []
    for set in xml.xpath('//set'):
      sets.append(int(set.attrib['id']))
      #print(etree.tostring(set.xpath('//[@id]')))
    print('read', len(sets),'sets from',inp[0] + ': ',end='') # only python3
    path = os.getcwd()
    dir = os.path.join(path, 'giffiles')
    if not os.path.exists(dir):
      os.mkdir(dir)
    for i in sets:
      print(i,' ',end='' if i < sets[-1] else '\n',flush=True)
      img, den = density_to_rotation_image(file, str(i), color='gray', angles_only=args.angles_only)
      img.savefig('giffiles/' + str(i).zfill(4) + '.png')
      plt.close(img)
    if args.animgif:
      # Build GIF
      with imageio.get_writer('mygif.gif', mode='I') as writer:
        for i in sets:
          image = imageio.imread('tmp/' + str(i).zfill(4) + '.png')
          writer.append_data(image)
  else:
   img.show()
   input("Press Enter to terminate.")


