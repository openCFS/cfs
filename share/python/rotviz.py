#!/usr/bin/env python
# this ia a tool for visualization of density.xml files with density and rotation design variables
import numpy as np
from numpy.linalg import norm
from numpy import sin
from numpy import cos
import os
import optimization_tools as ot
import cfs_utils as ut
from spaghetti import combine_angles
from PIL import Image, ImageDraw
import sys
import argparse
import glob
import os
from lxml import etree

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
def density_to_rotation_image(filename, set=None, color='gray', angles_only=False, fill=None, coarsen=False):
  if not ot.is_valid_density_file(filename):
    print("not a valid density file given!")
    sys.exit(1)

  # read density.xml which needs to contain design values for density and rotAngle
  xml = ut.open_xml(filename)
  if ut.has(xml, '//design[@name="density"]') and not angles_only:
    des, nr_vec = ot.read_multi_design(filename, 'density', 'rotAngle', matrix=False, set=set, attribute='physical')
    # extract density field and angle from matrix
    density = des[:,0]
    angle = des[:,1]
  # fallback if density does not exist, assume constant full material
  else:
    des, nr_vec = ot.read_multi_design(filename, 'rotAngle', matrix=False, set=set)
    # extract angle field from matrix
    angle = des[:,0]
    density = np.ones(angle.shape)

  # get dimensions
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  x, y, z = ot.read_mesh_info_xml(tree)
  print(x,y,z)

  if z > 1:
    sys.exit("can handle only 2D data")

  # could respect non-unit regions and out of bounds movement
  minimal = [0,0]
  min_dim = min((x,y))
  maximal = [x/min_dim,y/min_dim] # normalize smaller dimension to 1, as there is no other element information in .density.xml
  
  fig, sub = create_figure(1600, minimal, maximal)
  
  # element length, up to now assumed as uniform in x and y direction
  dx = max((1/x,1/y))

  max_val = 1
  min_val = 0
  sm = cmx.ScalarMappable(colors.Normalize(min_val, max_val), cmap=plt.get_cmap('gray' if color == 'grayscale' else color))

  count = 0
  rows = y
  cols = x
  if coarsen:
    assert(int(x/2) == x/2)
    assert(int(y/2) == y/2)
    rows = int(y/2)
    cols = int(x/2)
    dx *= 2
  for i in range(rows):
    for j in range(cols):

      if count < len(angle):
        if not coarsen:
          # we cut off density values larger than 1
          # these can occur in feature mapping when using certain functions for the smooth maximum
          rho = np.min((1.,density[count])) if not angles_only else 1
          theta = angle[count]
        else:
          if (2*i+1)*x+2*j >= len(angle):
            if fill is None:
              return (fig, sub)
            else:
              rho = 1
              theta = float(fill)
          else:
            rho = np.min((1.,(density[2*i*x+2*j] + density[(2*i+1)*x+2*j] + density[2*i*x+2*j+1] + density[(2*i+1)*x+2*j+1]) / 4)) if not angles_only else 1
            theta = combine_angles([angle[2*i*x+2*j],angle[(2*i+1)*x+2*j],angle[2*i*x+2*j+1],angle[(2*i+1)*x+2*j+1]],[density[2*i*x+2*j], density[(2*i+1)*x+2*j], density[2*i*x+2*j+1], density[(2*i+1)*x+2*j+1]])
          count += 3
      else:
        if fill is None:
          return (fig, sub)
        rho = 1
        theta = float(fill)

      if rho > 5e-2:
        # center of finite element
        center = [(j+.5)*dx, (i+.5)*dx]
        # we plot this vector centered at the point "center" above
        vec = [.3*dx*np.cos(theta), .3*dx*np.sin(theta)]
        #c = sm.to_rgba(1-rho)
        X = [center[0]-vec[0], center[0]+vec[0]]
        Y = [center[1]-vec[1], center[1]+vec[1]]
        sub.add_line(plt.Line2D(X,Y,color=str(1-rho)))
      count += 1

  return (fig, sub)





parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="the density.xml file to visualize")
parser.add_argument('--save', help="optional filename to write image")
parser.add_argument('--angles_only', help="shows angles without density scaling", action='store_true')
parser.add_argument('--set', help="optional label of set, default is the last one")
parser.add_argument('--animgif', help="write single animated gif for all sets to filename")
parser.add_argument('--savesets', help="write <savesets>_XXX.png files for all sets")
parser.add_argument('--fill', help="Fill incomplete data with default angle and full material. Default: leave empty", default=None)
parser.add_argument('--noticks', help="omit axis tick labels", action='store_true')
parser.add_argument('--coarsen', help="coarsen grid resolution by factor two", action='store_true')


args = parser.parse_args()

inp = args.input if len(args.input) != 1 else glob.glob(args.input[0]) # for Windows potentially globalize 

for file in args.input:
  img, den = density_to_rotation_image(file, args.set, color='gray', angles_only=args.angles_only, fill=args.fill, coarsen=args.coarsen)

  if args.noticks:
    plt.tick_params(
        axis='both',          # changes apply to x-axis and y-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelleft=False,   # labels along the left edge are off
        labelbottom=False) # labels along the bottom edge are off
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


