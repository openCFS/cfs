#!/usr/bin/env python

from optimization_tools import *
from PIL import Image, ImageDraw, ImageColor
import sys
import argparse

# you can approximate the stuff iteratively in python by
# Image.fromarray(255* data.T).transpose(Image.FLIP_TOP_BOTTOM).show()
# with data being an array

def refine(vals, size):
  new = numpy.zeros((2*size, 2*size), dtype="uint8")
  for i in range(size):
    for j in range(size):
      val = vals[i,j]
      new[2*i,2*j] = val 
      new[2*i+1,2*j] = val 
      new[2*i,2*j+1] = val 
      new[2*i+1,2*j+1] = val
  return new


#@ return image, density_array
def density_to_image(filename, set, design):
  if not is_valid_density_file(filename):
    print "not a valid density file given!"
    sys.exit(1)

  dens = read_density(filename, attribute = 'design' if design else 'physical', set=set)
  
  x, y, z = getDim(dens)
  
  if z > 1:
    sys.exit("can handle only 2D data")

  # the image needs to be transposed, as the first index is y (row) and the second the column!

  ret = numpy.zeros((y, x), dtype="uint8")
  
  # copy data from linear list
  for i in range(y):
    for j in range(x):
        ret[y-i-1][j] = 255 - int(255 * dens[j][i])

  return Image.fromarray(ret), dens


def print_grid_on_image(I, dens):
  x,y = dens.shape
  s = x
  
  # does not work atm!!
  draw = ImageDraw.Draw(I)
  xsize =  I.size[0]
  ysize =  I.size[1]
  for i in range(s):
    f = 0
    if(s == 120):
      f = 8
    if(s == 240):
      f = 4
    iii = (i+1)*f
    draw.line((0, iii, xsize, iii), fill="Black")
    draw.line((iii, 0, iii, ysize), fill="Black")


parser = argparse.ArgumentParser()
parser.add_argument("input", help="the density.xml file to visualize")
parser.add_argument('--save', help="optional filename to write image")
parser.add_argument('--design', help="show 'design' instead of 'physical'", action='store_true')
parser.add_argument('--grid', help="draw mesh lines", action='store_true')
parser.add_argument('--orgsize', help="suppress resizing", action='store_true')
parser.add_argument('--info', help="print some info about the density file and exit", action='store_true')
parser.add_argument('--set', help="optional label of set, default is the last one")
args = parser.parse_args()

if args.info:
  ids = read_set_ids(args.input)
  print 'number of sets in ' + args.input + ': ' + str(len(ids)) 
  if len(ids) > 0:
    print "first set: '" + ids[0] + "'"
  if len(ids) > 1:
    print "last set: '" + ids[-1] + "'"
    
  os.sys.exit()  

img,dens = density_to_image(args.input, args.set, args.design)
img.convert('L')

if args.grid:
  print_grid_on_image(img,dens)

if not args.orgsize:
  ix, iy = dens.shape
  f = 1600 / max(ix, iy)
  img = img.resize((f * ix, f * iy))
#I = I.rotate(90)
#I = I.transpose(0)

if args.save:
  print "saving image to file " + args.save
  img.save(args.save)
else:
  img.show()

#I.save("tmp1.png")
#I.resize((480,480)).save("tmp2.png")
