#!/usr/bin/env python
from optimization_tools import *
from PIL import Image, ImageDraw
import sys
import argparse
import glob

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
def density_to_image(filename, set, design, fillval=0.0):
  if not is_valid_density_file(filename):
    print("not a valid density file given!")
    sys.exit(1)


  if not design and not test_density_xml_attribute(filename, 'physical', set):
    print("the 'physical' design is not present, use non-physical 'design'")
    design = 'desgin'
  
  dens = read_density(filename, attribute = 'design' if design else 'physical', set=set, fill=fillval)

  x, y, z = getDim(dens)
  
  if z > 1:
    sys.exit("can handle only 2D data")

  # the image needs to be transposed, as the first index is y (row) and the second the column!

  ret = numpy.zeros((y, x), dtype="uint8")
  
  # copy data from linear list
  for i in range(y):
    for j in range(x):
      ret[y-i-1][j] = 255 - int(255 * min(dens[j][i],1))

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


# return image and density
def get_image(input, set, design, fill=0.0):
  img, dens = density_to_image(input, set, design, fill)
  img.convert('L')
  
  if args.grid:
    print_grid_on_image(img,dens)
  
  if args.orgsize:
    args.size = img.size
  else:
    if not args.size:
      args.size = [800, 800]

  print args.size
  img = img.resize(args.size)
  
  return img, dens


parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*',  help="the density.xml file to visualize or the files(s) for --saveall")
parser.add_argument('--save', help="optional filename to write image")
parser.add_argument('--saveall', help="saves all input files as png", action='store_true')
parser.add_argument('--design', help="show 'design' instead of 'physical'", action='store_true')
parser.add_argument('--grid', help="draw mesh lines", action='store_true')
parser.add_argument('--orgsize', help="suppress resizing", action='store_true')
parser.add_argument('--size', help="resize image", nargs='*', type=int)
parser.add_argument('--info', help="print some info about the density file and exit", action='store_true')
parser.add_argument('--set', help="optional label of set, default is the last one")
parser.add_argument('--tile', help="show periodic repetition of tile x tile patches", type=int)
parser.add_argument('--tileborder', help="show tile borders when repeating patches. works only in combination with --tile", action='store_true')
parser.add_argument('--fill', help="fill elements without density information with this pseudodensity value", type=float, default="0.0")

args = parser.parse_args()

if len(args.size) == 1:
    args.size = (args.size[0], args.size[0])
elif len(args.size) > 2:
    args.size = args.size[0:2]

input = args.input if len(args.input) > 0 else glob.glob("*.info.xml")
if not args.saveall:
  if not os.path.exists(input[0]):
    print("file '" + input[0] + "' not found")
    os.sys.exit()

if args.info:
  for file in input:
    ids = read_set_ids(file)
    print('number of sets in ' + file + ': ' + str(len(ids)))
    if len(ids) > 0:
      print("first set: '" + ids[0] + "'")
    if len(ids) > 1:
      print("last set: '" + ids[-1] + "'")
  
    print_design_info(file, 'design', args.set)
    print_design_info(file, 'physical', args.set)
  os.sys.exit()  

if args.saveall:
  for f in input:
    img, den = get_image(f, args.set, args.design)
    base = f[:-12] if f.endswith('.density.xml') else f
    print("save '" + base + ".png'")
    img.save(base + '.png')
else:
  for file in input:
    img, den = get_image(file, args.set, args.design, args.fill)

    if args.tile:
      assert(img.size[0] == img.size[1]) # extend if you need
      if args.orgsize:
        args.size = (img.size[0]*args.tile, img.size[1]*args.tile)
      else:
        if not args.size:
          args.size = [1000, 1000]
      img = img.resize((int(args.size[0]/args.tile), int(args.size[1]/args.tile)))
      nx, ny = img.size
      dat = numpy.array(img) 
      tiled = numpy.zeros((args.tile * ny, args.tile * nx), dtype="uint8")
      for i in range(args.tile):
        for j in range(args.tile):
          tiled[i*nx : (i+1)*nx , j*ny : (j+1)*ny] = dat
          if args.tileborder:
              if i > 0:
                tiled[i*nx,:] = 0
              if j > 0:  
                tiled[:,j*nx] = 0
      img = Image.fromarray(tiled)
    if args.save:
      print("saving image to file " + args.save)
      img.save(args.save)
    else:
     img.show()
              
