#!/usr/bin/env python
from optimization_tools import *
from PIL import Image
import argparse
import sys

# work in density files, e.g. perform threshold

parser = argparse.ArgumentParser()
parser.add_argument("input", help="input density.xml file")
parser.add_argument("output", help="output density.xml file")
parser.add_argument("--threshold", help="threshold for void material with input lower and 1 (default 0.5)", type=float)
parser.add_argument("--lower", help="scale input range to new lower range", type=float)
parser.add_argument('--show', help="show output as image", action='store_true')
parser.add_argument('--attribute', help="what to read from input ('design' or 'physical')", default="design")
parser.add_argument('--set', help="specifiy set to read if you don't want the last one")
args = parser.parse_args()


if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 
  
# usually 'design' or 'physical'  
dens = read_density(args.input, args.attribute, set=args.set)
# this is 'design'
des = dens if args.attribute == 'design' else read_density(args.input, 'design', set=args.set)

print "for '" + args.attribute + "' min=" + str(numpy.amin(dens)) + " max=" + str(numpy.amax(dens))

out = numpy.zeros(dens.shape)

if args.threshold:
    

if(args.show):
  x, y = dens.shape  
  ret = numpy.zeros((y, x), dtype="uint8")
  for i in range(y):
    for j in range(x):
      ret[y-i-1][j] = 255 - int(255 * dens[j][i])
    
  img = Image.fromarray(ret)
  img = img.resize((800, dens.shape[1] * 800/dens.shape[0]))
  img.show()
  
  