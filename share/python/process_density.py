#!/usr/bin/env python
from optimization_tools import *
from PIL import Image
import argparse
import sys

# work in density files, e.g. perform threshold

parser = argparse.ArgumentParser()
parser.add_argument("input", help="input density.xml file")
parser.add_argument("output", help="output density.xml file")
parser.add_argument("--threshold", help="threshold for void material with input lower and 1 (default 0.5)", type=float, default=0.5)
parser.add_argument("--vol", help="threshold to match volume", type=float)
parser.add_argument("--lower", help="scale input range to new lower range", type=float)
parser.add_argument('--show', help="show output as image", action='store_true')
parser.add_argument('--attribute', help="what to read from input ('design' or 'physical')", default="design")
parser.add_argument('--set', help="specifiy set to read if you don't want the last one")
parser.add_argument('--dim', help="dimension: 2D or 3D as integer",type=int)
args = parser.parse_args()


if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 
  
# usually 'design' or 'physical'  
dens = read_density(args.input, args.attribute, set=args.set)
# this is 'design'
des = dens if args.attribute == 'design' else read_density(args.input, 'design', set=args.set)
if args.dim == 3:
  x, y, z = dens.shape
else:
  x, y = dens.shape

if args.attribute <> 'design':
  print "for 'design' min=" + str(numpy.amin(des)) + " max=" + str(numpy.amax(des)) + " vol=" + str(numpy.sum(des) / (x*y*z if args.dim == 3 else x*y))
print "for '" + args.attribute + "' min=" + str(numpy.amin(dens)) + " max=" + str(numpy.amax(dens)) + " vol=" + str(numpy.sum(dens) / (x*y*z if args.dim == 3 else x*y))

lower = args.lower if args.lower else numpy.amin(des)

out = numpy.zeros(dens.shape)

if args.threshold:
  out = threshold_filter(dens, args.threshold, lower, 1.0)    
elif args.vol:
  out = auto_threshold_filter(dens, lower, args.vol, 1.0)  
elif args.lower:
  # is either an option for threshold or an action by itself 
  org_lower = numpy.amin(dens)
  assert(False)
  #print "scale data from lower " + str(org_lower) + " to lower " + str(args.lower)
else:
  print "no action selected"
  sys.exit(0)    

print "write '" + args.output + "' min=" + str(numpy.amin(out)) + " max=" + str(numpy.amax(out)) + " vol=" + str(numpy.sum(out) / (x*y))  
write_density_file(args.output, out)  
  
if(args.show) and args.dim != 3:
  ret = numpy.zeros((y, x), dtype="uint8")
  for i in range(y):
    for j in range(x):
      ret[y-i-1][j] = 255 - int(255 * out[j][i])
    
  img = Image.fromarray(ret)
  img = img.resize((800, dens.shape[1] * 800/dens.shape[0]))
  img.show()
  
  