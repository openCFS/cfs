#!/usr/bin/env python
from optimization_tools import *
from PIL import Image
try:
  from hdf5_tools import *
except:
  print("failed to import hdf5_tools, processing .h5 files won't work")    
import argparse
import sys

# work on density files, e.g. perform threshold

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="input density.xml file(s)")
parser.add_argument("--threshold", help="threshold for void material with input lower and 1", type=float)
parser.add_argument("--vol", help="threshold to match volume", type=float)
parser.add_argument("--lower", help="scale input range to new lower range", type=float)
parser.add_argument('--show', help="show output as image", action='store_true')
parser.add_argument('--attribute', help="what to read from input ('design' or 'physical')", default="design")
parser.add_argument('--set', help="specifiy set to read if you don't want the last one")
parser.add_argument('--output', help="optional output file name, for single input only")
args = parser.parse_args()

if args.output and len(args.input) > 1:
  print('error: --output only possible with single file input')
  sys.exit()

for input in args.input:
  if not os.path.exists(input):
    print('error: file not found: ' + input)
    sys.exit() 
  
  if input.endswith(".h5") or input.endswith(".h5ref") or input.endswith(".cfs"):
     f = h5py.File(input, 'r')
     dump_h5_meta(f)
     os.sys.exit()   
  # usually 'design' or 'physical'  
  dens = read_density(input, args.attribute, set=args.set)
  # this is 'design'
  des = dens if args.attribute == 'design' else read_density(input, 'design', set=args.set)
  
  nr = read_density(input, "nr", set=args.set) if dens.ndim == 1 else None
  
  if args.attribute != 'design':
    print("original 'design' min=" + str(numpy.amin(des)) + " max=" + str(numpy.amax(des)) + " vol=" + str(numpy.sum(des) / dens.size))
  print("original '" + args.attribute + "' min=" + str(numpy.amin(dens)) + " max=" + str(numpy.amax(dens)) + " vol=" + str(numpy.sum(dens) / dens.size))
  
  lower = args.lower if args.lower else numpy.amin(des)
  
  out = numpy.zeros(dens.shape)
  
  if args.threshold:
    out = threshold_filter(dens, args.threshold, lower, 1.0)    
  elif args.vol:
    out, _ = auto_threshold_filter(dens, lower, args.vol, 1.0)    
  elif args.lower:
    # is either an option for threshold or an action by itself 
    org_lower = numpy.amin(dens)
    assert(False)
    #print "scale data from lower " + str(org_lower) + " to lower " + str(args.lower)
  else:
    print("no action selected")
    sys.exit(0)    
  
  output = args.output if args.output and len(args.input) == 1 else input[:input.find('.density.xml')] + '-out.density.xml'
  
  print("write '" + output + "' min=" + str(numpy.amin(out)) + " max=" + str(numpy.amax(out)) + " vol=" + str(numpy.sum(out) / dens.size))
  if dens.ndim == 1:
    write_density_file(output, out, elemnr=nr)
  else:
    write_density_file(output, out)  
    
  if args.show and dens.ndim == 2:
    x, y = dens.shape
    ret = numpy.zeros((y, x), dtype="uint8")
    for i in range(y):
      for j in range(x):
        ret[y-i-1][j] = 255 - int(255 * out[j][i])
      
    img = Image.fromarray(ret)
    img = img.resize((800, y * 800/x))
    img.show()
    
  