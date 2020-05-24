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
parser.add_argument('--hist', help="show density as histogram", action='store_true')
parser.add_argument('--attribute', help="what to read from input ('design' or 'physical')", default="design")
parser.add_argument('--set', help="specifiy set to read if you don't want the last one")
parser.add_argument('--extrude', help="extrude 2d density to 3d. give number of slices", type=int)
parser.add_argument('--swap', help="swap dimensions of density field. useful for trelis element numbering 'xzy'", choices=['yx', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx'], type=str.lower)
parser.add_argument('--coarser', help="coarsen the density X times. useful in combination with extrude", type=int, default=0)
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
    print("original 'design' min=" + str(numpy.amin(des)) + " max=" + str(numpy.amax(des)) + " vol=" + str(numpy.sum(des) / des.size))
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
  elif args.extrude:
#    out = numpy.zeros((dens.shape[0], dens.shape[1], args.extrude))
    mx, my, _ = read_mesh_info(input, silent=True)
    print("original resolution: " + str(mx) + " x " + str(my))
    out = extrude(dens, args.extrude)
    print("new resolution: " + str(mx) + " x " + str(my) + " x " + str(args.extrude))
    if nr:
      nr = (nr + (mx*my) * numpy.arange(args.extrude)[:,None]).ravel()
  elif args.coarser:
    mx, my, _ = read_mesh_info(input, silent=True)
    print("original resolution: " + str(mx) + " x " + str(my))
  elif args.swap:
    print("swaping axes")
    out = des
  elif args.hist:
    import matplotlib.pyplot as plt
    plt.hist(dens.ravel(), bins = 10)
    plt.show()  
  else:
    print("no action selected")
    sys.exit(0)    
  
  output = args.output if args.output and len(args.input) == 1 else input[:input.find('.density.xml')] + '-out.density.xml'
  
  print("write '" + output + "' min=" + str(numpy.amin(out)) + " max=" + str(numpy.amax(out)) + " vol=" + str(numpy.sum(out) / out.size))
  if dens.ndim == 1:
    if args.swap:
      out = swap(out, args.swap)
    write_density_file(output, out, elemnr=nr)
  else:
    if args.swap:
      out = swap(out, args.swap)
    write_density_file(output, out)  
    for _ in range(args.coarser):
      coarsen_density(output, output)
      mx, my, mz = read_mesh_info(output, silent=True)
      print("new resolution: " + str(mx) + " x " + str(my) + " x " + str(mz))
    
  if args.show and dens.ndim == 2:
    x, y = dens.shape
    ret = numpy.zeros((y, x), dtype="uint8")
    for i in range(y):
      for j in range(x):
        ret[y-i-1][j] = 255 - int(255 * out[j][i])
      
    img = Image.fromarray(ret)
    img = img.resize((800, int(y * 800/x)))
    img.show()
    
  
