#!/usr/bin/env python
from optimization_tools import *
from PIL import Image
try:
  from hdf5_tools import *
except:
  print("failed to import hdf5_tools, processing .h5 files won't work")    
import argparse
import sys
import numpy as np
import glob
 
parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='+', help="input density.xml file(s)")
parser.add_argument('--data', help="what to read from input (design or physical)", default="design")
parser.add_argument('--set', help="specifiy set to read if you don't want the last one")
parser.add_argument("--threshold", help="threshold for void material with input lower and 1", type=float)
parser.add_argument("--vol", help="threshold to match volume", type=float)
parser.add_argument("--lower", help="scale input range to new lower range", type=float)
parser.add_argument('--extrude', help="extrude 2d density to 3d. give number of slices", type=int)
parser.add_argument('--swap', help="swap dimensions of density field. useful for trelis element numbering 'xzy'", choices=['yx', 'xzy', 'yxz', 'yzx', 'zxy', 'zyx'], type=str.lower)
parser.add_argument('--coarser', help="coarsen the density X times. useful in combination with extrude", type=int, default=0)
parser.add_argument('--hist', help="show density as histogram with given bars", type=int)
parser.add_argument('--checkerboard', help="visualize ceckerboard calculation", action='store_true')
parser.add_argument('--show', help="show output as image", action='store_true')
parser.add_argument('--nosave', help="do not write output file", action='store_true')
parser.add_argument('--save', help="optional output file name, for single input only")
args = parser.parse_args()

allinput = args.input if len(args.input) != 1 else glob.glob(args.input[0]) # for Windows potentially globalize 

if args.save and len(allinput) > 1:
  print('error: --save only possible with single file input')
  sys.exit()

for input in allinput:
  if not os.path.exists(input):
    print('error: file not found: ' + input)
    sys.exit()
  
  if input.endswith(".h5") or input.endswith(".h5ref") or input.endswith(".cfs"):
     f = h5py.File(input, 'r')
     dump_h5_meta(f)
     os.sys.exit()   
     
  if not test_density_xml_attribute(input, args.data):
    print("Error: <element .../> attribute '" + args.data + "' not found in '" + input + "'")
    os.sys.exit()   
  
  # usually 'design' or 'physical'  
  data = read_density(input, args.data, set=args.set)
      
  if len(allinput) == 1:    
    # this is 'design'
    plain    = data if args.data == 'design' else read_density(input, 'design', set=args.set)
    physical = data if args.data == 'physical' else (None if not test_density_xml_attribute(input, 'physical') else read_density(input, 'physical', set=args.set))
    extract_all_density_info(plain, physical, silent = False)   
  
  nr = read_density(input, "nr", set=args.set) if data.ndim == 1 else None
  lower = args.lower if args.lower else np.amin(data)
  
  out = None

  if args.checkerboard:
    cb, img = checkerboard(data, image = True)
    print('normalized checkerboard for ' + args.data + ': ', cb)
    x, y = data.shape
    img = img.resize((800, int(y * 800/x)), Image.NEAREST)
    img.show()
    os.sys.exit()
  elif args.threshold:
    out = threshold_filter(data, args.threshold, lower, 1.0)    
  elif args.vol:
    out, _ = auto_threshold_filter(data, lower, args.vol, 1.0)  
  elif args.lower:
    # is either an option for threshold or an action by itself 
    org_lower = np.amin(data)
    assert(np.amax(data) <= 1.00000001)
    out = args.lower + (1-args.lower) * (data - org_lower) / (1-org_lower)
    print("scale data from lower " + str(org_lower) + " to lower " + str(args.lower))
  elif args.extrude:
#    out = np.zeros((data.shape[0], data.shape[1], args.extrude))
    mx, my, _ = read_mesh_info(input, silent=True)
    print("original resolution: " + str(mx) + " x " + str(my))
    out = extrude(data, args.extrude)
    print("new resolution: " + str(mx) + " x " + str(my) + " x " + str(args.extrude))
    if nr:
      nr = (nr + (mx*my) * np.arange(args.extrude)[:,None]).ravel()
  elif args.coarser:
    mx, my, _ = read_mesh_info(input, silent=True)
    print("original resolution: " + str(mx) + " x " + str(my))
  elif args.swap:
    print("swaping axes")
    out = data
  elif args.hist:
    import matplotlib.pyplot as plt
    plt.hist(data.ravel(), bins = args.hist)
    plt.show()
  elif args.show and data.ndim == 2:
    out = data
  if not out is None and not args.nosave:
    save = args.save if args.save and len(allinput) == 1 else input[:input.find('.density.xml')] + '-out.density.xml'
    
    print("write '" + save + "' min=" + str(np.amin(out)) + " max=" + str(np.amax(out)) + " vol=" + str(np.mean(out)))
    if data.ndim == 1:
      if args.swap:
        out = swap(out, args.swap)
      write_density_file(save, out, elemnr=nr)
    else:
      if args.swap:
        out = swap(out, args.swap)
      write_density_file(save, out)
      for _ in range(args.coarser):
        coarsen_density(save, save)
        mx, my, mz = read_mesh_info(save, silent=True)
        print("new resolution: " + str(mx) + " x " + str(my) + " x " + str(mz))

  if args.show and data.ndim == 2:
    x, y = data.shape
    ret = np.zeros((y, x), dtype="uint8")
    for i in range(y):
      for j in range(x):
        ret[y-i-1][j] = 255 - int(255 * out[j][i])
      
    img = Image.fromarray(ret)
    img = img.resize((800, int(y * 800/x)), Image.NEAREST)
    img.show()
    
