#!/usr/bin/env python
import glob
import argparse
from cfs_utils import *

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="selection of the info.xml files to process (with wildcards), default is all")
parser.add_argument('--fast', help="don't generate images when they already exist", action='store_true')
parser.add_argument('--failsave', help="continue on errors", action='store_true')
args = parser.parse_args()

list = args.input if len(args.input) > 0 else glob.glob("*.info.xml")   

for f in list:
  problem = f[0:-len(".info.xml")]
  xml = open_xml(f)
   
  running = xpath(xml, "//cfsInfo/@status") == "running"
  label = problem if not running else problem + "_running"     

  if args.fast and not running and os.path.exists(label + ".2d.png") and os.path.exists(label + ".tiled.png") and os.path.exists(label + ".bloch.png"):
    continue
  
  try:
    execute("find_bloch_band_gaps.py " + problem + ".bloch.dat --gnuplot png --nicelabel > " + problem + ".plot")  
    execute("gnuplot -c " + problem + ".plot")
    execute("show_density.py " + problem + ".density.xml --save " + label + ".2d.png")
    execute("show_density.py " + problem + ".density.xml --tile 3 --save " + label + ".tiled.png")
  except:
    print('caught exception ' +  sys.exc_info()[0])
    if not args.failsafe:
      os.sys.exit(1)
    

