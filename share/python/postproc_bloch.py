#!/usr/bin/env python
import glob
import sys
import argparse
from cfs_utils import *
from concurrent.futures._base import _create_and_install_waiters

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="selection of the info.xml files to process (with wildcards), default is all")
parser.add_argument('--fast', help="don't generate images when they already exist", action='store_true')
parser.add_argument('--failsafe', help="continue on errors", action='store_true')
parser.add_argument('--sort', help="don't generate images but print results sorted by normalized band-gap", action='store_true')
args = parser.parse_args()

input = args.input if len(args.input) > 0 else glob.glob("*.info.xml")   
# tuples name without .info.xml and larges relative band gap by ev_x_max and ev_(x+1)_min, the two values and the lower ev
results = []

for f in input:
  problem = f[0:-len(".info.xml")]
  try:
    xml = open_xml(f)
     
    running = xpath(xml, "//cfsInfo/@status") == "running"
    label = problem if not running else problem + "_running"     
  
    if args.sort:
      iter = xml.xpath("//iteration[last()]")[0]

      for ev in range(1,30):
        if 'ev_'+ str(ev) + '_max' in iter.attrib and 'ev_'+ str(ev+1) + '_min' in iter.attrib:  
          lower = float(iter.attrib['ev_'+ str(ev) + '_max'])
          upper = float(iter.attrib['ev_'+ str(ev+1) + '_min'])
          assert(lower > 0)
          results.append((problem, (upper-lower) / (lower + .5*(upper-lower)), lower, upper, ev))
          break # no need to loop further
          
    else:        
      if args.fast and not running and os.path.exists(label + ".2d.png") and os.path.exists(label + ".tiled.png") and os.path.exists(label + ".bloch.png"):
        continue
      
      execute("find_bloch_band_gaps.py " + problem + ".bloch.dat --gnuplot png --nicelabel > " + problem + ".plot")  
      execute("gnuplot -c " + problem + ".plot")
      execute("show_density.py " + problem + ".density.xml --save " + label + ".2d.png")
      execute("show_density.py " + problem + ".density.xml --tile 3 --save " + label + ".tiled.png")
  except KeyboardInterrupt:
    os.sys.exit(1)
  except RuntimeError as re:
   print('caught RuntimeError working with ' + problem + ': ' + str(re))
   if not args.failsafe:
     os.sys.exit(1)
  except:
    print('caught exception working with ' + problem + ': ' + str(sys.exc_info()[0]))
    if not args.failsafe:
      os.sys.exit(1)
    
if args.sort:
  sorted = sorted(results, key=lambda x: x[1], reverse=False)
  for gap in sorted:
    lower = gap[2]
    upper = gap[3]
    rel = (upper-lower) / lower
    print('norm: {:4.2f}'.format(gap[1]) + ' rel: {:5.2f}'.format(rel) + ' {:7.1f}'.format(gap[2]) + ' ->{:7.1f}'.format(gap[3]) + ' ' + str(gap[4]) + '->' + str(gap[4]+1) + ' \t' + gap[0])
