#!/usr/bin/env python
import libxml2
import os.path
import glob
import argparse
from cfs_utils import *

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="selection of the info.xml files to process (with wildcards), default is all")
args = parser.parse_args()

list = args.input if len(args.input) > 0 else glob.glob("*.info.xml")   

for f in list:
  problem = f[0:-len(".info.xml")]
  doc = libxml2.parseFile(f)
  xml = doc.xpathNewContext()
   
  running = xpath(xml, "//cfsInfo/@status") == "running"
  label = problem if not running else problem + "_running"     
  
  execute("find_bloch_band_gaps.py " + problem + ".bloch.dat --gnuplot png --nicelabel > " + problem + ".plot")  
  execute("gnuplot -c " + problem + ".plot")
  execute("show_density.py " + problem + ".density.xml --save " + label + ".2d.png")
  execute("show_density.py " + problem + ".density.xml --tile 3 --save " + label + ".tiled.png")

