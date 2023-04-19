#!/usr/bin/env python
import argparse
import os
import sys
import shutil

# this tool checks if we are within a cfs test case, takes the current directory as name and copies
# name.h5ref, name.ref.info.xml and if available name.ref.plot.dat and name.ref.density.xml
# @return True if copied
def copy(source, target):
  if not os.path.exists(source):
    print("cannot find '" + source + "' -> ignore ")
    return False
  else:
    print("copy '" + source + "' -> '" + target + "'")
    shutil.copy(source, target)
    return True
  

parser = argparse.ArgumentParser(description = 'Within a cfs test-case the reference files are created')
parser.add_argument('--dry', help="doesn't do anything but shows what would be done", action='store_true')
args = parser.parse_args()

pwd = os.getcwd()
if 'TESTSUIT' not in pwd:
  print("your current path containts no 'TESTSUIT', therefor this is probably not a cfs test case")
  os.sys.exit()

name = os.path.basename(pwd)

if not copy(os.path.join('results_hdf5', name + '.cfs'), name + '.h5ref'): 
  copy(os.path.join('results_hdf5', name + '.h5'), name + '.h5ref') # fallback to old stuff

copy(name + '.info.xml', name + '.ref.info.xml')

copy(name + '.plot.dat', name + '.ref.plot.dat')

copy(name + '.density.xml', name + '.ref.density.xml')

copy(name + '.bloch.dat', name + '.ref.bloch.dat')

copy('nonlin.txt', 'nonlin.ref.txt')
