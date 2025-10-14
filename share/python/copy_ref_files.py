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
    return
  if not args.onlyclean and not args.dry:    
    print("copy '" + source + "' -> '" + target + "'")
    shutil.copy(source, target)
  else:
    print("would copy '" + source + "' -> '" + target + "' without " + ("--dry" if args.dry else "--onlyclean"))  

def clean(file):
  if os.path.isdir(file):
    if not os.listdir(file) == 0:
      if (args.dry is False and args.clean) or args.onlyclean:
        print("remove empty directory '" + file + "'")
        shutil.rmtree(file)    
      else:
        print("would remove empty directory '" + file + "' with --onlyclean or --clean and without --dry")  
    else:  
      print("directory '" + file + "' is not empty")
  else:
    if not os.path.exists(file):
      return
    if (args.dry is False and args.clean) or args.onlyclean:
      print("remove local file '" + file + "'")
      os.remove(file)    
    else:
      print("would remove local file '" + file + "' with  --onlyclean or --clean and without --dry")  
  

parser = argparse.ArgumentParser(description = 'Within a cfs test-case the reference files are created/overwritten')
parser.add_argument('--clean', help="remove known local files not mean for submission and write ref files", action='store_true')
parser.add_argument('--onlyclean', help="remove only known local files from local run", action='store_true')
parser.add_argument('--dry', help="doesn't do anything but shows what would be done", action='store_true')
args = parser.parse_args()

pwd = os.getcwd()
if 'TESTSUIT' not in pwd:
  print("your current path containts no 'TESTSUIT', therefor this is probably not a cfs test case")
  os.sys.exit()

name = os.path.basename(pwd)

# copy and clean respect --dry and --onlyclean
copy(os.path.join('results_hdf5', name + '.cfs'), name + '.h5ref')
copy(name + '.info.xml', name + '.ref.info.xml')
copy(name + '.plot.dat', name + '.ref.plot.dat')
copy(name + '.density.xml', name + '.ref.density.xml')
copy(name + '.bloch.dat', name + '.ref.bloch.dat')
copy(name + '.snopt', name + '.ref.snopt')
copy(name + '.mma', name + '.ref.mma')
copy('nonlin.txt', 'nonlin.ref.txt')

clean(os.path.join('results_hdf5', name + '.cfs'))
clean('results_hdf5')
clean(name + '.info.xml')
clean(name + '.plot.dat')
clean(name + '.density.xml')
clean(name + '.bloch.dat')
clean(name + '.snopt')
clean(name + '.mma')
