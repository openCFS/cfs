#!/usr/bin/env python
import argparse
import os
import sys

# Simply replace all dots but the last in a filename with an underscore.
# Makes valid latex figure names as latex interprets anything after the first dot as extension

parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='*', help="input files (wildcards ok) where all '.' but the last are replaced by '_'")
parser.add_argument('--dry', help="doesn't do anything but shows what would be done", action='store_true')
args = parser.parse_args()

for file in args.input:
  if not os.path.exists(file):
    print('error: file not found: ' + file)
    sys.exit() 
  
  ext = file.rfind('.')
  if ext > 1 and file[:ext].find('.') > 0:
    new = file[:ext].replace('.','_') + file[ext:]
    if args.dry:
      print(file + ' -> ' + new)
    else:
      print(new)
      os.replace(file, new)
    
  