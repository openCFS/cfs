#!/usr/bin/env python

# this easy tool helps to set the CFS/OMP/MKL_NUM_THREADS environment variables for openCFS

import sys
import os
import argparse

def env(name):
  v = os.environ.get('CFS_NUM_THREADS')
  return '' if v is None else v

parser = argparse.ArgumentParser(description='help to set the CFS/OMP/MKL_NUM_THREADS environment variables for openCFS')
parser.add_argument('threads', nargs='*', help="optional number of threads to be set, if not fiven shows current")
args = parser.parse_args()


if len(args.threads) == 0:
  print("CFS_NUM_THREADS='" + env('CFS_NUM_THREADS') + "'")
  print("OMP_NUM_THREADS='" + env('OMP_NUM_THREADS') + "'")
  print("MKL_NUM_THREADS='" + env('MKL_NUM_THREADS') + "'")
elif len(args.threads) == 1:
  n = args.threads[0]
  assert int(n) > 0
  print('copy and paste:')
  print('export CFS_NUM_THREADS=' + str(n) + '; OMP_NUM_THREADS=' + str(n) + '; MKL_NUM_THREADS=' + str(n))
else:
  print('Usage: threads.py [<num threads> for CFS/OMP/MKL_NUM_THREADS]')
  sys.exit(1)


  
