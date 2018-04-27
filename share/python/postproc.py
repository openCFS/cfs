#!/usr/bin/env python
import glob
import sys
import argparse
import collections
from cfs_utils import *

def add_key(xml, dic, query, key, quiet = True):
  value = xml.xpath(query)
  if len(value) == 1:
    dic[key] = float(value[0]) if value[0].isdigit() else value[0]

def read_general_info(xml, dic):
  add_key(xml, dic, '//cfsInfo/@status', 'status')
  add_key(xml, dic, '//cfsInfo/summary/timer/@wall', 'wall')
  add_key(xml, dic, '//cfsInfo/summary/timer/@cpu', 'cpu')
  add_key(xml, dic, '//cfsInfo/summary/memory/@peak', 'mem')
  add_key(xml, dic, '//cfsInfo/header/domain/@nx', 'nx')
  add_key(xml, dic, '//cfsInfo/header/progOpts/@problem', 'problem')
  
parser = argparse.ArgumentParser(description='General tool to analyze a bunch of .info.xml files')
parser.add_argument("input", nargs='*', help="selection of the info.xml files to process (with wildcards), default is all")
parser.add_argument("--sort", help="sort for the key")
parser.add_argument("--revsort", help="sort reversly for the key")
parser.add_argument('--failsafe', help="continue on errors", action='store_true')
args = parser.parse_args()

input = args.input if len(args.input) > 0 else glob.glob("*.info.xml")   
# tuples name without .info.xml and larges relative band gap by ev_x_max and ev_(x+1)_min, the two values and the lower ev
res = []

for f in input:
  problem = f[0:-len(".info.xml")]
  try:
    xml = open_xml(f)
    dic = collections.OrderedDict()
    
    read_general_info(xml, dic)
    res.append(dic)
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

# all keys
keys = []
size = {}
for dic in res:
  for k in dic:
    if not k in keys:
      keys.append(k)
      size[k] = max(len(str(k)),len('(99)')) 
    size[k] = max(size[k], len(str(dic[k])))  

if args.sort or args.revsort:
  res = sorted(res, key=lambda x: x[args.sort if args.sort else args.revsort], reverse=True if args.revsort else False)
# print header
print('#',end='')
for idx, k in enumerate(keys):
  s = '{:>' + str(size[k]) + '}'

  print(s.format('(' + str(idx+1) + ')') + ' ',end='')
print()
print('#',end='')
for k in keys:
  s = '{:>' + str(size[k]) + '}'
  print(s.format(k) + ' ',end='')
print()

# print res
for dic in res:
  print(' ',end='')
  for k in keys:
    s = '{:>' + str(size[k]) + '}'
    print(s.format(dic[k]) + ' ',end='')
  print()
    

