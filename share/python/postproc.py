#!/usr/bin/env python
import glob
import sys
import argparse
import collections
from cfs_utils import *
from lxml.etree import LxmlSyntaxError


##@param key if not given the attribute @attribute from the query is the name
def add_key(xml, dic, query, key = None, quiet = True):
  value = xml.xpath(query)
  if not key:
    key = query.split('@')[1] if len(query.split('@')) == 2 else query
  if len(value) == 1:
    if value[0].isdigit(): # failes for -3
      dic[key] = int(value[0])
    elif isfloat(value[0]):
      dic[key] = float(value[0])
    else:  
      dic[key] = value[0]

def read_general_info(xml, dic):
  add_key(xml, dic, '//cfsInfo/@status')
  add_key(xml, dic, '//cfsInfo/summary/timer/@wall')
  add_key(xml, dic, '//cfsInfo/summary/timer/@cpu')
  add_key(xml, dic, '//cfsInfo/summary/memory/@peak', 'mem')
  add_key(xml, dic, '//cfsInfo/header/domain/@nx')
  add_key(xml, dic, '//cfsInfo/header/progOpts/@problem')
  add_key(xml, dic, '//cfsInfo/header/@id')

def read_selected_opt(xml, dic):
  last_iter = xml.xpath('//optimization/process/iteration')
  if len(last_iter) == 0:
    return
  add_key(last_iter[-1], dic, '@number', 'iter')
  on = xml.xpath('//optimization/header/objective/@name')
  if len(on) == 1:
    add_key(last_iter[-1], dic, '@' + on[0], on[0])
  add_key(last_iter[-1], dic, '@alpha')
  add_key(last_iter[-1], dic, '@slack')
  for s in range(1, 20):
    add_key(last_iter[-1], dic, '@bandgap_' + str(s) + '_' + str(s+1))
  
  tmp = xml.xpath('//constraints')
  assert(len(tmp) == 1)
  for c in tmp[0]:
    name = c.get('name').replace(' ','') # there was a bug with 'pyhsical_ volume'
    if '(' in name: # skip 'slope_(node)'
      continue
    if name:
      add_key(last_iter[-1], dic, '@' + name)

def read_all_opt(xml, dic):
  iter = xml.xpath('//optimization/process/iteration')
  if len(iter) == 0:
    return
  last = iter[-1]

  for a in last.attrib:
    add_key(last, dic, '@' + a)  

  
def handle_exception(args, problem, exception, message):
  if not args.silentfailsafe: 
    print('caught ' + exception + ' working with ' + problem + ': ' + message)
  if not args.failsafe and not args.silentfailsafe:
    os.sys.exit(1)

parser = argparse.ArgumentParser(description='General tool to analyze a bunch of .info.xml files')
parser.add_argument("input", nargs='*', help="selection of the info.xml files to process (with wildcards), default is all")
parser.add_argument("--query", help="xpath query, e.g. //transferFunction/@param where the attribute becomes the key")
parser.add_argument("--alliter", help="read all attributes from iteration", action='store_true')
parser.add_argument("--sort", help="sort for the key")
parser.add_argument("--revsort", help="sort reversly for the key")
parser.add_argument("--extract", help="extract only a single column. With sort, two columns are written")
parser.add_argument('--failsafe', help="continue on errors", action='store_true')
parser.add_argument('--silentfailsafe', help="continue on errors w/o message", action='store_true')
args = parser.parse_args()

input = args.input if len(args.input) > 0 else glob.glob("*.info.xml")   
# tuples name without .info.xml and larges relative band gap by ev_x_max and ev_(x+1)_min, the two values and the lower ev
res = []

do_sort = True if args.sort or args.revsort else False 
sort_key = None if not do_sort else (args.sort if args.sort else args.revsort)

for f in input:
  problem = f[0:-len(".info.xml")]
  try:
    xml = open_xml(f)
    dic = collections.OrderedDict()

    if args.query:
      add_key(xml, dic, args.query, quiet = False)
     
    if args.alliter: 
      read_all_opt(xml, dic)
    else:
      read_selected_opt(xml, dic)    
    read_general_info(xml, dic)
    res.append(dic)
  except KeyboardInterrupt:
    os.sys.exit(1)
  except RuntimeError as re:
   handle_exception(args, problem, 'RuntimeError', str(re)) 
  except lxml.etree.XMLSyntaxError as se:
   handle_exception(args, problem, 'XMLSyntaxError', str(se))    
#  except:
#   handle_exception(args, problem, 'exception', str(sys.exc_info()[0]))    

# all keys
keys = []
size = {}
for dic in res:
  for k in dic:
    if not k in keys:
      keys.append(k)
      size[k] = max(len(str(k)),len('(99)')) 
    size[k] = max(size[k], len(str(dic[k])))  

# enrich dict for missing keys such that we can sort anyway
for dic in res:
  for k in keys:
    if not k in dic:
      dic[k] = 0

if do_sort:
  res = sorted(res, key=lambda x: x[sort_key], reverse=True if args.revsort else False)

# print header for gnuplot output
if args.extract:
  if do_sort:
    print('#' + sort_key + ' ' + args.extract)
  else:
    print('#' + args.extract)
else:    
  if args.query:
    print('#query: ' + args.query)
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

# print result
for dic in res:
  if args.extract:
    if do_sort:
      print(str(dic[sort_key]) + " ", end='')
    print(str(dic[args.extract]))
  else:    
    print(' ',end='')
    for k in keys:
      s = '{:>' + str(size[k]) + '}'
      print(s.format(dic[k]) + ' ',end='')
    print()
    

