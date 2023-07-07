#!/usr/bin/env python
import glob
import sys
import argparse
import collections
import re
from optimization_tools import *
from lxml.etree import LxmlSyntaxError


# try to convert a value to int or float, if it fails, leave it but make empty string to 0
def interpret_value(val):
  if val.isdigit(): # failes for -3
    return int(val)
  if isfloat(val):
    return float(val)
  if val.is_text:
    return val.__str__().split() # converts a vector or matrix into an array
  if val == '':
    return 0
  return val

##@param key if not given the attribute @attribute from the query is the name
def add_key(xml, dic, query, key = None, quiet = True):
  value = xml.xpath(query)
  if not key:
    key = query
    if query.count('@') == 1:
      key = query.split('@')[1]
    if query.count('@') > 1:
      key = label(query)
  if len(value) >= 1:
    val = interpret_value(value[0])
    if not isinstance(val, (list, tuple, np.ndarray)):
      dic[key] = val
    else:
      for idx, elem in enumerate(val):
        dic[key + '_' + str(idx)] = elem
    

def read_general_info(xml, dic):
  add_key(xml, dic, '//cfsInfo/@status')
  add_key(xml, dic, '//cfsInfo/summary/timer/@wall')
  add_key(xml, dic, '//cfsInfo/summary/timer/@cpu')
  add_key(xml, dic, '//cfsInfo/summary/memory/@peak', 'mem')
  add_key(xml, dic, '//cfsInfo/header/domain/@nx')
  add_key(xml, dic, '//cfsInfo/header/progOpts/@problem')
  split_problem(dic, dic['problem']) # id shall be last  
  add_key(xml, dic, '//cfsInfo/header/@id')

def read_perf(xml, dic, extend):
  add_key(xml, dic, '//cfsInfo/summary/timer/@wall')
  add_key(xml, dic, '//cfsInfo/summary/timer/@cpu')
  add_key(xml, dic, '//openmp/@CFS_NUM_THREADS','cfs')
  add_key(xml, dic, '//openmp/@OMP_NUM_THREADS','omp')
  if extend:
    add_key(xml, dic, '//cfsInfo/summary/memory/@peak', 'mem')
    add_key(xml, dic, '//grid/@elements') 
    add_key(xml, dic, '//feFunctions/@totalNumEqns') # fails for more than one PDE! 
  add_key(xml, dic, '//cfsInfo/header/progOpts/@problem')


def read_opt_issue(xml, dic):
  add_key(xml, dic, '//optimization/summary/@problem', 'opt_issue')

def read_selected_opt(xml, dic, first=False):
  iter = xml.xpath('//optimization/process/iteration')
  if len(iter) == 0:
    return

  idx = 0 if first else -1
  
  add_key(iter[idx], dic, '@number', 'iter')
  on = xml.xpath('//optimization/header/objective/@name')
  if len(on) == 1:
    add_key(iter[idx], dic, '@' + on[0], on[0])
  add_key(iter[idx], dic, '@alpha')
  add_key(iter[idx], dic, '@slack')
  for s in range(1, 20):
    add_key(iter[idx], dic, '@bandgap_' + str(s) + '_' + str(s+1))
  
  tmp = xml.xpath('//constraints')
  assert(len(tmp) == 1)
  for c in tmp[0]:
    name = c.get('name').replace(' ','') # there was a bug with 'pyhsical_ volume'
    if '(' in name: # skip 'slope_(node)'
      continue
    if name:
      add_key(iter[idx], dic, '@' + name)

def read_all_opt(xml, dic, first=False):
  iter = xml.xpath('//optimization/process/iteration')
  if len(iter) == 0:
    return
  last = iter[0 if first else -1]

  for a in last.attrib:
    add_key(last, dic, '@' + a)  

  
def handle_exception(args, problem, exception, message):
  if not args.silentfailsafe: 
    print('caught ' + exception + ' working with ' + problem + ': ' + message)
  if not args.failsafe and not args.silentfailsafe:
    os.sys.exit(1)
    
# read a bunch of .info.xml files    
def read_info_xml_input(args, input): 
  res = []   
  for f in input:
    problem = f[0:-len(".info.xml")]
    if not f.endswith(".info.xml"):
      handle_exception(args, problem, 'wrong filetype', '')
    try:
      xml = open_xml(f)
      dic = collections.OrderedDict()

      if args.skipaborted:
        status = xml.xpath('//cfsInfo/@status')
        if len(status) == 1 and status[0] == 'aborted':
          continue

      if args.query:
        for query in args.query:
          add_key(xml, dic, query, quiet = False)
        
      if args.perf:
        read_perf(xml, dic, args.extend)
      else:   
        if args.extend:
          read_all_opt(xml, dic, args.first)
        else:
          read_selected_opt(xml, dic, args.first)
        if args.issue:
          read_opt_issue(xml, dic)
        read_general_info(xml, dic)

      res.append(dic)
    except KeyboardInterrupt:
      os.sys.exit(1)
    except KeyError as re:
      handle_exception(args, problem, 'KeyError', str(re))
    except RuntimeError as re:
      handle_exception(args, problem, 'RuntimeError', str(re))
    except lxml.etree.XMLSyntaxError as se:
      handle_exception(args, problem, 'XMLSyntaxError', str(se))
  #  except:
  #   handle_exception(args, problem, 'exception', str(sys.exc_info()[0]))    
  return res  

def read_density_xml_input(args, input): 
  res = []   
  for f in input:
    problem = f[0:-len(".density.xml")]
    try:
      plain    = None if not test_density_xml_attribute(f, 'design') else read_density(f, 'design')
      physical = None if not test_density_xml_attribute(f, 'physical') else read_density(f, 'physical')
      
      dic = extract_all_density_info(plain, physical, silent = True)
      dic['problem'] = problem
      split_problem(dic, dic['problem']) # id shall be last
         
      res.append(dic)
    except KeyboardInterrupt:
      os.sys.exit(1)
    except RuntimeError as re:
      handle_exception(args, problem, 'RuntimeError', str(re)) 
  return res  


# read .dat files, e.g. plot.dat or grad.dat
# using plotviz might not be the fastest but who cares :)
def read_dat_input(args, input):
  import plotviz
  res = []   
  for f in input:
    dic = collections.OrderedDict()
    meta, data = plotviz.process(f)
    
    assert len(meta) == len(data[0])  
    for i, m in enumerate(meta):
      dic[m] = data[0 if args.first else -1][i] # already no more string if possible
    problem = f
    if problem.endswith('.plot.dat'):
      problem = problem[:-len('.plot.dat')]  
    if problem.endswith('.grad.dat'):
      problem = problem[:-len('.grad.dat')]  
    if problem.endswith('.dat'):
      problem = problem[:-len('.dat')]  
      
    dic['problem'] = problem
    split_problem(dic, problem) 
    res.append(dic)

  return res    

# tries to extract numerical values and their keys from the filename 
# filter-o_snopt-f_9.5 -> f:9.5
def split_problem(dic, name):
  # https://stackoverflow.com/questions/45001775/find-all-floats-or-ints-in-a-given-string/45001796
  name = str(name)
  ll = re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", name)
  sub = name
  for n in ll:
     pos = sub.find(n)
     key = sub[:pos]
     #print(key,key.split('-'),key.split('_'),[b.split('_') for b in key.split('-')],[b.split('-') for b in key.split('_')])
     #[b.split('-') for b in key.split('_')]
     token = [] # becomes ['penal', 'o', 'snopt', 'p', '']
     for s in key.split('_'):
       token += s.split('-')
     for t in reversed(token):
       if len(t) > 0:
         key = t
         break # first occurence is fine
    
     dic[key] = interpret_value(n)
     sub =  sub[pos+len(n):]

# tries to make a lable out of a xpath query. Partially copy & paste from study.py
def label(query):
  # clean the stuff as much as possible
  for s in ['/','[', ']', '@', '"', '=', 'item', 'result', 'data', 'value', 'type', 'name', 'param', 'constraint']:
    query = query.replace(s,'')     

  return query  



parser = argparse.ArgumentParser(description='General tool to analyze a bunch of .info.xml files')
parser.add_argument("input", nargs='+', help="selection of .info.xml, .density.xml or .plot.dat files to process (with wildcards)")
parser.add_argument("--query", nargs='+', help="xpath queries, e.g. //transferFunction/@param where the attribute becomes the key")
parser.add_argument("--perf", help="selected data for benchmarking/ performance analysis", action='store_true')
parser.add_argument("--extend", help="add more information for opt or perf", action='store_true')
parser.add_argument("--issue", help="print reason for terminating optimization", action='store_true')
parser.add_argument("--sort", help="sort for the key")
parser.add_argument("--revsort", help="sort reversly for the key")
parser.add_argument("--extract", help="extract only a single column. With sort, two columns are written")
parser.add_argument("--first", help="extract first iteration data instead of default last", action='store_true')
parser.add_argument('--failsafe', help="continue on errors", action='store_true')
parser.add_argument('--skipaborted', help="skip .info.xml files which are aborted", action='store_true')
parser.add_argument('--silentfailsafe', help="continue on errors w/o message", action='store_true')
args = parser.parse_args()

input = args.input if len(args.input) != 1 else glob.glob(args.input[0]) # for Windows potentially globalize 
 
res = []
if input[0].endswith('.info.xml'):
  res = read_info_xml_input(args, input)
elif input[0].endswith('.density.xml'):
  res = read_density_xml_input(args, input)  
else:
  res = read_dat_input(args, input)  

do_sort = True if args.sort or args.revsort else False 
sort_key = None if not do_sort else (args.sort if args.sort else args.revsort)

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
    for query in args.query:
      print('#query:',query)

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
