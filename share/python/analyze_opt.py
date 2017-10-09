#!/usr/bin/env python
import argparse
from optimization_tools import *
from csv import excel
from cfs_utils import dicts_to_gnuplot
import collections

def find_task(input):
  for file in input:
    try:
      xml = open_xml(file)
      task = xpath(xml, '//objective/@task')
      assert(task == 'minimize' or task == 'maximize')
      return task
    except:
      # do nothing
      continue
  assert(False)
    
## conditionally extract an attribute from iteration[last()] and add it to the dictionary    
def extract_float_attribute(dict, xml, attribute):
  
  query = '//iteration[last()]/@' + attribute
  if attribute in ["E_1", "E_2", "E_3", "v_21", "v_12"]:
    query = '//iteration[last()]/homogenizedTensor/orthotropy/@' + attribute
  
  if has(xml, query):
    dict[attribute] = float(xpath(xml, query))

## conditionally extract band gap data    
def extract_band_gap(dict, xml):
  for l in range(1,20):
    if has(xml, '//iteration[last()]/@bandgap_' + str(l) + '_' + str(l+1)):
      lower = float(xpath(xml, '//iteration[last()]/@ev_' + str(l) + '_max'))
      upper = float(xpath(xml, '//iteration[last()]/@ev_' + str(l+1) + '_min'))
      dict['rel_gap'] = (upper-lower)/lower # gap by lower
      dict['norm_gap'] = (upper-lower)/(.5 * (lower + upper)) # gap by center frequency
      dict['ev_' + str(l) + '_max'] = lower
      dict['ev_' + str(l+1) + '_min'] = upper
      return                                   
                                         
# process a file and inf args.info and not already recursice call the function for additional info output recursively
# @param out an array which will contain OrderedDict in the info case, otherwise it might be None
def process_file(file, args, out, recursive):
  if not os.path.exists(file):
    print('error: file not found: ' + file)
    sys.exit(3) 

  xml = None 
  try:
    xml = open_xml(file)
    # one could check for running here!
  except:
    print('caught exception working with ' + file + ': ' + str(sys.exc_info()[0]))
    if args.failsafe:
      return
    else:
      os.sys.exit(1)
  
  process = True

  q = xml.xpath('//optimization/summary/@problem')
  assert(len(q) <= 1)
  if len(q) == 0:
    if args.verbose:
      print("pass '" + file + "', no optimization summary")
    return
  msg = q[0]
  
  if args.restrict == "snopt_difficulties":
    process = True if 'difficulties' in msg else False
    if not process and args.verbose:
      print("skip '" + file + "' with summary '" + msg + "' as it has not snopt numerical difficulties")

  if process:
    # detect the order we have
    dic = collections.OrderedDict()
    base = file[:-len('.info.xml')]
    level = ord(base[-1]) - ord('a') if base[-2] == '_' and base[-1] >= 'a' and base[-1] <= 'z' else -1

    iter = int(xpath(xml, "//iteration[last()]/@number"))
    if not recursive and iter < args.min_iter:
      if args.verbose:
        print("skip '" + file + "' with " + str(iter) + " iterations")
      return

    # objective name to identify the value in 'iteration'. objective/@name added 22.2.2017 
    cost = xpath(xml, '//objective/@name') if has(xml,'//objective/@name') else xpath(xml, '//objective/@type')  
    value = float(xpath(xml, '//iteration[last()]/@' + cost))
    dic[cost] = value 
    extract_band_gap(dic, xml)
    extract_float_attribute(dic, xml, 'volume')
    extract_float_attribute(dic, xml, 'youngsModulusE1')
    extract_float_attribute(dic, xml, 'v_21')
    
    if has(xml, '//constraints/constraint[@type="curvature"][@design="profile"]/@bound_value'):
      dic['curv_prof'] = float(xpath(xml, '//constraints/constraint[@type="curvature"][@design="profile"]/@bound_value'))
    
    dic['iter'] = iter

    if not recursive and args.bound is not None:
      task = xpath(xml, '//objective/@task')
      assert(task == 'minimize' or task == 'maximize')
      if (task == 'maximize' and value < args.bound) or (task == 'minimize' and value > args.bound):
        if args.verbose:
          print("skip '" + file + "' with " + cost + "=" + str(value))
        return   

    if not recursive and not os.path.exists(base + '.density.xml'):
      print("'" + base + ".density.xml' does not exist")
      return   

    exe = xpath(xml, '//cfs/@exe')
    full = xpath(xml, '//progOpts/@meshFile') 
    mesh = full[full.rfind('/')+1:] # remove path stuff
    full = xpath(xml, '//progOpts/@parameterFile')
    prob = full[full.rfind('/')+1:] # remove path stuff
    new  = base + '_a' if level == -1 else base[:-1] + chr(ord('a') + level+1)
    
    cmd = exe + ' -m ' + mesh + ' -x ' + base + '.density.xml -p ' + prob + ' ' + new

    if args.warmstart:
      if args.warmstart == "qsub":
        print(generate_qsub_script("qsub_template.sh", cmd, new + '.sh', silent = True))
      else:
        print(cmd)
    else:
      if args.recursive:
        if level == 0:
          process_file(base[:-2] + '.info.xml', args, out, recursive=True) # base is ..._a
        if level > 0:  
          process_file(base[:-1] + chr(ord('a') + level-1) + '.info.xml', args, out, recursive=True)
          
      dic['level'] = level    
      dic['base'] = base
      if args.recursive:
        dic['issue'] = msg[msg.find('- ')+2:]       
      out.append(dic)    
      if not recursive:
        return out  
  
parser = argparse.ArgumentParser(description="The purpose of this tool is to analyze succeeded and failed optimization and optionally restart them with the current .density.xml")
parser.add_argument("input", nargs='*', help="info xml files (wildcards ok) to be ckecked for numerical aborts")
parser.add_argument('--restrict', help="optionally filter optimization results",  choices=["snopt_difficulties"])
parser.add_argument('--min_iter', help="minimal number of iterations such that we proceed", type=int, default = 0)
parser.add_argument('--bound', help="only process results not below or upper the bound, depending on the objective tast", type=float)
parser.add_argument('--verbose', help="additional information about skipped files", action='store_true')
parser.add_argument('--recursive', help="identify forerunner results", action='store_true')
parser.add_argument('--show_density', help="offer show_density.py calls", action='store_true')
parser.add_argument('--warmstart', help="generate warmstart calls for local run or RRZE HPC qsub submission. Reads info.xml for exec, mesh, ...", choices=["local", "qsub"])
parser.add_argument('--failsafe', help="continue on read error with next file", action='store_true')
parser.add_argument('--sort', help="how to sort the results (header label)", required=False)
args = parser.parse_args()

if len(args.input) == 0:
  print('usage: give .info.xml files as input')
  os.sys.exit(1)


if "qsub" == args.warmstart:
  if not os.path.exists('qsub_template.sh'):
    print("error: cannot find 'qsub_template.sh'")
    os.sys.exit(2)
  print('#!/bin/bash')
  sys.stdout.flush() # somehow necessary ?!
 
raw = [] # in the info case tuples (value, out) where out is a list of strings 
 
for file in args.input:
  res = process_file(file, args, [], recursive=False)
  if res: # skip None
    raw.append(res)
   
if not args.warmstart:

  if not args.recursive:
    # raw is a list of lists with one item each
    list = [item[0] for item in raw]    

    if args.sort:
      list = sorted(list, key=lambda k: k[args.sort])

    dicts_to_gnuplot(list)
  else:
    for entry in raw:
      dicts_to_gnuplot(entry)
      print()
      
      #if len(entry[1]) > 1:
      #  print() # empty line to separare blocks
      #for line in entry[1]:
      #  print(line)  
  