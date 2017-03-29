#!/usr/bin/env python
import argparse
from optimization_tools import *
from csv import excel


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
    
# process a file and inf args.info and not already recursice call the function for additional info output recursively
# @param out an array which will contain strings in the info case, otherwise it might be None

def process_file(file, args, out, recursive):
  if not os.path.exists(file):
    raise RuntimeError('error: file not found: ' + file)

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

  msg = ""
  q = xml.xpath('//optimization/summary/@problem')
  assert(len(q) <= 1)
  if len(q) == 1:
     msg = q[0]
  
  if args.restrict == "snopt_difficulties":
    process = True if 'difficulties' in msg else False
    if not process and args.verbose:
      print("skip '" + file + "' with summary '" + msg + "' as it has not snopt numerical difficulties")

  if process:
    # detect the order we have
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

    if not recursive and args.bound is not None:
      task = xpath(xml, '//objective/@task')
      assert(task == 'minimize' or task == 'maximize')
      if (task == 'maximize' and value < args.bound) or (task == 'minimize' and value > args.bound):
        if args.verbose:
          print("skip '" + file + "' with " + cost + "=" + str(value))
        return   

    if not recursive and (not os.path.exists(base + '.density.xml') and not os.path.exists(base + '.density.xml.gz')):
      print("'" + base + ".density.xml' or .density.xml.gz' does not exist")
      return   

    exe = xpath(xml, '//cfs/@exe')
    full = xpath(xml, '//progOpts/@meshFile') 
    mesh = full[full.rfind('/')+1:] # remove path stuff
    full = xpath(xml, '//progOpts/@parameterFile')
    prob = full[full.rfind('/')+1:] # remove path stuff
    new  = base + '_a' if level == -1 else base[:-1] + chr(ord('a') + level+1)
    if os.path.exists(base + '.density.xml'):
      cmd = exe + ' -m ' + mesh + ' -x ' + base + '.density.xml -p ' + prob + ' ' + new
    elif os.path.exists(base + '.density.xml.gz'):
      cmd = exe + ' -m ' + mesh + ' -x ' + base + '.density.xml.gz -p ' + prob + ' ' + new
    else:
      print("No densityz file found!")
      return
    if args.warmstart:
      if args.warmstart == "qsub":
        print(generate_qsub_script("qsub_template.sh", cmd, new + '.sh', silent = True))
      else:
        print(cmd)
    else:
      if not args.no_recursive:
        if level == 0:
          process_file(base[:-2] + '.info.xml', args, out, recursive=True) # base is ..._a
        if level > 0:  
          process_file(base[:-1] + chr(ord('a') + level-1) + '.info.xml', args, out, recursive=True)
      line = ('iter:' if not recursive else '  >>:') + '{:3}'.format(iter) + ' ' + cost + ': {:8.4f}'.format(value) + " issue: '" + msg[msg.find('- ')+2:] 
      if args.show_density:
        line += "' -> show_density.py " + file[:-9] + ".density.xml"
      else:
        line += "' problem: '" + file + "'"
      
      out.append(line)    
      if not recursive:
        return (value, out)  
  
parser = argparse.ArgumentParser(description="The purpose of this tool is to analyze succeeded and failed optimization and optionally restart them with the current .density.xml")
parser.add_argument("input", nargs='*', help="info xml files (wildcards ok) to be ckecked for numerical aborts")
parser.add_argument('--restrict', help="optionally filter optimization results",  choices=["snopt_difficulties"])
parser.add_argument('--min_iter', help="minimal number of iterations such that we proceed", type=int, default = 0)
parser.add_argument('--bound', help="only process results not below or upper the bound, depending on the objective tast", type=float)
parser.add_argument('--verbose', help="additional information about skipped files", action='store_true')
parser.add_argument('--no_recursive', help="skip printing information of forerunner results", action='store_true')
parser.add_argument('--show_density', help="offer show_density.py calls", action='store_true')
parser.add_argument('--warmstart', help="generate warmstart calls for local run or RRZE HPC qsub submission. Reads info.xml for exec, mesh, ...", choices=["local", "qsub"])
parser.add_argument('--failsafe', help="continue on read error with next file", action='store_true')
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
  try:
    res = process_file(file, args, [], recursive=False)
    if res: # skip None
      raw.append(res)
  except RuntimeError as re:
    print("Error processing ",file," -> ", str(re))
   
if not args.warmstart:
  list = sorted(raw, key=lambda x: x[0], reverse=False if find_task(args.input) == "maxizie" else True)   
  for entry in list:
    if not args.no_recursive and len(entry[1]) > 1:
      print() # empty line to separare blocks
    for line in entry[1]:
      print(line)  
  