#!/usr/bin/env python
import argparse
from optimization_tools import *

# the purpose of this tool is to restart a bunch of problems, starting from the old solution 

def process_file(file, args, out):
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

  q = xml.xpath('//optimization/summary/@problem')
  msg = q[0] if len(q) >= 1 else ""
  
  if args.restrict == "snopt_difficulties":
    process = True if 'difficulties' in msg else False
    if not process and args.verbose:
      print("skip '" + file + "' with summary '" + msg + "' as it has not snopt numerical difficulties")

  if process:
    # detect the order we have
    base = file[:-len('.info.xml')]
    assert(base[-2] != '_' or base[-1] <= 'y') # nothing over 'z'
    level = ord(base[-1]) - ord('a') if base[-2] == '_' else -1

    iter = int(xpath(xml, "//iteration[last()]/@number"))
    if iter < args.min_iter:
      if args.verbose:
        print("skip '" + file + "' with " + str(iter) + " iterations")
      return

    # objective name to identify the value in 'iteration'. objective/@name added 22.2.2017 
    cost = xpath(xml, '//objective/@name') if has(xml,'//objective/@name') else xpath(xml, '//objective/@type')  
    value = float(xpath(xml, '//iteration[last()]/@' + cost)) 

    if args.bound is not None:
      task = xpath(xml, '//objective/@task')
      assert(task == 'minimize' or task == 'maximize')
      if (task == 'maximize' and value < args.bound) or (task == 'minimize' and value > args.bound):
        if args.verbose:
          print("skip '" + file + "' with " + cost + "=" + str(value))
        return   

    if not os.path.exists(base + '.density.xml') and not os.path.exists(base + '.density.xml.gz'):
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
      print("No density file found!")
      return
    if args.qsub:
      print(generate_qsub_script("qsub_template.sh", cmd, new + '.sh', silent = True))
    else:
      print(cmd)
      
    return (value, out)  
  
parser = argparse.ArgumentParser(description="The purpose of this tool is to restart optimization problems with the current .density.xml")
parser.add_argument("input", nargs='*', help="info xml files (wildcards ok)")
parser.add_argument('--qsub', help="RRZE job submission system", action='store_true')
parser.add_argument('--restrict', help="optionally filter optimization results",  choices=["snopt_difficulties"])
parser.add_argument('--min_iter', help="minimal number of iterations such that we proceed", type=int, default = 0)
parser.add_argument('--bound', help="only process results (objective) not below or upper the bound, depending on minimize/maximize", type=float)
parser.add_argument('--verbose', help="additional information about skipped files", action='store_true')
parser.add_argument('--failsafe', help="continue on read error with next file", action='store_true')
args = parser.parse_args()

if len(args.input) == 0:
  print('usage: give .info.xml files as input')
  os.sys.exit(1)

if args.qsub:
  if not os.path.exists('qsub_template.sh'):
    print("error: cannot find 'qsub_template.sh'")
    os.sys.exit(2)
  print('#!/bin/bash')
  sys.stdout.flush() # somehow necessary ?!
 
raw = [] # in the info case tuples (value, out) where out is a list of strings 
 
for file in args.input:
  try:
    res = process_file(file, args, [])
    if res: # skip None
      raw.append(res)
  except RuntimeError as re:
    print("Error processing ",file," -> ", str(re))
   
