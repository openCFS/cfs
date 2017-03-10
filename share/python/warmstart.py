#!/usr/bin/env python
import argparse
from optimization_tools import *


# process a file and inf args.info and not already recursice call the function for additional info output recursively
# @param out an array which will contain strings in the info case, otherwise it might be None
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

  q = xml.xpath('//optimization/summary/@problem')
  assert(len(q) <= 1)
  if len(q) == 0:
    if args.verbose:
      print("pass '" + file + "', no optimization summary")
    return
  msg = q[0]
  
  if 'difficulties' in msg:
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

    if args.info:
      if not args.no_recursive:
        if level == 0:
          process_file(base[:-2] + '.info.xml', args, out, recursive=True) # base is ..._a
        if level > 0:  
          process_file(base[:-1] + chr(ord('a') + level-1) + '.info.xml', args, out, recursive=True)
      msg = ('iter:' if not recursive else '  >>:') + '{:3}'.format(iter) + ' ' + cost + ': {:8.4f}'.format(value) + " issue: '" + msg[msg.find('- ')+2:] + "' problem: '" + file + "'"
      out.append(msg)    
      if not recursive:
        return (value, out)  
    else:          
      if args.qsub:
        print(generate_qsub_script("qsub_template.sh", cmd, new + '.sh', silent = True))
      else:
        print(cmd)
  
parser = argparse.ArgumentParser(description="The purpose of this tool is to detect failed optimization and restart them with the current .density.xml")
parser.add_argument("input", nargs='*', help="info xml files (wildcards ok) to be ckecked for numerical aborts")
parser.add_argument('--min_iter', help="minimal number of iterations such that we proceed", type=int, default = 1)
parser.add_argument('--bound', help="only process results not below or upper the bound, depending on the objective tast", type=float)
parser.add_argument('--verbose', help="additional information about skipped files", action='store_true')
parser.add_argument('--info', help="only prints information about the status files", action='store_true')
parser.add_argument('--no_recursive', help="skip printing information of forerunner results", action='store_true')
parser.add_argument('--sort', help="sort info output by objective value", action='store_true')
parser.add_argument('--qsub', help="generate qsub scripts, requires 'qsub_template.sh'", action='store_true')
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
  res = process_file(file, args, [], recursive=False)
  if res: # skip None
    raw.append(res)
   
if args.info:
  list = sorted(raw, key=lambda x: x[0], reverse=False)   
  for entry in list:
    if not args.no_recursive:
      print() # empty line to separare blocks
    for line in entry[1]:
      print(line)  
  