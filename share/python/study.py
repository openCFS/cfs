#!/usr/bin/env python

# this tool shall help to perform standard parameter studies from command linee
# to be followed by postproc.py and plotviz.py
# there are some similarities with continuation.py and warmstart.py
# e.g. study.py cfs -m bulk2d_4.mesh -p mech.xml --base base --var '//cfs:constraint[@type="volume"]/@value' --range 0.1 0.9 0.15

import sys
import argparse
import numpy as np
from cfs_utils import *


# extract the range parameters
def make_range(par):
  if len(par) < 2 or len(par) > 3:
    print('Error: --range requires 2 or 3 parameters as  <start> <end> [<step>]', par)
    sys.exit()
  
  # https://stackoverflow.com/questions/46647744/checking-to-see-if-a-string-is-an-integer-or-float
  try:
    start = float(par[0])
    end   = float(par[1])
    step  = 1.0 if len(par) < 3 else float(par[2])
    
    if start.is_integer() and end.is_integer and step.is_integer():
      r = range(int(start), int(end+1), int(step))
      return [str(t) for t in r]
    else:
      r = np.arange(start, end+step/2, step)
      return [str(round(t,8)) for t in r]
  except ValueError:
    print('Error: cannot convert to floats:', par)
    
# tries to extract the label from the query or use optinal parameter
def label(query, given):
  if given:
    return given
  
  # clean the stuff as much as possible
  for s in ['/','[', ']', '@', '"', '=', 'cfs:', 'value', 'type', 'name', 'parameter', 'param', 'constraint']:
    query = query.replace(s,'')
  
  if query == 'volume':
    return 'vol'

  if query == 'optimizer':
    return 'opt'
  
  if query == 'transferFunction':
    return 'tf'
    
  
  # extend for typical usage

  return query  
      
txt =  'perform easy parameter studies. Afterwards use postproc.py and plotviz.py.\n'      
txt += 'example: study.py cfs_master_rel -m cantilever2d_100.mesh -p simp2d.xml -b var_power_law -v //cfs:transferFunction/@param -r 1 3.01 .2\n'
txt += 'When you are happy, add --execute'
parser = argparse.ArgumentParser(description=txt)

parser.add_argument('executable', help="cfs executable")
parser.add_argument('-m', '--mesh', help="the mesh file with extension")
parser.add_argument('-x', '--initial', help="optional density.xml(.gz) for initial design (with extension)")
parser.add_argument('-p', '--problem', help="the problem -xml file used for the study", required=True)
parser.add_argument('-b', '--base', help="the common base name for the study", required=True)

parser.add_argument('-v', '--var', help='the query for the variable, e.g. \'//cfs:constraint[@type="volume"]/@value\'', required=True)
parser.add_argument('-r', '--range', nargs='+', help='the range for the variable as <start> <end> [<step>] in Python style')
parser.add_argument('-c', '--choice', nargs='+', help='range for var as alternative to --range, e.g. 1 10 100 1000')
parser.add_argument('-l', '--label', help="short label for the var variable if automatism fails")

parser.add_argument('-v2', '--var2', help='optional outer loop of var/[range/choice]/label')
parser.add_argument('-r2', '--range2', nargs='+', help='optional outer loop of var/[range/choice]/label')
parser.add_argument('-c2', '--choice2', nargs='+', help='optional outer loop of var/[range/choice]/label')
parser.add_argument('-l2', '--label2', help="optional outer loop of var/[range/choice]/label")

parser.add_argument('-v3', '--var3', help='optional most outer loop of var2/[range2/choice2]/label2')
parser.add_argument('-r3', '--range3', nargs='+', help='optional most outer loop of var2/[range2/choice2]/label2')
parser.add_argument('-c3', '--choice3', nargs='+', help='optional most outer loop of var2/[range2/choice2]/label2')
parser.add_argument('-l3', '--label3', help='optional most outer loop of var2/[range2/choice2]/label2')

parser.add_argument('--execute', help='not just printing the commands but actually executing them', action='store_true')

parser.add_argument('--redirect-output', help='redirect output to files', action='store_true')

args = parser.parse_args()

if (args.range and args.choice) or (args.range == None and args.choice == None):
  print('Error: give either --range or --choice')
  sys.exit()

if (args.var2 and not (args.range2 or args.choice2)) or (not args.var2 and (args.range2 or args.choice2)) or (args.range2 and args.choice2):
  print('Error: optional outer loop with --var2 and either --range2 or --choice2')
  sys.exit()

if (args.var3 and not (args.range3 or args.choice3)) or (not args.var3 and (args.range3 or args.choice3)) or (args.range3 and args.choice3):
  print('Error: optional most outer loop with --var3 and either --range3 or --choice3')
  sys.exit()


cmd = args.executable 
if args.mesh:
  cmd += ' -m ' + args.mesh
if args.initial:
  cmd += ' -x ' + args.initial
cmd += ' ' 

try:
  xml = open_xml(args.problem)

  l3 = label(args.var3, args.label3) if args.var3 else None
  r3 = make_range(args.range3) if args.range3 else (args.choice3 if args.choice3 else [None]) # at least one entry to enter the loop  

  l2 = label(args.var2, args.label2) if args.var2 else None
  r2 = make_range(args.range2) if args.range2 else (args.choice2 if args.choice2 else [None])   

  l1 = label(args.var, args.label)
  r1 = make_range(args.range) if args.range else args.choice  

  for v3 in r3:
    p3 = args.base
    if v3: # skip if not given
      p3 += '-' + l3 + '_' + v3
      replace(xml, args.var3, v3)  

    for v2 in r2:
      p2 = p3
      if v2:
        p2 = p3 + '-' + l2 + '_' + v2
        replace(xml, args.var2, v2)
        
      for v1 in r1: # cannot be None
        problem = p2 + '-' + l1 + '_' + v1
        replace(xml, args.var, v1, False)
        #replace(xml, '//cfs:constraint[@type="globalBucklingLoadFactor"]/@parameter', v1, False)

        xml.write(problem + '.xml')

        # setting the innermost variable as id allows easy sorting with postproc.py
        # if id is not set, we cannot simply replace(xml, ..) it
        c = cmd + problem + ' --id ' + v1
        c = cmd + '-x ' + problem[0:6] + problem[12:] + '.density.xml ' + problem + ' --id ' + v1
        
        if args.redirect_output:
          c = c + ' > ' + problem + '.txt 2>&1'
        
        if args.execute:
          execute(c, output=True)
        else:
         print(c)        

except RuntimeError as re:  
  print('Error:', re)


