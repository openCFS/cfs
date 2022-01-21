#!/usr/bin/env python
import sys
import argparse
from cfs_utils import *
import numpy as np

## performs continuation by doubling filter/density/beta, starting from 1
# @param range_idx when we have a range with at least one doubled value. Then we add _a, _b from the second value on
# @param cnt counter of the calls (start with 0) for qsub only to define dependencies
# @param type either 'warmstart' or a query
# @param short_problem the short name of the problem where the variation will be added
# @param org_problem the name of the problem file including .xml 
def continuation(initial, cnt, type, old_var, var, mesh, short_problem, org_problem, executable, show, failsafe = False, range_idx = -1, qsub = None, plot = None, maxIter = None, cmdonly = False):

  assert(range_idx < 26) # is 25 is z
  # to make use of range_idx. In the -1 case we don't use this below anyway
  old_ord = '_' + chr(ord('a') + range_idx-2) if range_idx > 1 else '' # _a for ri=2
  new_ord = '_' + chr(ord('a') + range_idx-1) if range_idx > 0 else '' # _a for ri=1
  old_problem = None
  var_problem = None

  if type == 'warmstart':
    old_problem = short_problem + old_ord
    var_problem = short_problem + new_ord  
  else:
    var_name = '-' + query_name(type) + '_'  

    old_problem = short_problem + var_name + str(old_var) + old_ord
    var_problem = short_problem + var_name + str(var) + new_ord  
    
  # start is initial if given ofr old_var = -1 or nothing for the first run (old_var == -1)
  start = ""  
  if old_var == -1 and initial != None:
    start = "-x " + initial + " "
  if old_var > -1:   
    start = "-x " + old_problem + dens_ext + " "

  xml = open_xml(org_problem)     

  if type != 'warmstart': # nothing to be done for continuation
    r = replace(xml, type, str(var), unique = False) # for robust, beta is not unique
    if r == 0:
      raise RuntimeError(" no '" + type + "' found")               
  
  if maxIter is not None:
    rmi = replace(xml, "//cfs:optimizer/@maxIterations", str(maxIter))
    if rmi == 0:
      raise RuntimeError("maxIterations not found for optimizer")
  print(var_problem)
  xml.write(var_problem + ".xml")
  
  cmd = executable + " " + start + "-m " + mesh + " " + var_problem
  if qsub:
    # see http://beige.ucs.indiana.edu/I590/node45.html
    generate_qsub_script(qsub, cmd, var_problem + '.sh', silent = True)
    if old_var == -1:
      print(('CONT' + str(cnt) + '=$(qsub ' + var_problem + '.sh)'))
      print(('echo $CONT' + str(cnt)))
    else:  
      print(('CONT' + str(cnt) + '=$(qsub -W depend=afterok:$CONT' + str(cnt-1) + ' ' + var_problem + '.sh)'))
      print(('echo $CONT' + str(cnt)))
  else:
    if cmdonly:
      print(cmd)
    else:     
      execute(cmd, output=True, silent = failsafe)
      if plot:
        if old_var == -1:
          # add header line once
          if os.path.exists(var_problem + '.plot.dat'): 
            plot.write(first_line(var_problem + '.plot.dat', '\t var'))
        if os.path.exists(var_problem + '.plot.dat'):
          plot.write(last_line(var_problem + '.plot.dat', '\t' + str(var)))    
  if show and not cmdonly:
    execute("show_density.py " + var_problem + dens_ext + " --save " + var_problem + ".png")

  return var_problem

# helper which gives the query from prefined var or original query if not given
def query(var, query):
  assert not (var is None and query is None)
  assert not (var is not None and query is not None)
  if query:
    return query
  if var == 'penalty':
    return '//cfs:transferFunction/@param'
  if var == 'beta':
    return '//cfs:filter/cfs:density/@beta'
  print('Error: unknown predefined var ' + var)
  sys.exit(-1)

# try to extract a short name from the query.
# extend for your needs
def query_name(query):
  if 'transferFunction/@param' in query:
    return 'pen'
  if '@' in query: 
    query = query[query.index('@')+1:] # includes --var is 'beta'#
  if '=' in query: # for //cfs:constraint[@type='volume']/@value use the volume
    query[query.index('=')+1:query.index(']')]    
    # clean the stuff as much as possible
  for s in ["'", '/','[', ']', '@', '"', '=', 'item', 'result', 'data', 'value', 'type', 'name', 'param', 'constraint']:
    query = query.replace(s,'')     
  return query


# do we compress density.xml?
def compress(org_problem):
  xml = open_xml(org_problem)
  if has(xml, '//cfs:export/@compress'):
    return xpath(xml, '//cfs:export/@compress') == 'true'
  else:
    return False
      
      
txt =  'Make continuation, e.g. for Heaviside. In the simple case just replace cfs by continuation.py. See also warmstart.py.\n'
txt += 'example: continuation.py cfs_rel killme -p mech.xml -m bulk2d_10.mesh --var penalty --range 1 4  --cmdonly'
parser = argparse.ArgumentParser(description=txt)
parser.add_argument('exe', help="cfs executable")
parser.add_argument('name', nargs='?', help="optional problem name")
parser.add_argument('-m', '--mesh', help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml(.gz) for initial design (with extension)")
parser.add_argument('-p', '--problem', help="optional problem.xml")
parser.add_argument('--var', help="predefined query for continuation variable. Together with one of range/factor/list", choices=['beta', 'penalty'])
parser.add_argument('--query', help='alternative to --var an open query like //cfs:constraint[@type="volume"]/@value')
parser.add_argument('--range', nargs='+', help='start stop [inc] for step + inc', type=float)
parser.add_argument('--factor', nargs=3, help='start stop factor for step * factor',type=float)
parser.add_argument('--list', nargs='+', help='val1 val2 val3 ...',type=float)
parser.add_argument('--warmstart', nargs='*', help='instead of var and rang/factor/list run the problem from [start] to end', type=int)
parser.add_argument('--show', help="call automatically show_density.py", action='store_true')
parser.add_argument('--failsafe', help="ignore cfs exiting with error", action='store_true')
parser.add_argument('--noplot', help="don't create a meta plot file", action='store_true')
parser.add_argument('--cmdonly', help="command output, no execution", action='store_true')
parser.add_argument('--qsub', help="template file to generate dependent job scripts for RRZE HPC (e.g. 'qsub_template.sh'")

args = parser.parse_args()

if not args.name and not args.problem:
  print('Error: give name and/or problem (-p/--problem)')
  sys.exit(-1)

short_problem = args.name if args.name else args.problem[:-len('.xml')]
org_problem = short_problem + '.xml' if not args.problem else args.problem   
if not os.path.exists(org_problem):
  print("Error: problem file '" + org_problem + "' not found")
  sys.exit(-1)
 
# the density file extension
dens_ext = '.density.xml'
if compress(org_problem):
  dens_ext += '.gz'
  if args.initial and not args.initial.endswith('.gz'):
    # but it could be implemented
    print('Error: created .density.xml are compressed, but not the initial one')
    sys.exit(-1) 

# the plot dat file where we collect the last lineas. Not for qsub
plot = None

if args.qsub:
  if not os.path.exists(args.qsub):
    print('Error: qsub template file not found ' + args.qsub)
    sys.exit(-1)
  print('#!/bin/bash')  
elif not args.noplot and not args.cmdonly:
  plot = open(short_problem + ".dat", "w")

if args.warmstart and (args.var or args.query or args.range or args.list or args.factor):
  print('Error: do not give warmstart with var/query and range/factor/list')
  sys.exit(-1)

if args.warmstart:
  if len(args.warmstart) > 2:
    print('Error: --warmstart [start] end has one or two arguments')
    sys.exit(-1)  
  start = 1 if len(args.warmstart) == 1 else args.warmstart[0]
  end   = args.warmstart[0] if len(args.warmstart) == 1 else args.warmstart[1]

  vals = list(range(start, end+1))
  for idx, v in enumerate(vals):
    old_var = -1 if idx == 0 else v-1
    continuation(args.initial, cnt = idx, type = 'warmstart', old_var=old_var, var=v, range_idx = v, mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=args.exe, show=args.show, failsafe=args.failsafe, qsub=args.qsub, plot = plot, cmdonly = args.cmdonly)
else:
  if args.var and args.query:
    print("Error: don't give --var and --query the same time")
    sys.exit(-1)
  # no warmtart
  if not args.var and not args.query:
    print('Error: either --warmstart [S] N or --var/query *')     
    sys.exit(-1)
  if (args.range and (args.factor or args.list)) or (args.list and(args.factor or args.range)) or (not args.range and not args.factor and not args.list):
    print('Error: with var give exatly one of range, factor, list')
    sys.exit(-1)       

  vals = []
  if args.list:
    vals = args.list
  else:
   if args.range:
     if len(args.range) < 2 and len(args.range) > 3:
       print('Error: range has parameters <start> <end> [inc]')
     # add eps to have end included  
     vals = list(np.arange(args.range[0], args.range[1] + 1e-10, args.range[2] if len(args.range) > 2 else 1)) 
   else:  
      v, e, f = args.factor
      while (v <= e) if f > 1 else (v >= e):
        vals.append(v)
        v *= f
   # 1.9000000000000001 to 1.9       
   vals = [round(v,10) for v in vals]
   # test if we have integers only
   if np.sum([abs(math.modf(v)[0]) for v in vals]) == 0:
     vals = [int(v) for v in vals] 
   if not vals or len(vals) == 0 or len(vals) == 1:
     print("Error: value given by range/factor/list evaluations only to",vals)
     sys.exit(-1)
  #print(vals)  
  for i in range(len(vals)):
    ri = i if len(vals) != len(set(vals)) else -1
    continuation(args.initial, cnt = i, type = query(args.var, args.query), old_var=-1 if i == 0 else vals[i-1], var=vals[i], mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=args.exe, show=args.show, failsafe=args.failsafe, range_idx = ri, qsub=args.qsub, plot = plot, cmdonly = args.cmdonly)  
if plot:
  print("saving meta plot file '" + short_problem + ".dat'")
