#!/usr/bin/env python
import sys
import argparse
from cfs_utils import *

## performs continuation by doubling filter/density/beta, starting from 1
# @param range_idx when we have a range with at least one doubled value. Then we add _a, _b from the second value on
# @param cnt counter of the calls (start with 0) for qsub only to define dependencies
# @param short_problem the short name of the problem where the variation will be added
# @param org_problem the name of the problem file including .xml 
def continuation(initial, cnt, type, old_var, var, mesh, short_problem, org_problem, executable, show, failsafe = False, range_idx = -1, qsub = None, plot = None, maxIter = None):

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
    if type == 'curvature':
      var_name = '-curv_'
    elif type == 'rel_profile_bound':
      var_name = '-rpb_'
    elif type == 'rel_node_bound':
      var_name = '-rnb_'      
    elif type == 'alphaSlackQuotient':
      var_name = '-asq'
    else:
      var_name = '-' + type + '_'  

    old_problem = short_problem + var_name + str(old_var) + old_ord
    var_problem = short_problem + var_name + str(var) + new_ord  
    
  # start is initial if given ofr old_var = -1 or nothing for the first run (old_var == -1)
  start = ""  
  if old_var == -1 and initial != None:
    start = "-x " + initial
  if old_var > -1:   
    start = "-x " + old_problem + dens_ext

  xml = open_xml(org_problem)     

  # we assume one hit or three for robust
  if type == 'beta':
    # check shape map first
    rsm = replace(xml, "//cfs:shapeMap/@beta", str(var), unique = False)
    if rsm == 0:
      rdf = replace(xml, "//cfs:filter/cfs:density/@beta", str(var), unique = False)
      if rdf == 0:
        raise RuntimeError("beta not found for filter/density/@beta and shapeMap/@beta")
  elif type != 'warmstart': # nothing to be done for continuation
    query = None
    if type == 'curvature':     
      query = '//cfs:constraint[@type="curvature"]/@value'
    if type == 'rel_profile_bound': 
      query = '//cfs:shapeMap/@relative_profile_bound'
    if type == 'rel_node_bound': 
      query = '//cfs:shapeMap/@relative_node_bound'
    if type == 'alphaSlackQuotient':
      query = '//cfs:constraint[@type="alphaSlackQuotient"]/@value'
    assert(query != None)

    r = replace(xml, query, str(var), unique = False)
    if r == 0:
      raise RuntimeError(" no '" + type + "' found")               
  
  if maxIter is not None:
    rmi = replace(xml, "//cfs:optimizer/@maxIterations", str(maxIter))
    if rdf == 0:
      raise RuntimeError("maxIterations not found for optimizer")
  xml.write(var_problem + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + var_problem
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
    execute(cmd, output=True, silent = failsafe)
    if plot:
      if old_var == -1:
        # add header line once
        if os.path.exists(var_problem + '.plot.dat'): 
          plot.write(first_line(var_problem + '.plot.dat', '\t var'))
      if os.path.exists(var_problem + '.plot.dat'):
        plot.write(last_line(var_problem + '.plot.dat', '\t' + str(var)))    
  if show:
    execute("show_density.py " + var_problem + dens_ext + " --save " + var_problem + ".png")

  return var_problem


# do we compress density.xml?
def compress(org_problem):
  xml = open_xml(org_problem)
  if has(xml, '//cfs:export/@compress'):
    return xpath(xml, '//cfs:export/@compress') == 'true'
  else:
    return False
      
parser = argparse.ArgumentParser(description='Make continuation, e.g. for Heaviside. In the simple case just replace cfs by continuation.py. See also warmstart.py')
parser.add_argument('input', nargs='+', help="first arguement cfs executable, second argument problem name")
parser.add_argument('-m', '--mesh', help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml(.gz) for initial design (with extension)")
parser.add_argument('-p', '--problem', help="optional problem.xml if different from second input argument")
parser.add_argument('--var', help="on which variable continuation shall be perfomed on. beta for filter or shape map", choices=['beta', 'curvature', 'rel_profile_bound', 'rel_node_bound','alphaSlackQuotient', 'warmstart'], default='beta')
parser.add_argument('--start', help="initial variable or --range", type=float, default=1)
parser.add_argument('--end', help="last variable which will be calculated or --range", type=float)
parser.add_argument('--inc', help="variable increment b += inc*b. inc=1 doubles or --range", type=float)
parser.add_argument('--step', help="variable increment b += step. Alternative to inc", type=float)
parser.add_argument('--range', help='alternative to start, end, inc/step; e.g. --range "0.01, 0.05, 0.1" or even "0.1, 0.1, 0.1"!')
parser.add_argument('--show', help="call automatically show_density.py", action='store_true')
parser.add_argument('--failsafe', help="ignore cfs exiting with error", action='store_true')
parser.add_argument('--qsub', help="template file to generate depenend job scripts for RRZE HPC (e.g. 'qsub_template.sh'")

args = parser.parse_args()

if not len(args.input) == 2:
  print('Error: usage continuation.py <executabe> <problem> -m <mesh> [--start <start> --end <end> --inc/--step <value>] / [--range <range>]')
  sys.exit(-1)

exe = args.input[0]
short_problem = args.input[1]
org_problem = short_problem + '.xml' if not args.problem else args.problem   
if not os.path.exists(org_problem):
  print("Error: file '" + org_problem + "' not found")
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
else:
  plot = open(short_problem + ".dat", "w")

if args.var == 'warmstart':
  if args.range or args.inc or args.step:
    print('Error: for --warmstart give no --range and no --inc/--step')
    sys.exit(-1)
  if not args.end:
    print('Error: for --warmstart --end is mandatory')
    sys.exit(-1)  
  vals = list(range(int(args.start if args.start else 1), int(args.end+1)))
  for idx, v in enumerate(vals):
    old_var = -1 if idx == 0 else v-1
    continuation(args.initial, cnt = idx, type = args.var, old_var=old_var, var=v, range_idx = v, mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=exe, show=args.show, failsafe=args.failsafe, qsub=args.qsub, plot = plot)
        
elif args.range:
  if args.start != 1 or args.end or args.inc or args.step:
    print("Error: don't give --start, --end, --inc or --step together with --range")
    sys.exit(-1)     
  vals = eval(args.range)   
  if len(vals) == 0 or len(vals) == 1:
    print("Errr: given --range '" + args.range + "' evaluates to " + str(len(val)) + " values")
    sys.exit(-1)
  for i in range(len(vals)):
    ri = i if len(vals) != len(set(vals)) else -1
    continuation(args.initial, cnt = i, type = args.var, old_var=-1 if i == 0 else vals[i-1], var=vals[i], mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=args.executable, show=args.show, failsafe=args.failsafe, range_idx = ri, qsub=args.qsub, plot = plot)  
else:
  if (args.inc and args.step) or (not args.inc and not args.step):
    print("Error: when no --range is give, provide either --inc or --step and --start, --end.")
    sys.exit(-1)
  if not args.end:
    print('Error: when no range is given, --end is required')
    sys.exit(-1)  
  old = -1  
  var = args.start
  i = 0
  dig = 1 if args.var == 'beta' else 2
  if args.var == 'alphaSlackQuotient':
    # decrease initial value until end value is reached, restart with smaller inc step if optimizer didn't converge
    success = args.initial
    old = -1
    dig = 6
    while var >= args.end:
      var_problem = continuation(success, type = args.var, old_var=old, var=digits(var, dig), mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=args.executable, show=args.show, failsafe=args.failsafe, plot = plot)
      infoXmlName = var_problem + ".info.xml"
      if os.path.exists(infoXmlName):
        doc_info = libxml2.parseFile(infoXmlName)
        xml_info = doc_info.xpathNewContext()
        conv  = xpath(xml_info, "//break/@converged")
      else:
        print(infoXmlName + " does not exist.")
      success = var_problem + dens_ext
      var -= args.inc
  while var <= args.end:
    maxIter = None if var == args.end else 20
    continuation(args.initial, cnt = i, type = args.var, old_var=old, var=digits(var, dig), mesh=args.mesh, short_problem=short_problem, org_problem=org_problem, executable=exe, show=args.show, failsafe=args.failsafe, qsub=args.qsub, plot=plot,maxIter=maxIter)
    old = digits(var, dig)
    if args.inc:
      var += var * args.inc
    else:
      var += args.step
    i += 1

if plot:
  print("saving meta plot file '" + short_problem + ".dat'")
