#!/usr/bin/env python
import sys
import libxml2
import argparse

from cfs_utils import *

## performs heaviside continuation by doubling filter/density/beta, starting from 1
# @param range_idx when we have a range with at least one doubled value. Then we add _a, _b from the second value on 
def continuation(initial, type, old_var, var, mesh, short_problem, executable, show, failsafe = False, range_idx = -1):

  assert(range_idx < 26) # is 25 is z
  # to make use of range_idx. In the -1 case we don't use this anysway below
  old_var_str = str(old_var) + ('_' + chr(ord('a') + range_idx-2) if range_idx > 1 else '') # _a for ri=2 
  new_var_str = str(var) + ('_' + chr(ord('a') + range_idx-1) if range_idx > 0 else '') # _a for ri=1
  print range_idx, old_var_str, new_var_str
  var_name = '-' + type + '_'
  if type == 'curvature':
    var_name = '-curv_'
  if type == 'rel_profile_bound':
    var_name = '-rpb_'
  if type == 'rel_node_bound':
    var_name = '-rnb-'      
  var_problem = short_problem + var_name + new_var_str  
    
  # first without iniital  
  start = "" if old_var < 1 else "-x " + short_problem + var_name + old_var_str + ".density.xml"
  # now check for initial
  if old_var < 1 and initial <> None:
    assert(start == "")
    start = "-x " + initial    

  if not os.path.exists(short_problem + ".xml"):
    print "error: file '" + short_problem + ".xml' not found"
    os.sys.exit()   
  doc = libxml2.parseFile(short_problem + ".xml")
  xml = doc.xpathNewContext()
  xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

  # we assume one hit or three for robust
  if type == 'beta':
    res = xml.xpathEval("//cfs:filter/cfs:density/@beta")
    if len(res) == 0:
      res = xml.xpathEval("//cfs:shapeMap/@beta")
      if len(res) == 0:  
        raise RuntimeError("beta not found for filter/density/@beta and shapeMap/@beta")
    for data in res:
      data.setContent(str(var))
  else:
    query = None
    if type == 'curvature':     
      query = '//cfs:constraint[@type="curvature"]/@value'
    if type == 'rel_profile_bound': 
      query = '//cfs:shapeMap/@relative_profile_bound'
    if type == 'rel_node_bound': 
      query = '//cfs:shapeMap/@relative_node_bound'
    assert(query <> None)
               
    res = xml.xpathEval(query)
    if len(res) == 0:
      raise RuntimeError(" no '" + type + "' found")
    for data in res:
      by_nx = str(data).find('/nx') > 0 and str(var).find('/nx') == -1
      data.setContent(str(var) + ('/nx' if by_nx else ''))
  
  doc.saveFile(var_problem + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + var_problem
  execute(cmd, output=True, silent = failsafe)
      
  if show:
    execute("show_density.py " + var_problem + ".density.xml --save " + var_problem + ".png")

  return


parser = argparse.ArgumentParser(description='Make continuation, e.g. for Heaviside. In the simple case just replace cfs by continuation.py')
parser.add_argument('problem', help="the problem xml without extension where e.g. '-beta_x' will be added")
parser.add_argument('-m', "--mesh", help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml for initial design (with extension)")
parser.add_argument('--var', help="on which variable continuation shall be perfomed on. beta for filter or shape map", choices=['beta', 'curvature', 'rel_profile_bound', 'rel_node_bound'], default='beta')
parser.add_argument('--start', help="initial variable", type=float, default=1)
parser.add_argument('--max', help="maxmum variable which will be calculated", type=float, default=64)
parser.add_argument('--inc', help="variable increment b += inc*b. inc=1 doubles", type=float, default=1.0)
parser.add_argument('--range', help='alternative to start, max, inc i like --range "0.01, 0.05, 0.1" or even "0.1, 0.1, 0.1"!')
parser.add_argument('--executable', help="what to call for cfs", default='cfs_rel')
parser.add_argument('--noshow', help="suppress calling show_density.py, e.g. for 3d!", action='store_true')
parser.add_argument('--failsafe', help="ignore cfs exiting with error", action='store_true')

args = parser.parse_args()

if args.range:
  vals = eval(args.range)   
  if len(vals) == 0 or len(vals) == 1:
    print "given --range '" + args.range + "' evaluates to " + str(len(val)) + " values"
    sys.exit(-1)  
  for i in range(len(vals)):
    ri = i if len(vals) <> len(set(vals)) else -1  
    continuation(args.initial, type = args.var, old_var=-1 if i == 0 else vals[i-1], var=vals[i], mesh=args.mesh, short_problem=args.problem, executable=args.executable, show=not args.noshow, failsafe=args.failsafe, range_idx = ri)  
  
else:
  old = -1  
  var = args.start

  dig = 1 if args.var == 'beta' else 2
  while var <= args.max:
   continuation(args.initial, type = args.var, old_var=old, var=digits(var, dig), mesh=args.mesh, short_problem=args.problem, executable=args.executable, show=not args.noshow, failsafe=args.failsafe)
   old = digits(var, dig)
   var += var * args.inc
