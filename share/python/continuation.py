#!/usr/bin/env python
import sys
import argparse
from cfs_utils import *

## performs continuation by doubling filter/density/beta, starting from 1
# @param range_idx when we have a range with at least one doubled value. Then we add _a, _b from the second value on 
def continuation(initial, cnt, type, old_var, var, mesh, short_problem, executable, show, failsafe = False, range_idx = -1, qsub = None):

  assert(range_idx < 26) # is 25 is z
  # to make use of range_idx. In the -1 case we don't use this below anyway
  old_var_str = str(old_var) + ('_' + chr(ord('a') + range_idx-2) if range_idx > 1 else '') # _a for ri=2 
  new_var_str = str(var) + ('_' + chr(ord('a') + range_idx-1) if range_idx > 0 else '') # _a for ri=1
  #print range_idx, old_var_str, new_var_str
  var_name = '-' + type + '_'
  if type == 'curvature':
    var_name = '-curv_'
  if type == 'rel_profile_bound':
    var_name = '-rpb_'
  if type == 'rel_node_bound':
    var_name = '-rnb_'      
  old_problem = short_problem + var_name + old_var_str  
  if type == 'alphaSlackQuotient':
    var_name = '-asq'

  var_problem = short_problem + var_name + new_var_str  
    
  # start is initial if given ofr old_var = -1 or nothing for the first run (old_var == -1)
  start = ""  
  if old_var == -1 and initial <> None:
    start = "-x " + initial
  if old_var > -1:   
    start = "-x " + old_problem + ".density.xml"

  if not os.path.exists(short_problem + ".xml"):
    print "error: file '" + short_problem + ".xml' not found"
    os.sys.exit()
  
  xml = open_xml(short_problem + ".xml")     

  # we assume one hit or three for robust
  if type == 'beta':
    # check shape map first
    rsm = replace(xml, "//cfs:shapeMap/@beta", str(var), unique = False)
    if rsm == 0:
      rdf = replace(xml, "//cfs:filter/cfs:density/@beta", str(var), unique = False)
      if rdf == 0:
        raise RuntimeError("beta not found for filter/density/@beta and shapeMap/@beta")
  else:
    query = None
    if type == 'curvature':     
      query = '//cfs:constraint[@type="curvature"]/@value'
    if type == 'rel_profile_bound': 
      query = '//cfs:shapeMap/@relative_profile_bound'
    if type == 'rel_node_bound': 
      query = '//cfs:shapeMap/@relative_node_bound'
    if type == 'alphaSlackQuotient':
      query = '//cfs:constraint[@type="alphaSlackQuotient"]/@value'
    assert(query <> None)

    r = replace(xml, query, str(var), unique = False)
    if r == 0:
      raise RuntimeError(" no '" + type + "' found")               
  
  xml.write(var_problem + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + var_problem
  if qsub:
    # see http://beige.ucs.indiana.edu/I590/node45.html
    generate_qsub_script(qsub, cmd, var_problem + '.sh', silent = True)
    if old_var == -1:
      print('CONT' + str(cnt) + '=$(qsub ' + var_problem + '.sh)')
      print('echo $CONT' + str(cnt))
    else:  
      print('CONT' + str(cnt) + '=$(qsub -W depend=afterok:$CONT' + str(cnt-1) + ' ' + var_problem + '.sh)')
      print('echo $CONT' + str(cnt))
  else:
    execute(cmd, output=True, silent = failsafe)
      
  if show:
    execute("show_density.py " + var_problem + ".density.xml --save " + var_problem + ".png")

  return var_problem


parser = argparse.ArgumentParser(description='Make continuation, e.g. for Heaviside. In the simple case just replace cfs by continuation.py')
parser.add_argument('problem', help="the problem xml without extension where e.g. '-beta_x' will be added")
parser.add_argument('-m', "--mesh", help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml for initial design (with extension)")
parser.add_argument('--var', help="on which variable continuation shall be perfomed on. beta for filter or shape map", choices=['beta', 'curvature', 'rel_profile_bound', 'rel_node_bound','alphaSlackQuotient'], default='beta')
parser.add_argument('--start', help="initial variable", type=float, default=1)
parser.add_argument('--end', help="last variable which will be calculated", type=float, default=64)
parser.add_argument('--inc', help="variable increment b += inc*b. inc=1 doubles", type=float, default=1.0)
parser.add_argument('--range', help='alternative to start, max, inc i like --range "0.01, 0.05, 0.1" or even "0.1, 0.1, 0.1"!')
parser.add_argument('--executable', help="what to call for cfs", default='cfs_rel')
parser.add_argument('--noshow', help="suppress calling show_density.py, e.g. for 3d! (standard for qsub)", action='store_true')
parser.add_argument('--failsafe', help="ignore cfs exiting with error", action='store_true')
parser.add_argument('--qsub', help="template file to generate depenend job scripts for RRZE HPC (e.g. 'qsub_template.sh'")

args = parser.parse_args()

if args.qsub:
  if not os.path.exists(args.qsub):
    print('qsub template file not found ' + args.qsub)
    os.sys.exit(1)
  args.noshow = True
  
if args.range:
  vals = eval(args.range)   
  if len(vals) == 0 or len(vals) == 1:
    print "given --range '" + args.range + "' evaluates to " + str(len(val)) + " values"
    sys.exit(-1)  
  for i in range(len(vals)):
    ri = i if len(vals) <> len(set(vals)) else -1  
    continuation(args.initial, cnt = i, type = args.var, old_var=-1 if i == 0 else vals[i-1], var=vals[i], mesh=args.mesh, short_problem=args.problem, executable=args.executable, show=not args.noshow, failsafe=args.failsafe, range_idx = ri, qsub=args.qsub)  
else:
  old = -1  
  var = args.start
  i = 0
  dig = 1 if args.var == 'beta' else 2
  if args.var == 'alphaSlackQuotient':
    # decrease initial value until end value is reached, restart with smaller inc step if optimizer didn't converge
    inc = args.inc
    success = args.initial
    old = -1
    dig = 6
    while var >= args.end:
      var_problem = continuation(success, type = args.var, old_var=old, var=digits(var, dig), mesh=args.mesh, short_problem=args.problem, executable=args.executable, show=not args.noshow, failsafe=args.failsafe)
      infoXmlName = var_problem + ".info.xml"
      if os.path.exists(infoXmlName):
        doc_info = libxml2.parseFile(infoXmlName)
        xml_info = doc_info.xpathNewContext()
        conv  = xpath(xml_info, "//break/@converged")
      else:
        print infoXmlName + " does not exist."
      #old = digits(var, dig)
      #if conv == 'yes':
      success = var_problem + ".density.xml"
      #  inc *= 1.5
      var -= inc
      #else:
      #  inc /= 2.
      #  var += inc     
  while var <= args.end:
   continuation(args.initial, cmt = i, type = args.var, old_var=old, var=digits(var, dig), mesh=args.mesh, short_problem=args.problem, executable=args.executable, show=not args.noshow, failsafe=args.failsafe, qsub=args.qsub)
   old = digits(var, dig)
   var += var * args.inc
   i += 1
