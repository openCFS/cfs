#!/usr/bin/env python
import argparse
import sys
import os.path
from cfs_utils import *

# this script has two main purposes:
# * for test-suite: extract timer data from info.xml
# * for interactive use: perform (parametric) performance tests


## this collects the information from a timer
# the data is label, wall and cpu. Extend for calls when needed
# the root Timer has the label 'total'. wall and cpu are float
class Timer:
  def __init__(self, label, wall, cpu, counter = -1, sub = False):
    self.label = label
    self.wall = wall
    self.cpu = cpu
    self.cnt = counter
    self.sub = sub
    self.speed = cpu / wall if wall >= 3 else None
  
# helper to get extract a timer from an xml node
# extracts the data from a node. If no id attribute is given tries to be smart
def set_timer(node):  
  w = float(node.attrib['wall'])
  c = float(node.attrib['cpu'])
  i = int(node.attrib['calls'])
  # the label is more complicated. Older info.xml might not have label yet
  l = node.tag
  if l == 'timer':
    l = node.attrib['label'] if 'label' in node.attrib else None
    if l == None: # do our best
      l = node.getparent().tag if node.getparent().tag != "summary" else node.getparent().getparent().tag
      # is this the root node?
      if l == "cfsInfo":
        l = 'total'
        i = -1 # skip counter 
  # the optional attribute sub="true" indicates that this is sub-element 
  # and shall not be considered for missing_time
  s = node.attrib['sub'] if 'sub' in node.attrib else None
  return Timer(l, w, c, i, s == 'true')      
          
## extracts all timers from info.xml and give back as array of Timer objects
#@param gap shall a timer added as the last one which contains the gap from the first 
#           minus the sum of the rest
def read_info(xml, gap = False):
  res = []
  all = xml.xpath("//*[contains(local-name(),'timer')]") # sum stuff is renamed like snopt_timer
  for node in all:
    res.append(set_timer(node))  

  if gap:
    wall = res[0].wall  
    cpu  = res[0].cpu
    for t in res[1:]: # skip first
      if not t.sub:  
        wall -= t.wall # hope we stay non-negative :)
        cpu -= t.cpu
    res.append(Timer('not_measured', wall, cpu))

  return res

## execute cfs and return the timer with gap
def run_and_read(binary, mesh, xml, problem):
  execute(binary + " -m " + mesh + " -p " + xml + " " + problem, output=True) 
  info = open_xml(problem + '.info.xml')
  timer = read_info(info, gap=True)
  return timer

## find the minimum values for each entry with a list of a list of Timer objects
# @return a timer object with minimal values. Note that gap is not recomputed
def minimal_timer(timers):
  # we have not an object with an info.xml timer status but each info.xml is a list of Timer

  # check if all timers have same length
  assert(all(len(t) == len(timers[0]) for t in timers)), 'number of timers in info.xml files differ'

  minimal = []
  meta = len(timers)  
  item = len(timers[0])
  for e in range(item):
    l = timers[0][e].label 
    w = min([t[e].wall for t in timers])     
    c = min([t[e].cpu  for t in timers])
    minimal.append(Timer(l, w, c))
  return minimal  

## create gnuplot output from timer
# timer list of Timer objects
def gnuplot(timer, header=True):
  if header:
    line = '#'  
    for i in range(len(timer)):
      line += '(' + str(i*2+1) + '):' + timer[i].label + ' \t(' + str(i*2+2) + ')<-cpu \t'
    print(line)

  line = ''
  for t in timer:
    line += str(t.wall) + ' \t' + str(t.cpu) + ' \t'
  print(line)      
  
## print standard analysis
# @param timer a list of Timer or a list of a list of Timer, then the minimum is printed first
# @param brief shorter ouput (no percentage and change)
# @param wall show wall time
# @param cpu show cpu time
# @param cnt show counter
# @param ref reference file
# @param threshold show timer only if wall >= threshold 
def print_timer(timers, brief, wall, cpu, cnt=False, ref=None, threshold=0.0):
  
  if type(timers[0]) == list and len(timers) > 1:
    timer = [minimal_timer(timers)]
    meta = len(timers)
  else:
    timer = timers
    meta = 0  

  # header
  max_label = max([len(t.label) for t in timer[0]]) + 1 # add 1 for * in sub case
  head = 'TIMER (sec)'.ljust(max_label)
  if brief:
    if not ref:
      head += ':'
      head += '   WALL' if wall else ''
      head += ' ~' if wall and cpu else '' 
      head += '   CPU ' if cpu else ''
      for m in range(meta):
        head += '  | ' if m == 0 else '  : '
        head += 'WALL_{:d}'.format(m) if wall else ''
        head += ' ~ ' if wall and cpu else ''
        head += ' CPU_{:d}'.format(m) if cpu else ''
    else:
      head += ': MIN_CPU |    REF'    
  else:
    if not ref:
      if meta:
        head += ': '
        head += '____MIN_WALL___' if wall else ''
        head += ' ~ ' if wall and cpu else ''
        head += '____MIN_CPU____' if cpu else ''
      else:
        head += ': cnt :'
        head += '______WALL_____' if wall else ''
        head += ' ~ ' if wall and cpu else ''
        head += '______CPU______' if cpu else ''
#      head += ' : PAR' if not meta else ''
      for m in range(meta):
        head += '  |' if m == 0 else '  :'
        head += ' WALL_{:d} '.format(m) if wall else ''
        head += '   ~' if wall and cpu else ''
        head += '  CPU_{:d} '.format(m) if cpu else ''
        head += '  '
    else:
      head += ':  ___MIN_CPU____  |  ___CPU_REF____  '
  print(head)

  total_wall = max(timer[0][0].wall, 1e-3)
  total_cpu  = max(timer[0][0].cpu, 1e-3)
  if ref:
    total_cpu_ref  = max(ref[0].cpu, 1e-3)
    
  # timer
  for e,t in enumerate(timer[0]):
    if t.wall < threshold and t.wall >= 0:
      continue
    # format for time display
    format_wall = '{: 6.0f}' if t.wall >= 10000 else '{: 6.2f}'
    format_cpu = '{: 6.0f}' if t.cpu >= 10000 else '{: 6.2f}'
    l = t.label + ('*' if t.sub else '')
    line = l.ljust(max_label)
    line += ':'
    if cnt:
      if t.cnt >= 0:
        line += '{:=4d} :'.format(t.cnt)
      else:
        line += '     :'.format(t.cnt)  
    if wall and not ref:
      # wall time
      line += format_wall.format(t.wall)
      if not brief:
        line += ' [{:.1%}]'.format(t.wall/total_wall).rjust(9)
    line += ' ~ ' if wall and cpu else ''
    if cpu:
      # cpu time
      line += format_cpu.format(t.cpu)   
      if not brief: 
        line += '[{:.1%}]'.format(t.cpu/total_cpu).rjust(9)
    
    if not meta and t.speed and t.speed >= 1:
      line += ' : {:.1f}'.format(t.speed)
    
    rel_eps = 0.1
    if not ref:
      # iterate over all timers
      for m in range(meta):
        line += '  | ' if m == 0 else '  : '
        if wall:
          line += format_wall.format(timers[m][e].wall)
          if not brief:
            # show change if influence on total time is > 3%
            if t.wall > total_wall * 0.03:
              rel_diff = (timers[m][e].wall - t.wall)/t.wall if t.wall != 0.0 else 0.0
              line += '(+)' if rel_diff > rel_eps else '(-)' if rel_diff < -rel_eps else '   '
            else:
              line += '   '
        line += ' ~ ' if wall and cpu else ''
        if cpu:
          line += format_cpu.format(timers[m][e].cpu)
          if not brief:
            # show change if influence on total time is > 3%
            if t.cpu > total_cpu * 0.03:
              rel_diff = (timers[m][e].cpu - t.cpu)/t.cpu if t.cpu != 0.0 else 0.0
              line += '(+)' if rel_diff > rel_eps else '(-)' if rel_diff < -rel_eps else '   '
            else:
              line += '   '
    else:
      line += '  | ' + format_cpu.format(ref[e].cpu)
      if not brief:
        line += ' [{:.1%}]'.format(ref[e].cpu/total_cpu_ref).rjust(9)
    print(line)


# searches in list of timers for timers which have not the name 'timer' and are as such ignored
def check_invalid_timers(infos):
  for info in infos:
    elem = xml.xpath("//*[@cpu]") # sum stuff is renamed like snopt_timer
    for e in elem:
      if e.tag != 'timer':
        print("Warning: ignored element '" + str(e.tag) + "': ",info)


# compare two timers
# @param eps compares relative error against eps (optional, default 10%) 
def has_rel_error(timer_ref, timer, eps=0.1, skip_noise=None):
  assert(len(timer_ref) == len(timer))

  total_cpu_ref  = max(timer_ref[0].cpu, 1e-3)
  total_cpu  = max(timer[0].cpu, 1e-3)

  error = False
  for time_ref, time in zip(timer_ref, timer):
    if time.label == 'not_measured':
      continue

    if skip_noise and time_ref.cpu < skip_noise and time.cpu < skip_noise:
      continue

    good = time_ref.cpu/total_cpu_ref
    test = time.cpu/total_cpu

    diff = (test - good)/good if good != 0.0 else test
    if diff > eps:
      print(' * error: {:s} is {:.0f}% slower than reference'.format(time.label, diff*100))
      error = True
    elif diff < -eps:
      print(' * error: {:s} is {:.0f}% faster than reference'.format(time.label, abs(diff)*100))
      error = True

  if error:
    return True

  return False
      

parser = argparse.ArgumentParser(description='when called with .info.xml the timers of this file are read. Else a performance test is run with -m and -e')
parser.add_argument("input", nargs='*', help="the xml file to run or the info.xml file(s) to analyze (each with extension)")
parser.add_argument('--brief', help="brief analysis output to make it within the 1K cdash buffer", action='store_true')
parser.add_argument('--wall', help="show only wall times", action='store_true', default=False)
parser.add_argument('--cpu', help="show only cpu times", action='store_true', default=False)
parser.add_argument('--threshold', help="show only wall time above this seconds", type=float, default=0.0)
parser.add_argument('-m', '--mesh', help="for execution give a mesh file for calculation, alternatively 'mesh_type' and 'res'")
parser.add_argument('--executable', help="for execution what to call for cfs", default='cfs_rel')
parser.add_argument('-r', '--repeat', help="how often shall execution be repeated - default is 1", type=int, default=1)
parser.add_argument("--skip_noise", help="suppress too small times", type=float, default=1e-1)
parser.add_argument('--ref', help="reference file to compare timers with")
args = parser.parse_args()

cpu = False if args.wall else True
wall = False if args.cpu else True

if args.input[0].endswith('.info.xml'):
  if args.mesh or args.repeat != 1:
    print("error: when analysing an info.xml don't give mesh or repeat")
    sys.exit(1)  
  
  timers = []
  for info_file in args.input:
    xml = open_xml(info_file)
    timers.append(read_info(xml, gap=True))

  if args.ref:
    xml = open_xml(args.ref)
    timer_ref = read_info(xml, gap=True)
    timer = minimal_timer(timers)
    if has_rel_error(timer_ref, timer, skip_noise=args.skip_noise):
      print_timer(timers, args.brief, wall, cpu, False, timer_ref)
      os.sys.exit(1)
    else:
      print('++ Times are good. ++')        
      os.sys.exit(0)
  else:
    print_timer(timers, args.brief, wall, cpu, True, None, args.threshold)
    check_invalid_timers(args.input)

else: # the whole run stuff
  assert(args.input.endswith('.xml'))
  problem = args.input[:-4]
  if not args.mesh:
    args.mesh = problem + '.mesh'
  if args.repeat == 1:
    timer = run_and_read(args.executable, args.mesh, problem + '.xml', problem)  
    if args.ref:
      xml = open_xml(args.ref)
      timer_ref = read_info(xml, gap=True)
      if has_rel_error(timer_ref, timer_tst, skip_noise=args.skip_noise):
        print_timer(timers, args.brief, wall, cpu, False, timer_ref)
    else:
      print_timer(timer, args.brief, wall, cpu)
  else:
    timers = []  
    for i in range(args.repeat):
      timers.append(run_and_read(args.executable, args.mesh, problem + '.xml', problem + "-" + str(i+1)))
    if args.ref:
      print("error: comparison of reference file with multiple runs not implemented yet")
    else:
      print_timer(timers, args.brief, wall, cpu)