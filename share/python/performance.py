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
  def __init__(self, label, wall, cpu, sub = False):
    self.label = label
    self.wall = wall
    self.cpu = cpu
    self.sub = sub
    self.speed = cpu / wall if wall >= 3 else None
  
# helper to get extract a timer from an xml node
# extracts the data from a node. If no id attribute is given tries to be smart
def set_timer(node):  
  w = float(node.attrib['wall'])
  c = float(node.attrib['cpu'])
  # the label is more complicated. Older info.xml might not have label yet
  l = node.tag
  if l == 'timer':
    l = node.attrib['label'] if 'label' in node.attrib else None
    if l == None: # do our best
      # is this the root node?
      if node.getparent().getparent().tag == "cfsInfo":
        l = 'total'
      else:
        l = node.getparent().tag if node.getparent().tag != "summary" else node.getparent().getparent().tag
  # the optional attribute sub="true" indicates that this is sub-element 
  # and shall to be considered for missing_time
  s = node.attrib['sub'] if 'sub' in node.attrib else None
  return Timer(l, w, c, s == 'true')      
          
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
  # we have not an object with an info.xml timer status but wach info.xml is a list of Timer
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
# @timer a list of Timer or a list of a list of Timer, then the minimum is printed first 
def print_timer(timers, brief=False):
  timer = minimal_timer(timers) if type(timers[0]) == list else timers  
  meta = len(timers) if type(timers[0]) == list else 0
     
  max_label = max([len(t.label) for t in timer]) + 1 # add 1 for * in sub case
  if not brief:  
    print('TIMER (sec)'.ljust(max_label) + ': _____WALL____ ~ ______CPU______' + (' : PAR' if not meta else ''))
  else:
    print('TIMER (sec)'.ljust(max_label) + ': WALL ~    CPU')    
 
  total_wall = max(timer[0].wall, 1)
  total_cpu  = max(timer[0].cpu, 1e-3)
 
  for e in range(len(timer)):
     t = timer[e] 
     l = t.label + ('*' if t.sub else '')
     line = l.ljust(max_label) + ': {:4d}'.format(int(t.wall)) 
     if not brief:
       line += ' [{:.1%}'.format(t.wall/total_wall).rjust(8) + ']'
     line += ' ~ {:6.0f}'.format(t.cpu) if t.cpu >= 10000 else ' ~ {:6.1f}'.format(t.cpu)   
     if not brief: 
       line += '[{:.1%}'.format(t.cpu/total_cpu).rjust(8) + ']' 
     
     if not meta and timer[e].speed and timer[e].speed >= 1:
       line += ' : {:.1f}'.format(timer[e].speed)
     
     for m in range(meta):
       line += (' \t | ' if m == 0 else ' \t : ') + str(int(timers[m][e].wall)) + '\t (' + str(timers[m][e].cpu) + ')'
     print(line)
        
parser = argparse.ArgumentParser(description='when called with .info.xml the timers of this file are read. Else a performance test is run with -m and -e')
parser.add_argument("input", help="the xml file to run or the info.xml file to analyse (each with extension)")
parser.add_argument('--brief', help="brief analysis output to make it within the 1K cdash buffer", action='store_true')
parser.add_argument('-m', "--mesh", help="for execution give a mesh file for calculation, alternatively 'mesh_type' and 'res'")
parser.add_argument('--exec', help="for execution what to call for cfs", default='cfs_rel')

parser.add_argument('--repeat', help="how often shall execution be repeated - default is 1", type=int, default=1)
args = parser.parse_args()


if not os.path.exists(args.input):
  print("error: cannot open file '"  + args.input + "'")
  sys.exit(2)

if args.input.endswith('.info.xml'):
  if args.mesh or args.repeat != 1:
    print("error: when analysing an info.xml don't give mesh or repeat")
    sys.exit(1)  
  
  xml = open_xml(args.input)
  timer = read_info(xml, gap=True)
  print_timer(timer, args.brief)

else: # the whole run stuff
  if not args.mesh:
    print('error: give --mesh or --mesh_type and --res when not doing .info.xml analysis')
    sys.exit(1)            
  assert(args.input.endswith('.xml'))
  problem = args.input[:-4]
  if args.repeat == 1:
    timer = run_and_read(args.executable, args.mesh, problem + '.xml', problem)  
    print_timer(timer)
  else:
    timers = []  
    for i in range(args.repeat):
      timers.append(run_and_read(args.executable, args.mesh, problem + '.xml', problem + "-" + str(i+1)))  
    print_timer(timers)    