#!/usr/bin/env python
import argparse
import libxml2
import sys
import os.path
from numpy.distutils.cpuinfo import cpu

# this script has two main purposes:
# * for test-suite: extract timer data from info.xml
# * for interactive use: perform (parametric) performance tests


## this collects the information from a timer
# the data is label, wall and cpu. Extend for calls when needed
# the root Timer has the label 'total'. wall and cpu are float
class Timer:
  def __init__(self, label, wall, cpu):
    self.label = label
    self.wall = wall
    self.cpu = cpu
  
# helper to get extract a timer from an xml node
# extracts the data from a node. If no id attribute is given tries to be smart
def set_timer(node):  
  w = float(node.prop('wall'))
  c = float(node.prop('cpu'))
  # the label is more complicated. Older info.xml might not have label yet
  l = node.prop('label')
  if l == None: # do our best
    # is this the root node?
    if node.parent.parent.name == "cfsInfo":
      l = 'total'
    else:
      l = node.parent.name if node.parent.name <> "summary" else node.parent.parent.name
  return Timer(l, w,c)      
          
## extracts all timers from info.xml and give back as array of Timer objects
#@param gap shall a timer added as the last one which contains the gap from the first 
#           minus the sum of the rest
def extract_timers(xml, gap = False):
  res = []
  all = xml.xpathEval('//timer')
  for node in all:
    res.append(set_timer(node))  

  if gap:
    wall = res[0].wall  
    cpu  = res[0].cpu
    for t in res[1:]: # skip first
      wall -= t.wall # hope we stay non-negative :)
      cpu -= t.cpu
    res.append(Timer('gap', wall, cpu))

  return res

## create gnuplot output from timer
# timer list of Timer objects
def gnuplot(timer, header=True):
  if header:  
    print "#",
    for i in range(len(timer)):
      print '(' + str(i*2+1) + '):' + timer[i].label + ' \t(' + str(i*2+2) + ')<-cpu \t',
    print ""

  for t in timer:
    print str(t.wall) + ' \t' + str(t.cpu) + ' \t',
  print ""      
  
parser = argparse.ArgumentParser(description='with --analyse timer information from an info.xml is extracted.')
parser.add_argument("input", help="the xml file to run or the info.xml file to analyse")
parser.add_argument('--analyse', help="extract the timers from the 'input' info.xml file ", action='store_true')
args = parser.parse_args()


if not os.path.exists(args.input):
  print "cannot open file '"  + args.input + "'"
  sys.exit(1)

if args.analyse:
  info = libxml2.parseFile(args.input).xpathNewContext()
  timer = extract_timers(info, gap=True)
  
  print 'TIMER: WALL in s (CPU in s)'
  for t in timer:
    print t.label + ': ' + str(t.wall) + ' (' + str(t.cpu) + ')'
  