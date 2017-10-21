#!/usr/bin/env python

# The purpose of the tool is to extract region results for bloch analysis. 
# In the bloch case the results are written with step_val <wave-number>.<mode number>
# wave-number is 0-based and mode-number is 1-based  

import optimization_tools
import argparse
import os
import sys
import numpy
import copy
from lxml import etree

def get_modes(xml):
  print xml
  res = []

  for e in xml:
    l = e.get('step_val').split('.') # e.g. '28.07'
    assert(len(l) == 2)
    if float(l[0]) == 0:
      res.append(l[1])

  return res

# @param res are all result elements with their childrem item
def print_resuts(res, mode):
  total_modes = get_modes(res[0])
  modes = copy.deepcopy(total_modes)

  if mode <> None:
    if mode not in modes:
      print "error: the requested mode '" + mode + "' is not within " + str(modes)
      sys.exit()
    else:
      modes = [mode]

  header = "#(1)wave_vector"
  
  c = 2
  for e in res:
    name = e.get('data') + '_' + e.get('location')
    for m in modes:    
     header += " \t(" + str(c) + ")" + name + "_" + m
     c += 1
    
  print header
  
  
  # we need to store the data in a matrix first
  results = len(res)
  data = numpy.ones((len(res[0])/len(total_modes), results * len(modes)))
  #print data.shape
  
  # for all results
  c = 0
  for e in res:
    # for all items
    for t in e:
      wv, m = t.get('step_val').split('.')
      v = t.get('value').split(',')[0][1:] # real part of "(0.49968,6.6045e-20)", imaginary would be .split(',')[1][:-1]
      if m in modes:
        data[int(wv), len(modes) * c + modes.index(m)] = float(v)         
        # print wv, m, v
    c += 1
    
  for r in range(data.shape[0]):
    line = str(r) + " \t"
    for c in range(data.shape[1]):
      line += str(data[r,c]) + " \t"
    print line

def find_region_results(file, data, region):
  tree = etree.parse(file, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  
  dr = "" if data is None else '[@data="' + data + '"]'
  rr = "" if region is None else '[@location="' + region + '"]'
  res = root.xpath('//result[@definedOn="region"]' + dr + rr)
  if len(res) == 0:
    msg = 'no region results found'
    if data:
      msg += " with data restriction '" + data + "'"
    if region:
      msg += " with region restriction '" + region + "'"   
    print msg
    sys.exit()  
  
  return res

parser = argparse.ArgumentParser()
parser.add_argument('input', help="the info-xml file")
parser.add_argument('--data', help="optional restriction of region result name")
parser.add_argument('--region', help="optional restriction of region result region")
parser.add_argument('--mode', help="optional restriction to mode number (1-based)")
args = parser.parse_args()

if not os.path.exists(args.input):
  print 'error: file ' + args.input + ' does not exist.'
  sys.exit() 

res = find_region_results(args.input, args.data, args.region)
print_resuts(res, args.mode)