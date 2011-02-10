#!/usr/bin/env python

from cfs_utils import *
import sys

def generate_headers(iters, results):
  dynamic = iters[0].hasProp('frequency') <> None
  header = "#iter" + cond(dynamic, "\tfrequency", "")
  
  for i in range(len(results)):
    res = results[i]  
    name = res.prop('data') + "_" + res.prop('location')
    
    if dynamic:
      header += "\t" + name + "_real\t" + name + "_imag"
    else:
      header += "\t" + name
      
  return header   

def extract_data(iters, item_list):
  dynamic = iters[0].hasProp('frequency') <> None
 
  out = []

  # TODO traverse the items as they might not have a 1 spacing!
  assert(len(iters) == len(item_list[0]))
   
  for i in range(len(iters)):
    iter = iters[i]
    line = iter.prop('number')
    if dynamic:
      line += "\t" + iter.prop('frequency')
        
    for r in range(len(item_list)):
      item = item_list[r][i]
      val = item.prop('value')
      if dynamic:
        line += "\t" + toGnuPlot(val)
      else:
        line += "\t" + val
    
    out.append(line)
          
  return out   


if len(sys.argv) <> 2:
  print "usage extract_info_xml_results <xxx.info.xml>"
else:
  if not os.path.exists(sys.argv[1]):
    print "cannot open '" + sys.argv[1] + "'"
  doc = libxml2.parseFile(sys.argv[1]) 
  xml = doc.xpathNewContext()   
  results = xml.xpathEval("//sequence/result")
  iters = xml.xpathEval("//process/iteration")
  
  if len(iters) > 0 and len(results) > 0:
    print generate_headers(iters, results)
    
    item_list = []
    for i in range(len(results)):
      res = results[i]
      query = "//sequence/result[@data='" + res.prop('data') + "'][@location='" + res.prop('location') + "']/item"
      items = xml.xpathEval(query)
      item_list.append(items)
    
    data = extract_data(iters, item_list)
    for i in range(len(data)):
      print data[i]
    