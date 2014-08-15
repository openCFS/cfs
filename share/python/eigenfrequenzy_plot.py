#!/usr/bin/env python

# this script extracts the eigenfrequencies from a eigenfrequency analysis info.xml
from lxml import etree
import os
import sys


if len(sys.argv) <> 2:
  print "usage eigenfrequency_plot.py <xxx.info.xml>"
else:
  file = sys.argv[1]
  if not os.path.exists(file):
    print "cannot open '" + file + "'"
    
  tree = etree.parse(file, etree.XMLParser(remove_comments=True))
  
  root = tree.getroot()
  result  = root.xpath("/cfsInfo/analysis/eigenFrequency/result/mode")
  # print "reading set with id = " + sett.get("id")
  length = len(result)
  
  if length == 0:
    print "no calculated eigenfrequencies found"
  else:
    print '#eigenfrequencies from ' + file 
    print '#mode\tfreq\terrror'
    for i in range(length):
      print str(result[i].get('nr')) + '\t' + str(result[i].get('frequency')) + '\t' + str(result[i].get('errorbound')) 
  
    