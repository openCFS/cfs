#!/usr/bin/env python

# extracts a single set from a density.xml file (e.g. the last) and creates a new density file.
# practically this extracts the final result from an optimization where all iterations were written

from lxml import etree
import os
import sys
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("input", help="the input density.xml file")
parser.add_argument("output", help="the name of the density file to export")
parser.add_argument('--set', help="optional label of set, default is the last one")
args = parser.parse_args()

if not os.path.exists(args.input):
  print "input file '" + filename + "' doesn't exist"
  os.sys.exit()  
  
input = etree.parse(args.input)

header = input.getroot().xpath('//header')
if len(header) <> 1:
  print "cannot find exaxtly one 'header' in input '" + filename + "'"
  os.sys.exit()  

query = '//set[last()]' if args.set is None else '//set[@id="' + set + '"]'
set = input.getroot().xpath(query)
if len(header) <> 1:
  print "cannot find exaxtly one 'set' with name '" + args.set + "' in input '" + filename + "'"
  os.sys.exit()  

output = etree.Element("cfsErsatzMaterial")
output.append(header[0])
output.append(set[0])

etree.ElementTree(output).write(args.output)
print 'write ' + args.output