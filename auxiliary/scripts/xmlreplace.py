#!/usr/bin/env python
#
# xmlreplace.py
#
# This program reads in a xml-file, searches for an attribute/element (specified by
# a xpath-expression), replaces its value and saves the file.
#
# author: Andreas Hauck
# last modified: $date$

import libxml2
import sys

# Check command line and print usage
if len(sys.argv) != 4 and len(sys.argv) != 5 :
    print ""
    print "usage:"
    print "\t", sys.argv[0], " in-file xpath new-value [out-file]\n"
    print "\t", "in-file:", "\t", "name of the original xml file"
    print "\t", "xpath:", "\t\t", "xpath to the element. The prefix is 'cfs:'"
    print "\t", "new-value:", "\t", "new value of the found elements/attributes"
    print "\t", "out-file:", "\t", "optional name for the output file."
    print "\t\t\tIf no name is specified, the input file is overwritten"
    print ""
    print "This program searches in an xml input file of cfs for given element(s) / attribute(s)"
    print "and replaces their original value with a new one."
    print
    print "Examples:"
    print "\t", sys.argv[0],"cube.xml \"//cfs:transient/cfs:numSteps\"  100 new-cube.xml"
    print "\t", sys.argv[0],"cube.xml \"//cfs:mechanic//cfs:dirichletInhom/@value\"  1 new-cube.xml"
    print
    print "For further information about the xpath-syntax, please refer to"
    print "\t http://www.w3schools.com/xpath/xpath_syntax.asp"
    print
    
    sys.exit(1)

# assign command line arguments
fileName = sys.argv[1]
path = sys.argv[2]
replace = sys.argv[3]
if len(sys.argv) == 5 :
    newFileName = sys.argv[4]
else:
    newFileName = fileName
    
# setup parser interface
doc = libxml2.parseFile(fileName)
ctxt = doc.xpathNewContext()
ctxt.xpathRegisterNs('cfs', 'http://www.cfs++.org')

# look for elements
res = ctxt.xpathEval(path)
if len(res) == 0 :
    print "No elements found matching the expression"
    print "\t", path

for i in res:
  i.setContent(replace)
#  print i.content

# save file
doc.saveFile( newFileName)
