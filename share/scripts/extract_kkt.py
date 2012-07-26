#!/usr/bin/env python
# -*- coding: utf-8 -*-

import libxml2
import shutil
import sys

from cfs_utils import *

def main():

  if len(sys.argv) == 2:
    full = sys.argv[1]
    
    doc = libxml2.parseFile(full)
    xml = doc.xpathNewContext()
    
    max_iter = int(xpath(xml, "//subsolver[last()]/@iter"))
    
    print "#kkt\tdelta_obj\tsub_iters"
    for i in range(1,max_iter):
      kkt = xpath(xml, "//subsolver[@iter='" + str(i) + "']/kkt_out/@grad_max")
      old_obj = xpath(xml, "//iteration[@number='" + str(i-1) + "']/@compliance")
      new_obj = xpath(xml, "//iteration[@number='" + str(i) + "']/@compliance")
      iters = xpath(xml, "//subsolver[@iter='" + str(i) + "']/@iters")
      print kkt + "\t" + str(float(new_obj) - float(old_obj))

  else:
    print "usage: <file.info.xml>"

main()


