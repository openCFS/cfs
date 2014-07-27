#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This is a template for probing the most relevant ilupack parameters for
# a given problem.
# The problem itself is not given, but the adaption is easy 

import os
from cfs_utils import *

for m in ["20", "30", "40"]: # problem sizes
  for test in ["compliance", "eval_hom"]:
    for ordering in ["rcm", "amf", "amd", "mmd", "pq"]:
      for condest in ["2", "5", "10", "20", "50", "100", "200", "500"]:
        for solver in ["pcg", "sqmr", "bcgstab", "sbcg"]:
          for dropTolLU in ["1e-3", "1e-2", "1e-1"]:
            for dropTolSchur in ["1e-3", "1e-2", "1e-1"]:
              for amg in ["amli", "ilu", "mg"]:
                mesh = "periodic3d_" + m + ".mesh"
                problem = "p_" + test + "-m_" + m + "-o_" + ordering + "-c_" + condest + "-s_" + solver \
                        + "-amg_" + amg + "-dtlu_" + dropTolLU + "-dtschur_" + dropTolSchur
                xml_problem = problem + ".xml"
                print "execute " + problem + "\n"

                if os.path.exists(problem + ".info.xml"):
                  print "problem already solved"
                  continue
         
                doc = libxml2.parseFile(test + ".xml")
                xml = doc.xpathNewContext()
                xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

                replace(xml, "//cfs:ilupack/cfs:ordering", ordering)
                replace(xml, "//cfs:ilupack/cfs:condest", condest)
                replace(xml, "//cfs:ilupack/cfs:iterativeSolver/@solver", solver)
                replace(xml, "//cfs:ilupack/cfs:dropTolLU", dropTolLU)
                replace(xml, "//cfs:ilupack/cfs:dropTolSchur", dropTolSchur)
                replace(xml, "//cfs:ilupack/cfs:amg/@type", amg)
                if test == "eval_hom":
                  replace(xml, "//cfs:loadErsatzMaterial/@file", "dens_cube_" + m + ".density.xml") 
              
                doc.saveFile(xml_problem)
 
                # normally we break on an error but we want to continue if ilupack had to few iterations allowed
                try:
                  execute("cfs_master -m " + mesh + " " + problem)
                except RuntimeError as ex:
                  print "runtime exception: ", ex
