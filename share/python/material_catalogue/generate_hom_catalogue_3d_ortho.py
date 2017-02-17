#!/usr/bin/env python

# this script generates a qsub shell script for submitting homogenization jobs on (woody) cluster
# one homogenizaton job consists of generating a basecell with 'basecell.py' and the material tensor calculation with cfs 
import argparse
import numpy as np
import os.path
import cfs_utils

parser = argparse.ArgumentParser()
parser.add_argument("--steps", help="number of grid points in one direction", type=int,required=True)
parser.add_argument("--res", help="resolution of the base cell in x-direction", type=int,required=True)
parser.add_argument("--xml", help="specify xml name for homogenization, should end with *.xml",required=True)

args = parser.parse_args()
steps = args.steps
xml = args.xml

# sanity checks
if not args.xml.endswith(".xml"):
  print("xml file must end with .xml!")
  sys.exit(0)
  
assert(os.path.exists(xml))
assert(os.path.exists("qsub_template.sh"))  
assert(os.path.exists("mat.xml"))
assert(steps > 0)

jobfile = open("qsub_jobs.sh", "w")
jobfile.write("#!/bin/bash\n")

folder = "shell_scripts"
if os.path.exists(folder):
  os.system("rm -r -f " + folder)
os.mkdir(folder)

for i,x1 in enumerate(np.arange(0,1.1,1.0/float(steps))):
  for j,y1 in enumerate(np.arange(0,1.1,1.0/float(steps))):
    for k,z1 in enumerate(np.arange(0,1.1,1.0/float(steps))):
      # basecell.py only take values 0 < x < 1, thus project values at interval boundary
      x = x1
      if i == 0:
        x = 1e-3
      elif i == steps + 1:
        x = 0.99
      
      y = y1
      if j == 0:
        y = 1e-3
      elif j == steps + 1:
        y = 0.99
        
      z = z1
      if k == 0:
        z = 1e-3
      elif k == steps + 1:
        z = 0.99
      
      problem = str(i) + "-" + str(j) + "-" + str(k)
      mesh = problem + ".mesh"
      
      cmd = "basecell.py --res " + str(args.res) + " --x1 " + str(x) + " --y1 " + str(y) + " --z1 " + str(z)+" --target volume_mesh --beta 7 --eta 0.6  --interpolation heaviside --save " + problem
      cmd += "; cfs_rel -m " + mesh + " -p " + xml + " " + problem
      cmd += "&& rm " + mesh  
      
      out = cfs_utils.generate_qsub_script("qsub_template.sh", cmd, folder+"/"+problem+'.sh', silent = True)
      jobfile.write(out+"\n")
      
jobfile.close()      