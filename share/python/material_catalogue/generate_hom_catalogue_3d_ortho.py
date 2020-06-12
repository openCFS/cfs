#!/usr/bin/env python

# this script generates a qsub shell script for submitting homogenization jobs on (woody) cluster
# one homogenizaton job consists of generating a basecell with 'basecell.py' and the material tensor calculation with cfs 
import argparse
import numpy as np
import os.path
import cfs_utils

parser = argparse.ArgumentParser()
parser.add_argument("--steps", help="number of grid points in one direction", type=int, required=True)
parser.add_argument("--res", help="resolution of the base cell in x-direction", required=True)
parser.add_argument("--xml", help="specify xml name for homogenization, should end with *.xml", required=True)
parser.add_argument("--min_diam", help="min diameter of basecell (default=1e-3)", type=float, default=1e-3)
parser.add_argument("--max_diam", help="max diameter of basecell (default=0.9)", type=float, default=0.91)
parser.add_argument("--bend", help="bending of spline, 0 is linear, 1 is corner (default=0.5)", type=float, default=0.5)

args = parser.parse_args()
steps = args.steps
xml = args.xml
 
eps = 1e-6
assert(args.min_diam > -eps and args.min_diam < 1)
assert(args.max_diam > eps and args.max_diam < 1 + eps)
 
# sanity checks
if not args.xml.endswith(".xml"):
  print("xml file must end with .xml!")
  sys.exit(0)
 
assert(os.path.exists(xml))
assert(os.path.exists("qsub_template.sh"))  
assert(os.path.exists("mat.xml"))
assert(steps > 0)

cwd = os.getcwd()
dens_folder = cwd + "/density-xml"
if os.path.exists(dens_folder):
  os.system("rm -r -f " + dens_folder)
os.mkdir(dens_folder)
 
jobfile = open("submit_jobs.sh", "w")
jobfile.write("#!/bin/bash\n")
 
script_folder = cwd + "/shell_scripts"
if os.path.exists(script_folder):
  os.system("rm -r -f " + script_folder)
os.mkdir(script_folder)
 
cwd = os.getcwd()
 
# range_diam = args.max_diam - args.min_diam
range_diam = 1
assert(range_diam > 0)
print("range_diam", range_diam)

out = os.system("create_mesh.py --type bulk3d --res " + args.res)
assert(out == 0)
bulk_mesh = "bulk3d_" + args.res + ".mesh"
# for i,x1 in enumerate(np.arange(args.min_diam ,args.max_diam+eps,range_diam/float(steps))):
#   for j,y1 in enumerate(np.arange(args.min_diam ,args.max_diam+eps,range_diam/float(steps))):
#     for k,z1 in enumerate(np.arange(args.min_diam ,args.max_diam+eps,range_diam/float(steps))):
for i, x1 in enumerate(np.arange(0, 1.05, 1.0 / float(steps))):
  for j, y1 in enumerate(np.arange(0, 1.05, 1.0 / float(steps))):
    for k, z1 in enumerate(np.arange(0, 1.05, 1.0 / float(steps))):
      
      if i == 0 and j == 0 and k == 0 or i == steps and j == steps and k == steps:
      # if i != steps-1 and j != steps-1 and k != steps-1 or i == 0 or j == 0 or k == 0:
        continue
      
      problem = str(i) + "-" + str(j) + "-" + str(k)
      
      p2 = str(i) + "-" + str(k) + "-" + str(j)
      p3 = str(j) + "-" + str(i) + "-" + str(k)
      p3 = str(j) + "-" + str(k) + "-" + str(i)
      p4 = str(k) + "-" + str(i) + "-" + str(j)
      p5 = str(k) + "-" + str(j) + "-" + str(i)
      
      skip_problem = False
      # skip symmetrical problems
      for p in [problem, p2, p3, p4, p5]:
        problemname = cwd + "/" + script_folder + "/" + p + ".sh"
        if os.path.exists(problemname):
          print("skip ", p)
          skip_problem = True
          break
      
      if skip_problem:
        continue
      
      problemname = script_folder + "/" + problem
      
      cmd = None
      # don't use basecell script but use create_density to generate density files for one or two intersecting cylinders
      if i == 0 or j == 0 or k == 0 or i == steps or j == steps or k == steps:
        rad1 = x1 / 2.0
        rad2 = y1 / 2.0
        rad3 = z1 / 2.0
        str_rad1 = str(round(rad1, 3))
        str_rad2 = str(round(rad2, 3))
        str_rad3 = str(round(rad3, 3))
        cmd = "create_density.py --res " + args.res + " --dim 3 --lower 1e-9 --three_cylinders --radius " + str_rad1 + "," + str_rad2 + "," + str_rad3 + " --save " + dens_folder + "/" + problem + ".density.xml"
        cmd += " && cfs_rel -d -m " + bulk_mesh + " -x " + dens_folder + "/" + problem + ".density.xml -p homogenize.xml " + problem
        cmd += " && rm " + dens_folder + "/" + problem + ".density.xml"
      else:
        mesh = problem + ".mesh"
        cmd = "basecell.py --res " + str(args.res) + " --x1 " + str(x1) + " --y1 " + str(y1) + " --z1 " + str(z1) + " --target volume_mesh --beta 7 --eta 0.6  --interpolation heaviside --save " + problem
        cmd += " --bend " + str(args.bend) + "  --stiffness_as_diameter &&"
        cmd += " cfs_rel -d -m " + mesh + " -p " + cwd + "/" + xml + " " + problem
        cmd += "&& rm " + mesh + " && rm " + problem + ".vtp"
      
      out = cfs_utils.generate_qsub_script("qsub_template.sh", cmd, problemname + ".sh", silent=True)
      jobfile.write(out + "\n")
      
jobfile.close()
os.system("chmod 755 submit_jobs.sh")      
