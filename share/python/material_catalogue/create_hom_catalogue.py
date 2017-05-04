#!/usr/bin/env python

import os.path
import os
from cfs_utils import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction", type=int)
parser.add_argument("dimension", help="Dimension of the problem", type=int, default=2)
parser.add_argument("res", help="resolution of the cell problem in x-direction", type=int)
parser.add_argument("folder",help="specify the folder of the h5 files")
parser.add_argument("mesh",help="name of the mesh in the folder meshes/")
parser.add_argument("shape", help="choose between frame or cross or frame with smooth corners",choices=['frame', 'cross','frame_modified','frame_w_triangles', 'auxetic','ortho_bc'])
parser.add_argument("--filt",help="filter for cell problem on or off",default="off")
parser.add_argument("--void",help="density of the void material",type=float,default=1e-9)
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")
parser.add_argument("--big", help="max. number of parallel cfs calculations, if turned on mtx files and vec files are not saved")
parser.add_argument("--penalization", help="creates a penalized material catalogue in the interval [0, 1/steps_p], step_p has to be given")
parser.add_argument("--gmsh",help="folder with stp files of geometry")
parser.add_argument("--qsub",help="option set, if job is submitted by qsub to e.g. woody")


args = parser.parse_args()
stp = args.stp
dim = args.dimension
res = args.res
folder = args.folder
mesh = args.mesh
shape = args.shape
filt = args.filt
void = args.void

pwd = os.path.dirname(os.path.abspath(__file__))
design = (' --design '+str(args.design)) if args.design else ''
penalization = (' --penalization '+str(args.penalization)) if args.penalization else ""
gmsh = (' --gmsh '+str(args.gmsh)) if args.gmsh else ''

execute("python calculate-crosses.py " +str(stp)+' '+str(dim)+' '+ str(res)+' '+str(folder)+ ' --hom ' + str(mesh)+' --shape '+str(shape)+ ' --filter '+str(filt)+' --void_material '+str(void)+ design + penalization + gmsh)



if args.big: #and args.penalization and args.design:
  execute("python evaluate.py "+str(stp)+' '+str(dim)+' '+str(res)+' '+str(folder) + ' --big '+ str(args.big) + design + penalization)
elif args.qsub:
  os.chdir(str(pwd)+'/'+str(folder))
  fobj = open("jobs")
  jobfile = open("qsub_jobs", "w")
  for line in fobj:
    problem = str(line).split()
    out = generate_qsub_script("../qsub_template.sh", line, problem[1] + '.sh', silent = True)
    jobfile.write(out+"\n")
  fobj.close()
  jobfile.close()
  os.chdir(str(pwd))
  #execute("python evaluate.py "+str(stp)+' '+str(dim)+' '+str(res)+' '+str(folder)+design+ penalization)
else:
  os.chdir(str(pwd)+'/'+str(folder))
  execute("bash jobs")
  os.chdir(str(pwd))
  execute("python evaluate.py "+str(stp)+' '+str(dim)+' '+str(res)+' '+str(folder)+design+ penalization)