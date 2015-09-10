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
parser.add_argument("msfem", help="msfem basis functions are evaluated, insert E-modulus and Poisson ratio")
parser.add_argument("shape", help="choose between frame or cross or frame with smooth corners",choices=['frame', 'cross','frame_modified','frame_w_triangles'])
parser.add_argument("--filt",help="filter for cell problem on or off",default="off")
parser.add_argument("--void",help="density of the void material",type=float,default=1e-9)
parser.add_argument("--epsilon",help="number of repetitions of base cells in msfem cell problem",type=int)
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")
parser.add_argument("--big", help="max. number of parallel cfs calculations, if turned on mtx files and vec files are not saved")
parser.add_argument("--oversampling", help="oversampling: name of the mesh with size minres/epsilon including only one base cell")



args = parser.parse_args()
stp = args.stp
dim = args.dimension
res = args.res
folder = args.folder
mesh = args.mesh
msfem = args.msfem
shape = args.shape
filt = args.filt
void = args.void

design = (' --design '+str(args.design)) if args.design else ''
oversampling =(' --oversampling '+str(args.oversampling)) if args.oversampling else ""

pwd = os.path.dirname(os.path.abspath(__file__))

execute("python calculate-crosses.py " +str(stp)+' '+str(dim)+' '+ str(res)+' '+str(folder)+ ' --msfem ' + str(mesh)+' --shape '+str(shape)+ ' --filter '+str(filt)+' --void_material '+str(void)+' --epsilon '+str(args.epsilon) + design +oversampling)

if args.big:
    execute("python evaluate.py "+str(stp)+' '+str(dim)+' '+str(res)+' '+str(folder)+ ' --msfem ' + str(msfem)+' --mesh '+str(mesh)+' --big '+str(args.big)+ ' --epsilon '+str(args.epsilon) +oversampling)
else:
  os.chdir(str(pwd)+'/'+str(folder))
  execute("bash jobs")
  os.chdir(str(pwd))
  execute("python evaluate.py "+str(stp)+' '+str(dim)+' '+str(res)+' '+str(folder)+ ' --msfem ' + str(msfem)+ ' --epsilon '+str(args.epsilon) +  design + oversampling)




