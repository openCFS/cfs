#!/usr/bin/env python

import sys, os, string
import numpy as np
import scipy.io
from numpy import dot, empty, size
from numpy import sin
from numpy import cos
from numpy import sqrt
from lxml import etree
import argparse
import scipy, scipy.sparse, scipy.sparse.linalg 
from hdf5_tools import *
from optimization_tools import *
from mesh_tool import *
import matviz_rot
from decimal import *
import os.path
import os
from cfs_utils import *
from subprocess import PIPE
import subprocess
import time


def submit_job_max_len(job_list, max_processes):
  # run multiple processes at the same time. number limited by max_processes
  sleep_time = 0.1
  FNULL = open(os.devnull, 'w')
  processes = list()
  for command in job_list:
    print('running {n} processes. Submitting {proc}'.format(n=len(processes),
        proc=str(command)))
    processes.append(subprocess.Popen(command, shell=False, stdout=FNULL,
      stdin=PIPE,close_fds=True))
    while len(processes) >= max_processes:
      time.sleep(sleep_time)
      processes = [proc for proc in processes if proc.poll() is None]
  while len(processes) > 0:
    time.sleep(sleep_time)
    processes = [proc for proc in processes if proc.poll() is None]

def CalcElementMat(R,Kmat,dof,dim,l=None,rep=None):
  # calculate MSFEM element matrix by applying restriction matrices overR to the global FEM
  # stiffness matrix of the cell problem
  
  # uniform quad element
  dx = 1.
  dy = 1.
  
  #Kglob is not necessary anymore
  #K = Kglob(E, nu, dx, dy, dim, dof, mesh)
  
  # stiffness matrix is read from CFS cell problems
  K = scipy.io.mmread(Kmat)
  # scipy.io.savemat('K.mat', mdict={'K': K})
  #scipy.io.savemat('tmp.mat', mdict={'K': tmp})
  
  eps = 1e-6
  if rep == None:
    # calculate MSFEM macro element matrix without oversampling
    elemMat = R * K * R.T    
  else:
    # calculate MSFEM macro element matrix with oversampling
    overR = np.matrix(np.zeros((dof, dim * (int(l/rep + 0.5 + eps)+1)**2)))
    # extract deformations from interior base cell from oversampled results
    for i in range(dof):
      count = int(l/rep *(l+1) + l/rep + 0.5 + eps)
      count2 = 0
      for j in range(int(l/rep + 0.5 + eps)+1):
        overR[i,dim*count2:dim*(count2 + int(l/rep + 0.5 + eps)+1)] = R[i,dim*count:dim*(count+int(l/rep + 0.5 + eps)+1)]
        count += (l + 1)
        count2 += int(l/rep + 0.5 + eps)+1
    elemMat = overR * K * overR.T
    #print "R = "+str(R)
    #print "overR = "+str(overR)
  return elemMat

def get_R(h5file,i,R,dof,dim,l,over = False,rep=None,C=None):
  # insert solution of the MSFEM 2D cell problem in the restriction matrix R 
  # h5file: h5 file with the cell problem solution for problem number i
  # dof: degrees of freedom of the cell problem
  # l: number of elements in x-direction for the cell problem
  # rep: number of base cell repetitions in x-direction for MSFEM oversampling 
  # matrix C necessary for MSFEM oversampling
  
  f = h5py.File(h5file, 'r')
  tmp = read_displacement(f)
  count = 0
  for j in range((l+1)**2):
    R[i, count] = tmp[j][0]
    R[i, count + 1] = tmp[j][1]
    count += 2
  f.close()
  if over:
    # oversampling: matrix C needs to be created by the values of the interior corner nodes 
    if l % rep  != 0:
      print("Warning: Res divided by eps should be an integer. Otherwise rounding errors could occur.")
    if not args.triangle:
      #0-based node numbers of the corners of the center element
      coord = [int((l) / rep * (l+1.) + l/rep + 0.5 + 1e-6),int((l) / rep * (l+1.) + l/rep + 0.5 + 1e-6) + int(l/rep + 0.5 + 1e-6),int(2.*l/rep *(l+1) + 2.*l/rep + 0.5 + 1e-6),int(2.*l/rep *(l+1) + 1.*l/rep + 0.5 + 1e-6)]
      count = 0
      for j in range(int(dof/dim)):
        C[i,count] = tmp[coord[j]][0]
        C[i,count+1] = tmp[coord[j]][1]
        count += 2
    else:
      print("Warning: Oversampling is not implemented for triangles, yet.")     
  return R
parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction", type=int)
parser.add_argument("dimension", help="Dimension of the problem", type=int, default=2)
parser.add_argument("res", help="resolution of the cell problem in x-direction", type=int)
parser.add_argument("folder",help="specify the folder of the h5 files")
parser.add_argument("--mesh", help="mesh file")
parser.add_argument("--msfem", help="msfem basis functions are evaluated, insert E-modulus and Poisson ratio", default="")
parser.add_argument("--triangle", help="triangle msfem elements on/off", default="")

# not necessary normally, just experimental
parser.add_argument("--force_msfem",help="option calculates the force catalogue for MSFEM")

parser.add_argument("--sparse",help ="sparse mesh is used for msfem calculations")

#for debugging
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")

parser.add_argument("--big", help="specify number of cfs.rel runs in parallel, if option is turned on, mtx files and vec files are not saved")
parser.add_argument("--epsilon", help="number of frames/crosses in the cell problem", type=int, default=1)
parser.add_argument("--oversampling", help="name of the mesh with size minres/epsilon including only one base cell")
parser.add_argument("--penalization", help="creates a penalized material catalogue in the interval [0, 1/steps_p], step_p has to be given",type=int)

# skip reading .info.xml for parameter at bound 0 or 1 and add trivial coefficients eps or +1 instead
parser.add_argument("--skip_bounds", help="skip searching for parameter values 0 or 1 and take trivial (1e-6 or 1) coeffs", action='store_true', default=False)



args = parser.parse_args()
getcontext().prec = 16

steps = args.stp
folder = args.folder
dim = args.dimension
rep = args.epsilon
if dim == 2:
  if args.msfem:
    # setup for MSFEM element matrix data file in 2D
    n = args.msfem.split(',')
    E = float(n[0])
    nu = float(n[1])
    if args.triangle:
      # triangle elements for MSFEM
      dof = 6
      filename = "detailed_stats_" + str(folder)
    else:
      # quadrilateral elements for MSFEM
      dof = 8
      filename = "detailed_stats_" + str(folder)
    if args.force_msfem:
      filename = "force_msfem" + str(folder)
    out = open(filename, "w")
    out.write(str(1. / args.stp) + '  ' + str(1. / args.stp) + ' ')    
    if not args.force_msfem:
      for i in range((dof + 1) * dof / 2):
        out.write('0.0  ')
      out.write('\n\n')
  else:
    # setup for homogenization tensor file in 2D
    filename = "detailed_stats_" + str(folder)
    out = open(filename, "w")
    out.write('  ' + str(steps) + '   ' + str(steps) + '   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 \n')
    if args.big:
      joblist = ()
      # Read homogenized material tensors from cell problems in 2D and create detailed_stats table
      for line in open(folder+"/jobs", 'r'):
        joblist += ((line.strip()).split(), )
      pwd = os.path.dirname(os.path.abspath(__file__))
      os.chdir(str(pwd)+'/'+str(folder))
      # start calculation of the Homogenization cell problems, args.big specifies the number of problems run at the same time
      submit_job_max_len(joblist, max_processes=int(args.big))
      os.chdir(str(pwd))
elif dim == 3:
  # setup for homogenization tensor file in 3D
  filename = "detailed_stats_" + str(folder)
  out = open(filename, "w")
  out.write('  ' + str(steps) + '   ' + str(steps) + '  ' + str(steps) + '   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00    0.000000e+00 0.000000e+00 0.000000e+00\n')
  filename = "detailed_stats_vol" + str(folder)
  out_vol = open(filename, "w")
  out_vol.write('  ' + str(steps) + '   ' + str(steps) + '  ' + str(steps) + '   0.000000e+00\n')
  if args.big:
    joblist = ()
    # Read homogenized material tensors from cell problems in 2D and create detailed_stats table
    for line in open(folder+"/jobs", 'r'):
      joblist += ((line.strip()).split(), )
    pwd = os.path.dirname(os.path.abspath(__file__))
    os.chdir(str(pwd)+'/'+str(folder))
    # start calculation of the Homogenization cell problems, args.big specifies the number of problems run at the same time
    submit_job_max_len(joblist, max_processes=int(args.big))
    os.chdir(str(pwd))  
  
if dim == 2:
  x = 0
  # loops over all discretized design variables
  while x < steps + 1:
    if args.msfem:
      y = 0
      while y < steps + 1:
        if args.penalization:
          steps_p = args.penalization
          x_tmp = x
          y_tmp = y
          x = float(x) / steps_p if x > 0 else 0
          y = float(y)/ steps_p if y > 0 else 0
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))  
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        if args.sparse and x == 0 and y == 0:
          continue
        if args.big:
          index = 0
          pwd = os.path.dirname(os.path.abspath(__file__))
          os.chdir(str(pwd)+'/'+str(folder))
          joblist = ()
          print(os.path.dirname(os.path.abspath(__file__)))
          # loop for all MSFEM cell problems
          for i in range(dof):
            if i % 2 == 0:
              joblist += (('cfs.rel -m ' + str(args.mesh) + ' -x ' + str(x) + "-" + str(y) + "_msfem.dens.xml "+ str(x) + "-" + str(y) + '_msfem' + str(index) + '_x').split(), )
            else:
              joblist += (('cfs.rel -m ' + str(args.mesh) + ' -x ' + str(x) + "-" + str(y) + "_msfem.dens.xml "+ str(x) + "-" + str(y) + '_msfem' + str(index) + '_y').split(), )
              index += 1
          if args.oversampling:
            joblist += (('cfs.rel -m ' + str(args.oversampling) + ' -x ' + str(x) + "-" + str(y) + "_msfem_oversample.dens.xml "+ str(x) + "-" + str(y) + '_msfem').split(), )
          else:
            joblist += (('cfs.rel -m ' + str(args.mesh) + ' -x ' + str(x) + "-" + str(y) + "_msfem.dens.xml "+ str(x) + "-" + str(y) + '_msfem').split(), )
          #print ('cfs.rel -m ' + str(args.mesh) + ' -x ' + str(x) + "-" + str(y) + "_msfem.dens.xml "+ str(x) + "-" + str(y) + '_msfem' + str(index) + '_x').split()
          #joblist = ("sleep 5".split(),)
          #joblist = (['cfs.rel', '-m', '/home/snthgues/meshes/10_ref.mesh', '-x', '0-0_msfem.dens.xml', '0-0_msfem3_y'],)
          
          # start calculation of the MSFEM cell problems, args.big specifies the number of problems run at the same time
          submit_job_max_len(joblist, max_processes=int(args.big))
          os.chdir(str(pwd))
        index = 0
        h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem0_x.h5"
        f = h5py.File(h5file, 'r')
        if args.force_msfem:
          mesh = create_mesh_from_hdf5(f, ['mech'],['bottom','top','left','right'])
        f.close()
        R = np.matrix(np.zeros((dof, dim * (args.res+1)**2)))
        
        #oversampling: get matrix C for fixing the macroscopic boundary conditions 
        if args.oversampling:
          C = np.matrix(np.zeros((dof,dof)))
        
        for i in range(dof):
          if i % 2 == 0:
            h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem" + str(index) + "_x.h5"
          else:
            h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem" + str(index) + "_y.h5"
          if args.sparse and x == 0 and y == 0:
            continue
          # calculate restriction matrix for MSFEM 2D
          R = get_R(h5file,i,R,dof,dim,args.res,True if args.oversampling else False, rep if args.oversampling else None, C if args.oversampling else None)
          if i % 2 != 0:
            index += 1    
               
        if args.force_msfem:
          print('WARNING: force has to be set manually in the code')
          f = force(dim,mesh,R,'bottom')
          out.write(str(x) + ' ' + str(y) + ' ')
          for i in range(dof):
            out.write(str(f[i]) + ' ')
          out.write('\n')
          #scipy.io.savemat('R.mat', mdict={'R': R})
          print('Calculation of force for st1,st2 = ' + str(x) + '\t' + str(y) + ' done \n')          
        else:
          if args.oversampling:
            #print "matrix "
            #for i in range(dof):
            #  for j in range(dof):
            #    print str(C[i, j]) + ' '
            #  print "\n"
            
            # inverting is ok, C is only a small matrix, 8x8 or 6x6
            Cinv = np.linalg.inv(C)
            R = np.dot(Cinv,R)
            # get macroscopic MSFEM element matrix 
            A = CalcElementMat(R,folder+'/'+str(x)+'-'+str(y)+'_msfem_0.mtx',dof,dim,args.res,rep) #if args.big else '_msfem_iter_0_excite_0.mtx',overR)
          else:
            # get macroscopic MSFEM element matrix 
            A = CalcElementMat(R,folder+'/'+str(x)+'-'+str(y)+'_msfem_0.mtx',dof,dim) #if args.big else '_msfem_iter_0_excite_0.mtx')
          out.write(str(x) + ' ' + str(y) + ' ')
          for i in range(dof):
            for j in range(i, dof):
              out.write(str(A[i, j]) + ' ')
          out.write('\n')
          if args.big:
            execute('rm -f '+folder+'/'+str(x)+'-'+str(y)+'_msfem_0.mtx')
            execute('rm -f '+folder+'/'+str(x)+'-'+str(y)+'_msfem_0.vec')
            #print 'rm -f '+folder+'/'+str(x)+'-'+str(y)+'_msfem_0.mtx'
          # scipy.io.savemat('A.mat', mdict={'A': A})
          # for i in range(dof):
          #    print '  ' + str(A[i, :])
          # scipy.io.savemat('R.mat', mdict={'R': R})
          print('Calculation of element matrix for st1,st2 = ' + str(x) + '\t' + str(y) + ' done \n')          
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          if args.penalization:
            x = x_tmp
            y = y_tmp
          y += 1 
    else:
      y = 0
      while y < steps + 1:
        if args.penalization:
          steps_p = args.penalization
          x_tmp = x
          y_tmp = y
          x = float(x) / steps_p if x > 0 else 0
          y = float(y)/ steps_p if y > 0 else 0
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))  
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        problem = (str(x_tmp) +"-" + str(y_tmp) + "_" + str(1./steps_p)) if args.penalization else str(x) + "-" + str(y)  
        infoxml = str(folder) + "/" + problem + ".info.xml"
        # print infoxml
        if os.path.isfile(infoxml):      
          doc = lxml.etree.parse(infoxml, lxml.etree.XMLParser(remove_comments=True))
          #xml = doc.xpathNewContext()
          # complex values!  
          print(infoxml + ' -> ')
          matrix = xpath(doc, "//homogenizedTensor/tensor/real/text()")	
          res = list(map(float, matrix.split())) # convert list with string elements to list with float elements	
          ts = np.asarray(res)            # convert list to array
          out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + ts[0] + ' ' + ts[1] + ' ' + ts[4] + ' ' + ts[8] + '\n')
          if x != y:
            out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
        else:
          print('file ' + infoxml + ' not found')
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          if args.penalization:
            x = x_tmp
            y = y_tmp
          y += 1
    if not args.design:
      x += 1
elif dim == 3:
  # Read homogenized material tensors from cell problems in 3D and create detailed_stats table
  x = 0
  while x < steps + 1:
    y = 0
    while y < steps + 1:
      z= 0
      while z < steps + 1:
        vol = None
        sampl = np.array([x,y,z])
        # only on parameter is 0
        if False:
          print("skip ",x,y,z)
#           z += 1
#           continue
        else:  
          if args.penalization:
            steps_p = args.penalization
            x_tmp = x
            y_tmp = y
            z_tmp = z
            x = float(x) / steps_p if x > 0 else 0
            y = float(y)/ steps_p if y > 0 else 0
            z = float(z)/ steps_p if z > 0 else 0
            if args.design:
              tmp = args.design.split(',')
              x = float(steps * float(tmp[0]))
              y = float(steps * float(tmp[1]))
              z = int(steps * float(tmp[2]))  
          elif args.design:
            tmp = args.design.split(',')
            x = int(steps * float(tmp[0]))
            y = int(steps * float(tmp[1]))
            z = int(steps * float(tmp[2]))
          infoxml = str(folder) + "/" + str(x) + "-" + str(y) + "-" + str(z) + ".info.xml"
          if np.sum(sampl == 0) == 2 and x == 0:
            if y == 0:
              infoxml = str(folder) + "/" + str(z) + "-0-0.info.xml"
            else:
              assert(z == 0)
              infoxml = str(folder) + "/" + str(y) + "-0-0.info.xml"
          if np.sum(sampl == 0) == 1 and x == 0:
            infoxml = str(folder) + "/" + str(z) + "-" + str(y) +  "-0.info.xml"
          elif np.sum(sampl == 0) == 1 and y == 0:        
            infoxml = str(folder) + "/" + str(x) + "-" + str(z) +  "-0.info.xml"
          # print infoxml
          if os.path.isfile(infoxml):      
            doc = lxml.etree.parse(infoxml, lxml.etree.XMLParser(remove_comments=True))
            #xml = doc.xpathNewContext()
            # complex values!  
            print(infoxml + ' -> ')
            matrix = xpath(doc, "//homogenizedTensor/tensor/real/text()")    
            res = list(map(float, matrix.split())) # convert list with string elements to list with float elements    
            ts = np.asarray(res)
            print(infoxml + ' -> '+ str(ts))
            
            if np.sum(sampl == 0) > 0:
              assert(dim == 3)
              tens = ts.reshape(6,6) 
              Q = None # rotation matrix
              
              # two params = 0, find out which one is != 0 in order to rotate tensor properly
              # we only have homogenized tensors for s1 != and s2=s3=0
              if np.sum(sampl == 0) == 2:
                if  x == 0:
                  if y == 0:
                    # x1 -> z1: rotate around y-axis by 90° 
                    Q = matviz_rot.get_rot_6x6(0,np.pi/2,0)
                  else:
                    assert(z == 0)  
                    # x1 -> y1: rotate around z-axis by 90°
                    Q = matviz_rot.get_rot_6x6(0,0,np.pi/2)
                else: # x != 0
                  Q = np.eye(6)
                
              # one param s_{1,2,3} is 0
              else:
                assert(np.sum(sampl == 0) == 1)
                ts_vals = []
                if x == 0:
                  print("\n here x:",x,y,z)
                  # s1 =0 , s2 > 0, s3 > 0
                  # keep y1 and rotate x1 onto z1
                  
                  # x1 -> z1: rotate around y-axis by 90° 
                  Q = matviz_rot.get_rot_6x6(0,np.pi/2,0)
                elif y == 0:
                  print("\n here y:",x,y,z)
                  # s1 > 0, s2 = 0, s3 > 0
                  # keep x1 and rotate y1 onto z1
                  # x1 -> z1: rotate around x-axis by 90°
                  Q = matviz_rot.get_rot_6x6(np.pi/2,0,0)
                else: 
                  print("\n here else:",x,y,z)
                  #s1 > 0, s2 > 0, s3 = 0
                  # just take identity matrix -> no rotation
                  Q = np.eye(6)
                
              assert(Q is not None)
              rot_tens =  np.dot(Q, np.dot(tens, Q.transpose()))
              out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + str(rot_tens[0,0]) + ' ' + str(rot_tens[1,0]) + ' ' + str(rot_tens[2,0]) + ' ' + str(rot_tens[1,1]) + ' ' + str(rot_tens[1,2]) + ' ' + str(rot_tens[2,2]) + ' ' + str(rot_tens[3,3]) + ' ' + str(rot_tens[4,4]) + ' ' + str(rot_tens[5,5]) + "\n")  
              # we have to use .density.xml files, thus read volume from optimization part
              vol = xpath(doc, "//iteration/@volume")
              out_vol.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + vol + '\n')    
                  
            else:  # neither s1, s2 or s3 = 0
              out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + str(ts[0]) + ' ' + str(ts[1]) + ' ' + str(ts[2]) + ' ' + str(ts[7]) + ' ' + str(ts[8]) + ' ' + str(ts[14]) + ' ' + str(ts[21]) + ' ' + str(ts[28]) + ' ' + str(ts[35]) + '\n')
              vol = xpath(doc, "//domain/@structure_volume")
              
              if vol != "cannot_determine":
                vol = float(vol)
                out_vol.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + str(vol) + '\n')
          
          # if x != y:
          #  out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
          else:
            print('file ' + infoxml + ' not found')
            
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
          z = steps + 1
        else:
          if args.penalization:
              x = x_tmp
              y = y_tmp
              z = z_tmp
          z += 1
      if not args.design:
        y += 1
    if not args.design:
      x += 1
out.close()
out_vol.close()


# rotate 3d tensor around one axis by given angle
# @tensor: 6x6 numpy array
# @axis: 0,1 or 2
# @angle: rotation angle in degrees
# returns rotated tensor
def rotate_3d_tensor(tensor,axis,angle):
  assert(tensor.shape == (6,6))
  assert(-1e-6 <= angle <= 360+1e-6)
  assert(0 <= axis <= 2)
  rad = np.radians(angle)
  Q = None
  if axis == 0:
    Q = matviz_rot.get_rot_6x6(rad,0,0)
  elif axis == 1:
    Q = matviz_rot.get_rot_6x6(0,rad,0)
  else: # axis == 2
    Q = matviz_rot.get_rot_6x6(0,0,rad)
    
  assert(Q is not None)
  print("rotated tensor by ", angle, " around axis ",axis)
  
  return np.dot(Q, dot(tensor, Q.transpose()))


# FUNCTIONS NOT USED CURRENTLY

def force(dim,mesh,R,region):
  # OLD code: scaled force for MSFEM
  f = scipy.sparse.lil_matrix((dim * len(mesh.nodes), 1))
  # volume force
  #for i in range(len(mesh.nodes)):
  #  f[2*i+1,0] = -1./(51.*51.)
  # boundary force
  for i in range(len(mesh.bc)):
    if mesh.bc[i][0] == region:
      for j in range(len(mesh.bc[i][1])):
        index = mesh.bc[i][1][j]
        f[2*index+1,0] = -1/51.
  return R*f
  
def Kelem(E, nu, dx, dy, dof, dim, points=None):
  # Calculate FEM element matrix
  if dof == 8:
    b = [-math.sqrt(1. / 3.), math.sqrt(1. / 3.)]
  elif dof == 6:
    # are not used for linear triangle elements
    b = [0.5 - 0.5 * math.sqrt(1. / 3.), 0.5 + 0.5 * math.sqrt(1. / 3.)]
    
  C = E / (1. - nu ** 2) * np.matrix([[1., nu, 0.], [nu, 1., 0.], [0., 0., (1. - nu) / 2. ]]);
  KE = np.matrix(np.zeros((dof, dof)))
  if dof == 8:
    for i in range(dim):
      for j in range(dim):
        (B, detJ) = CalcBMat(b[i], b[j], dx, dy, dof, points)
        KE += B.T * C * B
  elif dof == 6:
    (B, detJ) = CalcBMat(1., 1., dx, dy, dof, points)
    KE = 0.5 * B.T * C * B
  # only valid for regular, uniform grids
  # print KE
  KE *= detJ
  return KE 
def CalcBMat(x, y, dx, dy, dof, points=None):
  # Calculate B-Matrix for transformation
  
  #quadrilateral elements 2D
  if dof == 8:
    # Derivatives of basis functions in reference element
    N1dx = -0.25 * (1. - y)
    N2dx = 0.25 * (1. - y)
    N3dx = 0.25 * (1. + y)
    N4dx = -0.25 * (1. + y)
    N1dy = -0.25 * (1. - x)
    N2dy = -0.25 * (1. + x)
    N3dy = 0.25 * (1. + x)
    N4dy = 0.25 * (1. - x)
    
    # matrix B in reference element coordinates
    # c1 = 2. / dx
    # c2 = 2. / dy
        # B = np.matrix([[c1 * N1dx, 0., c1 * N2dx, 0., c1 * N3dx, 0., c1 * N4dx, 0.],
    #         [0., c2 * N1dy, 0., c2 * N2dy, 0., c2 * N3dy, 0., c2 * N4dy],
    #         [c2 * N1dy, c1 * N1dx, c2 * N2dy, c1 * N2dx, c2 * N3dy, c1 * N3dx, c2 * N4dy, c1 * N4dx]])
    # B = np.matrix([[c1 *N1dx,c1*N2dx,c1*N3dx,c1*N4dx,0.,0.,0.,0.],[0.,0.,0.,0.,c2*N1dy,c2*N2dy,c2*N3dy,c2*N4dy],[c2*N1dy,c2*N2dy,c2*N3dy,c2*N4dy,c1*N1dx,c1*N2dx,c1*N3dx,c1*N4dx]])
    J11 = 0.25 * (-(1. - y) * points[0, 0] + (1. - y) * points[1, 0] + (1. + y) * points[2, 0] - (1. + y) * points[3, 0])
    J12 = 0.25 * (-(1. - y) * points[0, 1] + (1. - y) * points[1, 1] + (1. + y) * points[2, 1] - (1. + y) * points[3, 1])
    J21 = 0.25 * (-(1. - x) * points[0, 0] - (1. + x) * points[1, 0] + (1. + x) * points[2, 0] + (1. - x) * points[3, 0])
    J22 = 0.25 * (-(1. - x) * points[0, 1] - (1. + x) * points[1, 1] + (1. + x) * points[2, 1] + (1. - x) * points[3, 1])
    detJ = J11 * J22 - J12 * J21
    # print 'detJ = ' + str(detJ)
    # B = (1. / detJ) * np.matrix([[J22 * N1dx - J12 * N1dy, J22 * N2dx - J12 * N2dy, J22 * N3dx - J12 * N3dy, J22 * N4dx - J12 * N4dy, 0., 0., 0., 0.], [0., 0., 0., 0., -J21 * N1dx + J11 * N1dy, -J21 * N2dy + J11 * N2dy, -J21 * N3dx + J11 * N3dy, -J21 * N4dx + J11 * N4dy], [-J21 * N1dx + J11 * N1dy, -J21 * N2dx + J11 * N2dy, -J21 * N3dx + J11 * N3dy, -J21 * N4dx + J11 * N4dy, J22 * N1dx - J12 * N1dy, J22 * N2dx - J12 * N2dy, J22 * N3dx - J12 * N3dy, J22 * N4dx - J12 * N4dy]])
    B = (1. / detJ) * np.matrix([[J22 * N1dx - J12 * N1dy, 0., J22 * N2dx - J12 * N2dy, 0., J22 * N3dx - J12 * N3dy, 0., J22 * N4dx - J12 * N4dy, 0.], [0., -J21 * N1dx + J11 * N1dy, 0., -J21 * N2dy + J11 * N2dy, 0., -J21 * N3dx + J11 * N3dy, 0., -J21 * N4dx + J11 * N4dy], [-J21 * N1dx + J11 * N1dy, J22 * N1dx - J12 * N1dy, -J21 * N2dx + J11 * N2dy, J22 * N2dx - J12 * N2dy, -J21 * N3dx + J11 * N3dy, J22 * N3dx - J12 * N3dy, -J21 * N4dx + J11 * N4dy, J22 * N4dx - J12 * N4dy]])
  # triangluar elements
  elif dof == 6:
    N1dx = -1.
    N1dy = -1.
    N2dx = 1.
    N2dy = 0.
    N3dx = 0.
    N3dy = 1.
    x21 = points[1, 0] - points[0, 0]
    x31 = points[2, 0] - points[0, 0]
    x13 = points[0, 0] - points[2, 0]
    x32 = points[2, 0] - points[1, 0]
    y12 = points[0, 1] - points[1, 1]
    y31 = points[2, 1] - points[0, 1]
    y23 = points[1, 1] - points[2, 1]
    y21 = points[1, 1] - points[0, 1]
    # print 'points ' + str(points.shape)
    J = np.matrix([[-points[0, 0] + points[1, 0], -points[0, 0] + points[2, 0]], [-points[0, 1] + points[1, 1], -points[0, 1] + points[2, 1]]])
    # print 'J= ' + str(J)
    detJ = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    # print 'triangle detJ =' + str(detJ)
    # B = (1. / detJ) * np.matrix([[y23, y31, y12, 0., 0., 0.], [0., 0., 0., x32, x13, x21], [x32, x13, x21, y23, y31, y12]])
    B = (1. / detJ) * np.matrix([[y23, 0., y31, 0., y12, 0.], [0., x32, 0., x13, 0., x21], [x32, y23, x13, y31, x21, y12]])
    # print 'B= ' + str(B)
  else:
    print('error: degree of freedom not defined')
  return (B, detJ)

def Kglob(E, nu, dx, dy, dim, dof, mesh):
  # calculate global FEM stiffness matrix
  # E material matrix, nu poisson's ratio; dx,dy length of the reference cell
  K = scipy.sparse.lil_matrix((dim * len(mesh.nodes), dim * len(mesh.nodes))) 
  points = np.matrix(np.zeros((dof / dim, dim)))
  for t in range(len(mesh.elements)):
    e = mesh.elements[t]
    edof = np.array([0] * dof)
    count = 0
    for k in range(dof / dim):
      edof[count] = 2 * e.nodes[k]
      edof[count + 1] = 2 * (e.nodes[k]) + 1
      count += 2
      points[k, 0] = mesh.nodes[(e.nodes[k])][0]
      points[k, 1] = mesh.nodes[(e.nodes[k])][1]
    KE = Kelem(E, nu, dx, dy, dof, dim, points)
    for i in range(dof):
      for j in range(dof):
        K[edof[i], edof[j]] += e.density * KE[i, j]
  # scipy.io.savemat('KE.mat', mdict={'KE': KE})
  return K

# # This rotates a 3*3 2D tensor via the third direction. As in Richter and CFS
def get_rot_3x3(angle): 
  R = numpy.zeros((2, 2))

  R[0][0] = cos(angle)
  R[0][1] = sin(angle)
  R[1][0] = -sin(angle)
  R[1][1] = cos(angle)
  
  Q = numpy.zeros((3, 3))
  
  Q[0][0] = R[0][0] * R[0][0];
  Q[0][1] = R[0][1] * R[0][1];
  Q[0][2] = 2.0 * R[0][0] * R[0][1];

  Q[1][0] = R[1][0] * R[1][0];
  Q[1][1] = R[1][1] * R[1][1];
  Q[1][2] = 2.0 * R[1][0] * R[1][1];

  Q[2][0] = R[0][0] * R[1][0];
  Q[2][1] = R[0][1] * R[1][1];
  Q[2][2] = R[0][0] * R[1][1] + R[0][1] * R[1][0];
  
  return Q
