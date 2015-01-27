#!/usr/bin/env python

import sys, os, string
import numpy as np
import scipy.io
from numpy import dot
from numpy import sin
from numpy import cos
from numpy import sqrt
from lxml import etree
import argparse
import scipy, scipy.sparse, scipy.sparse.linalg 
from hdf5_tools import *
from optimization_tools import *
from mesh_tool import *
from decimal import *

def force(dim,mesh,R,region):
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
    print 'error: degree of freedom not defined'
  return (B, detJ)

def Kglob(E, nu, dx, dy, dim, dof, mesh):
  K = scipy.sparse.lil_matrix((dim * len(mesh.nodes), dim * len(mesh.nodes))) 
  points = np.matrix(np.zeros((dof / dim, dim)))
  for t in range(len(mesh.elements)):
    e = mesh.elements[t]
    edof = np.array([0] * dof)
    count = 0
    for k in range(dof / dim):
      # edof[k] = e.nodes[k]
      # edof[k + dof / dim] = e.nodes[k] + len(mesh.nodes)
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

def CalcElementMat(dim, dof, E, nu, R, mesh,Kmat):
  dx = 1.
  dy = 1.
  #K = Kglob(E, nu, dx, dy, dim, dof, mesh)
  K = scipy.io.mmread(Kmat)
  # scipy.io.savemat('K.mat', mdict={'K': K})
  #scipy.io.savemat('tmp.mat', mdict={'K': tmp})
  elemMat = R * K * R.T
  return elemMat

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
parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction", type=int)
parser.add_argument("dimension", help="Dimension of the problem", type=int, default=2)
parser.add_argument("folder",help="specify the folder of the h5 files")
parser.add_argument("--msfem", help="msfem basis functions are evaluated, insert E-modulus and Poisson ratio", default="")
parser.add_argument("--triangle", help="triangle msfem elements on/off", default="")
parser.add_argument("--force_msfem",help="option calculates the force catalogue for MSFEM")
parser.add_argument("--sparse",help ="sparse mesh is used for msfem calculations")
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")

args = parser.parse_args()
getcontext().prec = 16

steps = args.stp
folder = args.folder
dim = args.dimension
if dim == 2:
  if args.msfem:
    n = args.msfem.split(',')
    E = float(n[0])
    nu = float(n[1])
    if args.triangle:
      dof = 6
      filename = "detailed_stats_" + str(folder)
    else:
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
    filename = "detailed_stats_" + str(folder)
    out = open(filename, "w")
    out.write('  ' + str(steps) + '   ' + str(steps) + '   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 \n')
elif dim == 3:
  filename = "detailed_stats_" + str(folder)
  out = open(filename, "w")
  out.write('  ' + str(steps) + '   ' + str(steps) + '  ' + str(steps) + '   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00    0.000000e+00 0.000000e+00 0.000000e+00\n')
  

  # out.write('# x   y 11           12           22           33\n')
if dim == 2:
  #if args.msfem:
    #if dof == 8:
      #h5file = str(folder) + "/1-0_msfem0_x.h5"
    #elif dof == 6:
      #h5file = str(folder) +"/1-0_msfem0_x.h5"
    #f = h5py.File(h5file, 'r')
    #mesh = create_mesh_from_hdf5(f, ['mech'],[])
    #f.close()
  x = 0
  while x < steps + 1:  #x = 5
    if args.msfem:
      y = 0
      while y < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        #y = 0
        # if True:
        index = 0
        h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem0_x.h5"
        if args.sparse and x == 0 and y == 0:
          continue
        f = h5py.File(h5file, 'r')
        mesh = create_mesh_from_hdf5(f, ['mech'],['bottom','top','left','right'])
        R = np.matrix(np.zeros((dof, dim * len(mesh.nodes))))
        f.close()
        for i in range(dof):
          if i % 2 == 0:
            h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem" + str(index) + "_x.h5"
          else:
            h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem" + str(index) + "_y.h5"
          if args.sparse and x == 0 and y == 0:
            continue
          f = h5py.File(h5file, 'r')
          tmp = read_displacement(f)
          # if dof == 8:
          #  index1 = [0, 4, 1, 5, 2, 6, 3, 7]
          # elif dof == 6:
          #  index1 = [0, 3, 1, 4, 2, 5]
          # for j in range(len(mesh.nodes)):
          #    R[index1[i], j] = tmp[j][0]
          #    R[index1[i], j + len(mesh.nodes)] = tmp[j][1]
          count = 0
          for j in range(len(mesh.nodes)):
            R[i, count] = tmp[j][0]
            R[i, count + 1] = tmp[j][1]
            count += 2
          if i % 2 != 0:
            index += 1
          f.close()
        # just a sample density xml file
        #if dof == 8:
        #  h5file = str(folder) + "/" + str(x) + "-" + str(y) + "_msfem0_x.h5"
        #else:
        #  print 'Option not implemented yet, for triangles'
        #if x== 0 and y == 0:
        #  continue
        if args.force_msfem:
          print 'WARNING: force has to be set manually in the code'
          f = force(dim,mesh,R,'bottom')
          out.write(str(x) + ' ' + str(y) + ' ')
          for i in range(dof):
            out.write(str(f[i]) + ' ')
          out.write('\n')
          #scipy.io.savemat('R.mat', mdict={'R': R})
          print 'Calculation of force for st1,st2 = ' + str(x) + '\t' + str(y) + ' done \n'          
        else:
          A = CalcElementMat(dim, dof, E, nu, R, mesh,folder+'/'+str(x)+'-'+str(y)+'_msfem_iter_0_excite_0.mtx')
          out.write(str(x) + ' ' + str(y) + ' ')
          for i in range(dof):
            for j in range(i, dof):
              out.write(str(A[i, j]) + ' ')
          out.write('\n')
          # scipy.io.savemat('A.mat', mdict={'A': A})
          # for i in range(dof):
          #    print '  ' + str(A[i, :])
          # scipy.io.savemat('R.mat', mdict={'R': R})
          print 'Calculation of element matrix for st1,st2 = ' + str(x) + '\t' + str(y) + ' done \n'          
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          x += 1
          y += 1   

    else:
      y = 0
      while y < x + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))  
        infoxml = str(folder) + "/" + str(x) + "-" + str(y) + ".info.xml"
        # print infoxml
        if os.path.isfile(infoxml):      
          doc = libxml2.parseFile(infoxml)
          xml = doc.xpathNewContext()
          # complex values!  
          print infoxml + ' -> '
          tensortext = xpath(xml, "//homogenizedTensor/tensor/real")	
          print infoxml + ' -> ' + tensortext	
          ts = tensortext.split()
          # print ts[0]
          # print ts[1]
          # print ts[4]
          # print ts[8]
    
          out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + ts[0] + ' ' + ts[1] + ' ' + ts[4] + ' ' + ts[8] + '\n')
          if x != y:
            out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
    
          if False:
            tensor = numpy.zeros((3, 3))
            tensor[0, 0] = ts[0]
            tensor[0, 1] = ts[1]
            tensor[0, 2] = ts[2]
            tensor[1, 0] = ts[3]
            tensor[1, 1] = ts[4]
            tensor[1, 2] = ts[5]
            tensor[2, 0] = ts[6]
            tensor[2, 1] = ts[7]
            tensor[2, 2] = ts[8]
            Q = get_rot_3x3(numpy.pi / 2.0)
            rott = dot(Q, dot(tensor, Q.transpose()))
            out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' {:.6e} {:.6e} {:.6e} {:.6e}\n'.format(rott[0, 0], rott[0, 1], rott[1, 1], rott[2, 2]))
        else:
          print 'file ' + infoxml + ' not found'
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          x += 1
          y += 1
elif dim == 3:
  x = 0
  while x < steps + 1:
    y = 0
    while y < steps + 1:
      z= 0
      while z < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
          z = int(steps * float(tmp[2]))
        infoxml = str(folder) + "/" + str(x) + "-" + str(y) + "-" + str(z) + ".info.xml"
        # print infoxml
        if os.path.isfile(infoxml):      
          doc = libxml2.parseFile(infoxml)
          xml = doc.xpathNewContext()
          # complex values!  
          print infoxml + ' -> '
          tensortext = xpath(xml, "//homogenizedTensor/tensor/real")  
          print infoxml + ' -> ' + tensortext  
          ts = tensortext.split()
          # print ts[0]
          # print ts[1]
          # print ts[4]
          # print ts[8]
          out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + ts[0] + ' ' + ts[1] + ' ' + ts[2] + ' ' + ts[7] + ' ' + ts[8] + ' ' + ts[14] + ' ' + ts[21] + ' ' + ts[28] + ' ' + ts[35] + '\n')
          # if x != y:
          #  out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
        else:
          print 'file ' + infoxml + ' not found'
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
          z = steps + 1
        else:
          x += 1
          y += 1
          z += 1
# print "number of tensors computed = " + str(count) + " / " + str(countall) + " (= {0:2.2f}%)".format(100.0*float(count)/float(countall))
# out.write("number of tensors computed = " + str(count) + "\n")

# out.write(str(steps) + ' ' + str(steps) + ' 5.698006e+00 1.994302e+00 5.698006e+00 1.851852e+00')

out.close()
