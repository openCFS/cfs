#!/usr/bin/env python
# this ia a tool for spaghetti visualization and calculation.
import numpy as np
from numpy.linalg import norm
import os
import cfs_utils as ut
from itertools import product
from multiprocessing import Pool
import time as ti

# for gradient check
import scipy.optimize as sciopt
import scipy.sparse as scispr


# interactive
if __name__ == '__main__':
  import optimization_tools as ot
  import argparse
  # noinspection PyUnresolvedReferences
  import vtkmodules.vtkInteractionStyle
  # noinspection PyUnresolvedReferences
  import vtkmodules.vtkRenderingOpenGL2
  from vtkmodules.vtkCommonColor import vtkNamedColors
  from vtkmodules.vtkCommonTransforms import vtkTransform
  from vtkmodules.vtkInteractionWidgets import vtkOrientationMarkerWidget
  from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
  from vtkmodules.vtkCommonCore import vtkPoints
  from vtkmodules.vtkIOXML import vtkXMLPolyDataWriter
  from vtkmodules.vtkCommonDataModel import (
      vtkCellArray,
      vtkPolyData
  )
  from vtkmodules.vtkFiltersCore import vtkGlyph3D
  from vtkmodules.vtkFiltersSources import (
    vtkCylinderSource,
    vtkSphereSource,
    vtkParametricFunctionSource,
    vtkArrowSource
  )
  from vtkmodules.vtkCommonComputationalGeometry import vtkParametricTorus
  from vtkmodules.vtkRenderingCore import (
      vtkActor,
      vtkPolyDataMapper,
      vtkRenderWindow,
      vtkRenderWindowInteractor,
      vtkRenderer
  )
  from vtk import (
    vtkOBJExporter,
    vtkX3DExporter
  )

# to conveniently access global values
class Global:
  # the default value are for easy debugging via import spaghetti
  def __init__(self):
    self.shapes = []         # array of Spaghetti
    self.rhomin = 1e-4      # void
    self.boundary = 'poly'   # up to now only 'poly' and 'linear'
    self.p = 10              # penalty parameter
    self.transition = 0.05   # parameter for boundary: 2*h
    self.h = 0.5*self.transition
    self.n = np.array([40,40,40],dtype=int)       # [nx, ny, nz]
    self.box_lower = np.array([0,0,0])
    self.box_upper = np.array([1,1,1])
    self.dx = (self.box_upper[0]-self.box_lower[0])/self.n[0]
    self.mid_lower = self.box_lower + 0.5*self.dx*np.ones((3))
    self.order = 1
    self.silent = False
    self.vtk_lists = False
    self.num_threads = int(os.environ.get('OMP_NUM_THREADS'))

  # total number of variables (not cached)
  def total(self):
    return sum([s.num_optvar for s in self.shapes])

  # give optimzation variables of all spaghetti as defined array such that we can easily differentiate
  # @return p_x,p_y,q_x,q_y,p,a_1,...a_n-1,p_x,p_y,q_x,q_y,p,a_1,...a_n-1,...
  def var_all(self):
    vars = []
    for i in range(len(self.shapes)):
      vars.extend(self.shapes[i].optvar())
    return vars

glob = Global()

# FE for integration list
class FE:
  def __init__(self, i,j,k):
    self.i = i
    self.j = j
    self.k = k
    self.shapes_idx_list = []
    self.shapes_full = []

# @settings dict of key/string from openCFS or from command line
# @design tupel with design names as strings, usually only 'density'
# @dict dictionary transparently given from the xml file to python
def cfs_init(settings, design, dict):

  glob.rhomin = float(settings['rhomin'])
  glob.rhomax = float(settings['rhomax'])-glob.rhomin
  glob.boundary = settings['boundary']
  glob.transition = float(settings['transition'])
  glob.h = .5*float(glob.transition)
  glob.combine = settings['combine']

  glob.box_lower = np.array(eval(settings['box_lower']))
  glob.box_upper = np.array(eval(settings['box_upper']))

  glob.n = np.array(eval(settings['n']),dtype=int)

  glob.dx = (glob.box_upper[0]-glob.box_lower[0])/glob.n[0]
  glob.mid_lower = glob.box_lower + 0.5*glob.dx*np.ones((3))

  if 'p' in dict:
    glob.p = dict['p']
  if 'vtk_lists' in dict:
    glob.vtk_lists = dict['vtk_lists']
  if 'silent' in dict:
    glob.silent = int(dict['silent'])

  glob.design = design
  glob.opts = dict


## set spaghetti. Initially create it, otherwise update
# @param id the 0-based index of the spaghetti, corresponds index within glob.shapes
# @param points list of points which are lists
# @param r are radii as list of size segments -1
# @param p and alpha are scalar 
def cfs_set_spaghetti(id, points, r, p, alpha):
  #print('cfs_set_spaghetti id=',id,'points=',points,'r=',r,'p=',p,'alpha=',alpha)
  
  assert id < int(len(glob.shapes) + 1) # we may add another spaghetti

  px = [v[0] for v in points]
  py = [v[1] for v in points]
  pz = [v[2] for v in points]

  if id >= len(glob.shapes):
    base = sum([len(s.optvar()) for s in glob.shapes])
    # def __init__(self,id,base,radius,P,Q,a,p):
    glob.shapes.append(Spaghetti(id, base, px,py,pz,r,p,alpha))
    print('cfs_set_spaghetti: create ', glob.shapes[-1], 'id=',id,'points=',points,'r=',r,'p=',p,'alpha=',alpha)
  else:
    glob.shapes[id].set(px,py,pz,r,p,alpha)
    if not glob.silent:
      print('cfs_set_spaghetti: update ', glob.shapes[id])

## return 1D-array with mapped design
def cfs_map_to_design():
  assemble_fe_field()
  start = ti.time()
  with Pool(glob.num_threads) as pool:
    rho = pool.map(integrate_fe, range(np.prod(glob.n)))
    pool.close()
    pool.join()
  #integrate_fe_field()
  if not glob.silent:
    print('time for integration:', ti.time()-start)
  start = ti.time()
  glob.rho = np.zeros(glob.n)
  glob.grad_rho = np.zeros((glob.n[0],glob.n[1],glob.n[2],glob.total()))
  for fe_num in range(np.prod(glob.n)):
    glob.rho[np.unravel_index(fe_num,glob.n)] = rho[fe_num][0]
    glob.grad_rho[np.unravel_index(fe_num,glob.n)] = rho[fe_num][1]
  if not glob.silent:
    print('time for unraveling:', ti.time()-start)
  if glob.vtk_lists:
    write_vtk_lists(glob.vtk_lists)
  return np.reshape(glob.rho, np.prod(glob.n), order='F')+glob.rhomin

def cfs_info_field_keys():
  if not glob.silent:
    print('Warning: cfs_info_field_keys only in dummy implementation')
  return []

# as we cannot create a numpy array in C (it should work but fails in reality) we get it here.
# it shall have the size of rho as a 1D array  
def cfs_get_drho_vector():
  size = len(glob.design) * np.prod(glob.n)
  assert size >= 1
    
  if not glob.silent:  
    print('cfs_get_drho_vector: returns zero vector of size ', size)

  return np.zeros(size)

## get the gradient of a density based function w.r.t. all shape variables.
# we calculate it for all and don't skip fixed variables.
# The element wise drho_vec * d_shape can be queried via cfs_get_field_info
# @param drho is a 1D np.array with d func/d rho - this is to be handled via chainrule.
# @param label the name of the current drho_vec's function to allow for cfs_get_field_info  
# @return d func/ d shape  
def cfs_get_gradient(ddes_vec, label):
  if glob.grad_rho is None:
    cfs_map_to_design()

  var_total = glob.total()
  sens = np.zeros(var_total)   # return vector 
  # drho_vec represents a 2D matrix, glob.grad_field is three dimensional with the last dimension the variable 
  for v in range(var_total):
    ds = glob.grad_rho[:,:,:,v].reshape(np.prod(glob.n), order='F') # the shape sensitivity as vector
    # sens_field = element wise the sensitivity of the function times the shape variables.
    sens_field = ddes_vec*ds
    sens[v] = sens_field.sum() # equals for i in range(nx): for j in range(ny): for k in range(nz): sens +=  drho[i,j,k] * glob.grad_field[i,j,k][v]
  #sens = np.sum(np.expand_dims(ddes_vec,1)*glob.grad_rho[:,:,:,v].reshape((np.prod(glob.n),var_total), order='F'),0)
  return sens 


class Spaghetti:
  # @param id starting from 0
  # @param base sum of all len(optvar()) of all previous shapes. Starts with 0
  def __init__(self,id,base,Px,Py,Pz,r,p,alpha):
    self.id = id
    self.base = base
    self.set(Px,Py,Pz,r,p,alpha)

  # does not only set the variables but also computes the internal helpers (C, ...) and deletes glob index fields
  def set(self, Px,Py,Pz,r,p,alpha, reset_fields=True):

    self.profile = p
    self.alpha = alpha
    self.P = np.array([Px,Py,Pz]).transpose()
    self.radii = np.concatenate([np.zeros(1),np.array(r),np.zeros(1)])
    self.num_optvar = self.P.size + self.radii.size # -2+1+1 because of non-design radii extension, p, alpha

    # helper
    self.m = len(Px) - 1 # number of segments
    self.n0_arc = np.zeros((self.m, 3))
    self.projection = np.zeros((self.m,3))
    self.fe_list_full = []
    self.fe_list_intermediate = []

    if reset_fields:
      glob.rho = None
      glob.grad_rho = None

    self.sens = {
            'grad_P': np.zeros((len(r)+2, self.num_optvar, 3)),
            'grad_t': np.zeros((len(r)+1, self.num_optvar, 3)),
            'grad_gamma': np.zeros((len(r)+2, self.num_optvar)), # extended with empty first gradient to be consistent with paper notation
            'grad_C': np.zeros((len(r)+2, self.num_optvar, 3)),
            'grad_n0_arc': np.zeros((len(r)+1, self.num_optvar, 3))} # extended with empty first gradient to be consistent with indices of arc midpoints C

    for i in range(self.m+1):
      self.sens['grad_P'][i][3*i:3*(i+1),:] = [[1,0,0], [0,1,0], [0,0,1]]

    # list of segment directions as normalized vectors T
    # Careful: segments are 0-based here and 1-based in paper!
    self.t = []
    for i in range(self.m): # segment index from 0 to m-1
      svec = self.P[i+1]-self.P[i]
      nsvec = norm(svec)
      self.t.append(svec/nsvec)
      nsvec3inv = 1/(nsvec**3)
      #self.sens['grad_t'][i][3*i:3*(i+1),:] = -self.sens['grad_P'][i][3*i:3*(i+1),:]/nsvec - np.multiply(svec,np.transpose(np.broadcast_to(svec,(3,3))))*nsvec3inv
      #self.sens['grad_t'][i][3*(i+1):3*(i+2),:] = self.sens['grad_P'][i+1][3*(i+1):3*(i+2),:]/nsvec - np.multiply(svec,np.transpose(np.broadcast_to(svec,(3,3))))*nsvec3inv
      self.sens['grad_t'][i][3*i,0] = -(svec[1]**2+svec[2]**2)*nsvec3inv
      self.sens['grad_t'][i][3*i,1] = svec[0]*svec[1]*nsvec3inv
      self.sens['grad_t'][i][3*i,2] = svec[0]*svec[2]*nsvec3inv
      self.sens['grad_t'][i][3*(i+1),:] = -self.sens['grad_t'][i][3*i,:]
      self.sens['grad_t'][i][3*i+1,0] = svec[1]*svec[0]*nsvec3inv
      self.sens['grad_t'][i][3*i+1,1] = -(svec[0]**2+svec[2]**2)*nsvec3inv
      self.sens['grad_t'][i][3*i+1,2] = svec[1]*svec[2]*nsvec3inv
      self.sens['grad_t'][i][3*(i+1)+1,:] = -self.sens['grad_t'][i][3*i+1,:]
      self.sens['grad_t'][i][3*i+2,0] = svec[2]*svec[0]*nsvec3inv
      self.sens['grad_t'][i][3*i+2,1] = svec[2]*svec[1]*nsvec3inv
      self.sens['grad_t'][i][3*i+2,2] = -(svec[0]**2+svec[1]**2)*nsvec3inv
      self.sens['grad_t'][i][3*(i+1)+2,:] = -self.sens['grad_t'][i][3*i+2,:]

    # For the arcs the center points C^i
    # the idea is to find the angle gamma between two segments, have vector B at half the angle
    # now scale radius/sin(gamma/2) along B0 to find C
    self.C = []
    self.C.append(self.P[0]) # add P^0 as C^0
    # note there is also list straight segment parts L from P->K_1s, K_1e->K_2s, ...->Q. For one seg this is P->Q
    # L is exclusively used in visualization
    self.L = []
    self.L.append([self.P[0]]) # we will append K1
    # list of arc is only meant as a helper for plotting the arcs
    self.arc = []
    self.gamma = [0] # artificial gamma_0, should not be used

    self.sens['grad_C'][0][0:3,:] = [[1,0,0], [0,1,0], [0,0,1]] # add derivatives for C^0 = P^0 and C^m = P^m 
    self.sens['grad_C'][self.m][self.m*3:(self.m+1)*3,:] = [[1,0,0], [0,1,0], [0,0,1]]
    for i in range(1,self.m):
      r = self.radii[i]
      Pi = self.P[i]
      tim1 = self.t[i-1]     # normalized segment before P (0-based in code!!)
      ti = self.t[i]    # normalized segment after P

      # get angle between two segments, standard form for angle between two vectors
      # cos(gamma) = -ti*tiplus1/ (||ti|| * ||tiplus1||)
      # If a=0 -> gamma = 0 and numerics goes crazy for the AD case -> see solution there
      cosa = np.dot(-tim1,ti)
      cosa = np.clip(cosa,-1,1)
      gamma = np.arccos(cosa)/2
      self.sens['grad_C'][i] = self.sens['grad_P'][i]
      self.gamma.append(gamma)
      B = ti - tim1
      nB2 = np.dot(B,B)
      nB = np.sqrt(nB2)
      if nB < 1e-14:
        # set grad_gamma only in else case to prevent division by 0
        # grad_gamma only appears in the derivatives of arc segments, which vanish for ti = tiplus1
        # should there still be an integration point in the numerically small remaining area,
        # the spaghetti will be straight and dC^i/da^i = dP^i/da^i without the curve should be good enough as approximation
        scaling = 0 # -> Ci = Pi
        B0 = B
      else:
        self.sens['grad_gamma'][i] = (np.dot(self.sens['grad_t'][i-1][:],ti) + np.dot(self.sens['grad_t'][i][:],tim1))/(2*np.sqrt(1-np.dot(tim1,ti)**2))
        scaling = r/np.sin(gamma) # would be r/sin(0)
        B0 = B/nB
        #self.sens['grad_C'][i] += r*(self.sens['grad_t'][i][:]-self.sens['grad_t'][i-1][:])/(np.sin(gamma)*nB) - r*(np.expand_dims(np.dot(self.sens['grad_t'][i],ti) - np.dot(self.sens['grad_t'][i-1][:],ti),axis=1))*B/(np.sin(gamma)*nB2**(1.5)) - r*np.cos(gamma)*B*np.expand_dims(self.sens['grad_gamma'][i],axis=1)/(np.sin(gamma)**2*nB)
        #self.sens['grad_C'][i] += r/(np.sin(gamma))*(self.sens['grad_t'][i][:]-self.sens['grad_t'][i-1][:])/nB - r/(np.sin(gamma))*B/(2*nB**3)*2*(self.sens['grad_t'][i][:]-self.sens['grad_t'][i-1][:]) - r*np.cos(gamma)*B*np.expand_dims(self.sens['grad_gamma'][i],axis=1)/(np.sin(gamma)**2*nB)
        for idx in range(self.P.size):
          self.sens['grad_C'][i][idx,:] = self.sens['grad_C'][i][idx,:] + self.radii[i]*(
            (self.sens['grad_t'][i][idx,:]-self.sens['grad_t'][i-1][idx,:])/(np.sin(gamma)*norm(ti-tim1))
            -np.cos(gamma)/(np.sin(gamma)**2)*self.sens['grad_gamma'][i][idx]*(ti-tim1)/norm(ti-tim1)
            -(ti-tim1)/(np.sin(gamma)*norm(ti-tim1)**3)*np.dot(ti-tim1,self.sens['grad_t'][i][idx,:]-self.sens['grad_t'][i-1][idx,:]))
        self.sens['grad_C'][i][self.P.size+i-1,:] = (ti-tim1)/(norm(ti-tim1)*np.sin(gamma))

      # from function of hypotenuse of a right-angled triangle
      C = Pi + scaling * B0
     # for namedvar in self.optvar_names():
     #   idx = self.namedidx[namedvar]
     #   self.sens['grad_C'][i][idx,:] = self.sens['grad_H'][i][idx,:]-

      # projection onto segment
      K1 = Pi-tim1*np.dot(-tim1,C-Pi)
      Pim1 = self.P[i-1]
      #K1 = Pi + np.dot(C-Pi,Pim1-Pi)/np.dot(Pim1-Pi,Pim1-Pi)*(Pim1-Pi)
      K2 = Pi+ti*np.dot(ti,C-Pi)
      Pim1 = self.P[i+1]
      #K2 = Pi + np.dot(C-Pi,Pim1-Pi)/np.dot(Pim1-Pi,Pim1-Pi)*(Pim1-Pi)

      self.C.append(C)
      self.L[-1].append(K1) # finish last (initial) L
      self.L.append([K2]) # to have added K1 for next arc in next iteration or P^m after the loop

      # only as helper for plotting the stuff.
      self.arc.append((C,gamma, K1, K2, r))

    self.C.append(self.P[-1]) # append P^m as C^m
    self.L[-1].append(self.P[-1])

    assert len(self.C)-2 == len(self.L)-1 # one segment more than arcs
    assert len(self.L) == len(self.t)   # L is part of the T vector

    for i in range(1,self.m):
      Pi = self.P[i]
      normal_arc = np.cross(-self.t[i-1],self.t[i])
      self.n0_arc[i] = normal_arc/norm(normal_arc)
      grad_normal_arc = np.cross(-self.sens['grad_t'][i-1],self.t[i])+np.cross(-self.t[i-1],self.sens['grad_t'][i])
      self.sens['grad_n0_arc'][i] = grad_normal_arc/norm(normal_arc) - normal_arc/norm(normal_arc)**3*np.dot(grad_normal_arc, np.expand_dims(normal_arc, axis=1))

  # give optimzation variables as defined array such that we can easily differentiate
  # @return profile,P^0_x,P^0_y,P^0_z,P^1_x,P^1_y,P^1_z,...,P^m_x,P^m_y,P^m_z,r_1,...r_m-1
  def optvar(self):
    r = np.zeros(self.num_optvar)
    #print(np.ravel(self.P, order='C'))
    r[0:self.P.size] = np.ravel(self.P, order='C')
    r[self.P.size:-2] = self.radii[1:-1]
    r[-2] = self.profile
    r[-1] = self.alpha

    return r
  
  # return spaghetti parameters from ordered array of variables as above
  def varopt(self, var):
    P = np.reshape(var[0:self.P.size], self.P.shape, order='C')
    radii = var[self.P.size:-2]
    profile = var[-2]
    alpha = var[-1]
    return P, radii, profile, alpha

  # search closest distance to all segments and arcs from point X
  def dist(self, X, where = None, what = 'all'):

    profile = self.profile # profile, we are negative inside and 0 at the boundary
    m = self.m # number of segments, one arc less

    minimal = (None,-1)

    # first m segments 0...m-1
    if not where or where == 's':
      for i in range(len(self.t)):
        C0 = self.C[i] # point on the line perpendicular to the start of the straight segment i (P^0 for first segment)
        C1 = self.C[i+1] # point on the line perpendicular to the end of the straight segment i (P^m for last segment)
        t = self.t[i]
        Pi = self.P[i]
        # we are in range if (X-C^i+1) @ T^i <= 0 and (X-C^i) @ T^i >= 0.
        if (X-C1) @ t <= 0 and (X-C0) @ t >= 0:
          # (X-g) @ M is distance to segment. Positive on M side, negative on the other side. Instead of g, we can use H
          minimal = idx_min(minimal, (norm(np.cross((X-Pi),t)) - profile,i)) # don't forget abs!

    # arcs
    if not where or where == 'c':
      # n-1 arcs n ... 2n-2
      for i in range(1, m): # arcs only are around interior C
        C = self.C[i]
        r = self.radii[i]
        # check if we are in the cone (base is the center)
        v1 = self.t[i-1]
        v2 = -self.t[i]
        XC = X - C

        if v1 @ XC > 0 and v2 @ XC > 0: # we are in the cone

          #assert len(self.E) == len(self.C)-2
          #E = self.E[i]
          normal = np.cross(v1,v2)
          normal = normal/norm(normal)
          X_projected = X - np.dot(X-self.P[i],normal)*normal # project to arc plane
          # pythagoras with distance to arc and distance of projection
          d = np.sqrt((r-norm(C-X_projected))**2+norm(X-X_projected)**2) - profile

          #print('dist cone', X,C,XC,norm(XC),abs(norm(XC) - r),d)
          minimal = idx_min(minimal, (d, m+i-1))

      assert minimal[1] <= 2*m-2

      # distance to start point
      minimal = idx_min(minimal,(norm(X-self.P[0])-profile,2*m-1))
      # distance to end point
      minimal = idx_min(minimal, (norm(X-self.P[m])-profile,2*m)) # faster would be root after min()

    if what == 'index':
      return minimal[1]
    elif what == 'distance':
      return (minimal[0], np.zeros(self.num_optvar)) # add fake gradient to be compatible with fast_dist
    elif what == 'gradient' or what == 'didxgrad':
      idx = minimal[1]
      grad = np.zeros(self.num_optvar)
      if idx == 2*m-1:
        grad[0:3] = (self.P[0]-X)/norm(X-self.P[0]) if norm(X-self.P[0]) > 1e-15 else np.ones(3)
      elif idx == 2*m:
        grad[self.m*3:(self.m+1)*3] = (self.P[m]-X)/norm(X-self.P[m]) if norm(X-self.P[m]) > 1e-15 else np.ones(3)
      elif idx < 0:
        grad = np.zeros(self.num_optvar)
      elif idx < m:
        Pi = self.P[idx]
        t = self.t[idx]
        grad = np.matmul(np.cross((X-Pi),self.sens['grad_t'][idx]) - np.cross(self.sens['grad_P'][idx],t),np.cross((X-Pi),t))/norm(np.cross((X-Pi),t))
      elif idx < 2*m-1:
        i = idx - m + 1
        C = self.C[i]
        r = self.radii[i]
        Pi = self.P[i]
        n0_arc = self.n0_arc[i]
        X_projected = X - np.dot(X-Pi,n0_arc)*n0_arc # project to arc plane
        for varidx in range(self.num_optvar-1):
          grad_X_projected = np.matmul(self.sens['grad_P'][i][varidx,:], n0_arc)*n0_arc + np.matmul(Pi-X, self.sens['grad_n0_arc'][i][varidx,:])*n0_arc + np.matmul(Pi-X, n0_arc)*self.sens['grad_n0_arc'][i][varidx,:]
          grad_r = 1 if varidx == self.P.size+i-1 else 0
          grad[varidx] = 0.5/np.sqrt((r-norm(C-X_projected))**2+norm(X-X_projected)**2)*(2*(r-norm(C-X_projected))*(grad_r - 1/norm(C-X_projected)*(np.matmul(C-X_projected,self.sens['grad_C'][i][varidx,:] - grad_X_projected))) - 2*np.matmul(X-X_projected, grad_X_projected))

      grad[-2] = -1
      if what == 'didxgrad':
        return (minimal, grad)
      return (minimal[0], grad)
    elif what == 'all':
      return minimal
    else:
      assert(False)


  # search closest distance by part index
  def dist_by_idx(self, X, indexlist, what='gradient'):

    profile = self.profile # profile, we are negative inside and 0 at the boundary
    m = self.m # number of segments, one arc less

    minimal = (None,-3)

    for idx in indexlist:

      if idx == 2*m-1: # distance to start point
        minimal = idx_min(minimal,(norm(X-self.P[0])-profile,idx))
      elif idx == 2*m: # distance to end point
        minimal = idx_min(minimal, (norm(X-self.P[m])-profile,idx)) # faster would be root after min()

      if idx < m:
        C0 = self.C[idx] # point on the line perpendicular to the start of the straight segment i (P^0 for first segment)
        C1 = self.C[idx+1] # point on the line perpendicular to the end of the straight segment i (P^m for last segment)
        t = self.t[idx]
        Pi = self.P[idx]
        # we are in range if (X-C^i+1) @ T^i <= 0 and (X-C^i) @ T^i >= 0.
        if ((X-C1) @ t <= 0 and (X-C0) @ t >= 0):
          # (X-g) @ M is distance to segment. Positive on M side, negative on the other side. Instead of g, we can use H
          minimal = idx_min(minimal, (norm(np.cross((X-Pi),t)) - profile,idx)) # don't forget abs!
        #else:
         # minimal = idx_min(minimal, (norm(np.cross((X-Pi),t)) - profile,idx)) # don't forget abs!
          #print('outside of segment', minimal, indexlist)
      elif idx < 2*m-1:
        i = idx - m + 1
        C = self.C[i]
        r = self.radii[i]
        # check if we are in the cone (base is the center)
        v1 = self.t[i-1]
        v2 = -self.t[i]
        XC = X - C

        if (v1 @ XC > 0 and v2 @ XC > 0): # we are in the cone

          #assert len(self.E) == len(self.C)-2
          #E = self.E[i]
          normal = np.cross(v1,v2)
          normal = normal/norm(normal)
          X_projected = X - np.dot(X-self.P[i],normal)*normal # project to arc plane
          # pythagoras with distance to arc and distance of projection
          d = np.sqrt((r-norm(C-X_projected))**2+norm(X-X_projected)**2) - profile

          #print('dist cone', X,C,XC,norm(XC),abs(norm(XC) - r),d)
          minimal = idx_min(minimal, (d, idx))

    minimal2 = self.dist(X)
    if minimal2[1] != minimal[1] and minimal2[0] >= -glob.h and minimal2[0] < glob.h:
      print('dist and dist_by_idx do not match! id=',self.id, 'dist=', minimal2[0], self.nice_idx(minimal2[1]), 'dist_by_idx=', minimal[0], self.nice_idx(minimal[1]))

    if minimal[0] is None:
      d = self.dist(X)
      if d[0] >= -glob.h and d[0] <= glob.h:
        print('Distance is None at point', X, 'with indices', indexlist, '. Actual distance is', d)
      return (glob.transition, np.zeros(self.num_optvar))

    if what == 'gradient' or what == 'didxgrad':
      idx = minimal[1]
      grad = np.zeros(self.num_optvar)
      if idx == 2*m-1: # starting point
        grad[0:3] = (self.P[0]-X)/norm(X-self.P[0]) if norm(X-self.P[0]) > 1e-15 else np.ones(3)
      elif idx == 2*m: # end point
        grad[self.m*3:(self.m+1)*3] = (self.P[m]-X)/norm(X-self.P[m]) if norm(X-self.P[m]) > 1e-15 else np.ones(3)
      elif idx < 0:
        grad = np.zeros(self.num_optvar)
      elif idx < m: # segment idx
        Pi = self.P[idx]
        t = self.t[idx]
        grad = np.matmul(np.cross((X-Pi),self.sens['grad_t'][idx]) - np.cross(self.sens['grad_P'][idx],t),np.cross((X-Pi),t))/norm(np.cross((X-Pi),t))
      elif idx < 2*m-1: # arc idx-m+1
        i = idx - m + 1
        C = self.C[i]
        r = self.radii[i]
        Pi = self.P[i]
        n0_arc = self.n0_arc[i]
        X_projected = X - np.dot(X-Pi,n0_arc)*n0_arc # project to arc plane
        for varidx in range(self.num_optvar-1):
          grad_X_projected = np.matmul(self.sens['grad_P'][i][varidx,:], n0_arc)*n0_arc + np.matmul(Pi-X, self.sens['grad_n0_arc'][i][varidx,:])*n0_arc + np.matmul(Pi-X, n0_arc)*self.sens['grad_n0_arc'][i][varidx,:]
          grad_r = 1 if varidx == self.P.size+i-1 else 0
          grad[varidx] = 0.5/np.sqrt((r-norm(C-X_projected))**2+norm(X-X_projected)**2)*(2*(r-norm(C-X_projected))*(grad_r - 1/norm(C-X_projected)*(np.matmul(C-X_projected,self.sens['grad_C'][i][varidx,:] - grad_X_projected))) - 2*np.matmul(X-X_projected, grad_X_projected))

      grad[-2] = -1
      if what == 'didxgrad':
        return (minimal, grad)
      return (minimal[0], grad)
    elif what == 'all':
      return minimal
    elif what == 'index':
      return minimal[1]
    elif what == 'distance':
      return (minimal[0], np.zeros(self.num_optvar)) # add fake gradient to be compatible with fast_dist
    else:
      assert(False)


  # return indices of potential segments cutting element and whether element is full material
  def fe_get_indices(self, X, dx, tr_half):
    idx_intermediate = [self.id]

    lax = 0.51*dx if glob.order == 2 else 0.01*dx # relax index search to distance to integration points
    profile = self.profile # profile, we are negative inside and 0 at the boundary
    m = self.m # number of segments, one arc less

    # first m segments 0...m-1
    for i in range(len(self.t)):
      Cip = self.C[i] # point on the line perpendicular to the start of the straight segment i (P^0 for first segment)
      Cin = self.C[i+1] # point on the line perpendicular to the end of the straight segment i (P^m for last segment)
      ti = self.t[i]
      Pi = self.P[i]
      # we are in range if (X-C^i+1) @ T^i <= 0 and (X-C^i) @ T^i >= 0. Relax to include whole element
      if ((X-lax*ti)-Cin) @ ti <= 0 and ((X+lax*ti)-Cip) @ ti >= 0:
        # (X-g) @ M is distance to segment. Positive on M side, negative on the other side. Instead of g, we can use H
        dist = norm(np.cross((X-Pi),ti)) - profile
        if dist < -tr_half - lax:
          return [], True
        elif dist < tr_half + lax:
          idx_intermediate.append(i)


    # arcs
    # n-1 arcs n ... 2n-2
    for i in range(1, m): # arcs only are around interior C
      C = self.C[i]
      r = self.radii[i]
      # check if we are in the cone (base is the center)
      tip = self.t[i-1]
      tin = -self.t[i]

      if tip @ ((X+lax*tip)-C) > 0 and tin @ ((X+lax*tin)-C) > 0: # we are in the cone or just outside

        #assert len(self.E) == len(self.C)-2
        #E = self.E[i]
        normal = np.cross(tip,tin)
        normal = normal/norm(normal)
        X_projected = X - np.dot(X-self.P[i],normal)*normal # project to arc plane
        # pythagoras with distance to arc and distance of projection
        d = np.sqrt((r-norm(C-X_projected))**2+norm(X-X_projected)**2) - profile

        if d < - tr_half - lax:
          return [], True
        elif d < tr_half + lax:
          idx_intermediate.append(m+i-1)

    # distance to start point
    d = norm(X-self.P[0])-profile
    if d < - tr_half - lax:
      return [], True
    elif d < tr_half + lax:
      idx_intermediate.append(2*m-1)
    # distance to end point
    d = norm(X-self.P[m])-profile
    if d < - tr_half - lax:
      return [], True
    elif d < tr_half + lax:
      idx_intermediate.append(2*m)

    if len(idx_intermediate) == 1:
      idx_intermediate = []

    return idx_intermediate, False


  # return bounding box of spaghetti_shape
  def get_bounding_box(self):
    bbox_min = np.min(self.P, axis=0) - self.profile - 1.1*glob.h
    bbox_max = np.max(self.P, axis=0) + self.profile + 1.1*glob.h
    return bbox_min, bbox_max

  # setup FE list with necessary distance computations for this shape
  def setup_fe_list(self):
    self.fe_list_intermediate = []
    self.fe_list_full = []
    if self.alpha < 1e-10:
      return
    bbox_min, bbox_max = self.get_bounding_box()
    fe_min = np.maximum(np.array(np.floor((bbox_min-glob.box_lower)/glob.dx)-1,dtype=int),0)
    fe_max = np.minimum(np.array(np.ceil((bbox_max-glob.box_lower)/glob.dx+1),dtype=int),np.array(glob.n,dtype=int))
    if not glob.silent:
      print(fe_min, fe_max)
    for i in range(fe_min[0], fe_max[0]):
      for j in range(fe_min[1], fe_max[1]):
        for k in range(fe_min[2], fe_max[2]):
          Xmidijk = glob.mid_lower + glob.dx*np.array((i,j,k))
          shape_idx_list, full = self.fe_get_indices(Xmidijk, glob.dx, glob.h)
          if full:
            self.fe_list_full.append([i,j,k])
          elif shape_idx_list:
           self.fe_list_intermediate.append([i,j,k,shape_idx_list])
    if not glob.silent:
      print('shape id:', self.id, 'number of full elements:', len(self.fe_list_full), 'number of intermediate elements:', len(self.fe_list_intermediate))
    #return [self.fe_list_full, self.fe_list_intermediate]

  # get constraint value to prevent arcs moving too close together and overlapping
  def get_constraint_arc_overlap(self, segnum):
    assert(len(self.C)>segnum)
    return np.dot(self.P[segnum]-self.P[segnum-1], self.C[segnum]-self.C[segnum-1])

  # get gradient value
  def get_gradient_arc_overlap(self, segnum):
    assert(len(self.radii)>segnum)
    return np.dot(self.sens['grad_P'][segnum]-self.sens['grad_P'][segnum-1], self.C[segnum]-self.C[segnum-1]) + \
      np.dot(self.sens['grad_C'][segnum]-self.sens['grad_C'][segnum-1], self.P[segnum]-self.P[segnum-1])


  # for debug purpose, the part idx human reable printed
  # the parts are 1-based!
  def nice_idx(self, idx):
    m = self.m
    i = idx
    assert i >= -2 and i <= m + m-1 + 2 # m segments, m-1 arcs, P, Q
    if i == -2:
      return "outside"
    if i == -1:
      return "inside"
    if i < m:
      return "seg_" + str(i)
    if i < m + m-1:
      return "arc_" + str(i - m)
    if i == m + m-1:
      return "P0"
    if i == m + m-1 + 1:
      return "Pm"
    if i == -3:
      return "not set"
    assert False

    # for gradient check only
  def func(self, var,args):
    #print(args)
    P, radii, profile, alpha = self.varopt(var)
    #print("varopt P: ", P)
    #print(P[0::3][0])
    self.set(P[:,0],P[:,1],P[:,2],radii,profile,alpha, reset_fields=False)
    if str(args[-1]).startswith('P'):
      #print(args)
      #print(self.P)
      return self.P[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('t'):
      return self.t[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('gamma'):
      return self.gamma[int(args[-1][-1])]
    elif str(args[-1]).startswith('C'):
      return self.C[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('n0_arc'):
      return self.n0_arc[int(args[-1][-1])][args[0]]
    elif str(args[-1]) == 'distance':
      dist = self.dist(args[0], what='distance')
    elif str(args[-1]) == 'boundary':
      dist = self.alpha*boundary(self.dist(args[0], what='distance'))
      return dist
    #res = self.fast_dist(args[0], args[1])
    #return res[0]
    return dist[0]

  # for gradient check only
  def grad(self, var,args):
    P, radii, profile, alpha = self.varopt(var)
    self.set(P[:,0],P[:,1],P[:,2],radii,profile,alpha, reset_fields=False)
    if str(args[-1]).startswith('P'):
      return self.sens['grad_P'][int(args[-1][-1])][:,args[0]]
    elif str(args[-1]).startswith('t'):
      return self.sens['grad_t'][int(args[-1][-1])][:,args[0]]
    elif str(args[-1]).startswith('gamma'):
      return self.sens['grad_gamma'][int(args[-1][-1])]
    elif str(args[-1]).startswith('C'):
      return self.sens['grad_C'][int(args[-1][-1])][:,args[0]]
    elif str(args[-1]).startswith('n0_arc'):
      return self.sens['grad_n0_arc'][int(args[-1][-1])][:,args[0]]
    elif str(args[-1]) == 'distance':
      dist = self.dist(args[0], what='gradient')
    elif str(args[-1]) == 'boundary':
      dist = boundary(self.dist(args[0], what='gradient'))
    #res = self.fast_dist(args[0], args[1])
    #print(res)
    return dist[1]

def setup_fe_list_glob(shape):
  shape.setup_fe_list()
  return [shape.fe_list_intermediate, shape.fe_list_full]

# returns the boundary modelling.
# @param dist positive inside, negative outside, boundary at 0
# @param transition is 2*h in feature mapping paper
# @return if outside |transition| -> 0 or 1
def boundary(dist, derivative=False, alpha=None):
  rm = glob.rhomin
  rmx = glob.rhomax
  # in the non-derivative case the first (or in the vtk case, only) value used. Whatever dist[1] is otherwise?!
  phi = -dist if not derivative and type(dist) is float or type(dist) is np.float64 else  -dist[0] # positive inside, negative outside
  h = glob.h
  #phi = phi-h
  if phi <= -h:
    rho = rm
  elif phi >= h:
    rho = rmx
  elif glob.boundary == 'linear':
    rho = .5*((rmx-rm)*phi/h+rmx+rm)
  elif glob.boundary == 'poly':
    rho = 3.0/4.0*(rmx - rm) * (phi/h - phi**3/(3*h**3)) + .5 * (rmx+rm)
  else:
    print("Error: boundary type '" + glob.boundary + "' not implemented!")
    os.sys.exit()
  grad = dist[1]
  if alpha is not None:
    rho_orig = rho
    rho = alpha*rho
    grad = alpha*dist[1]
  if derivative != True:
    return rho
  else:
    if phi <= -h:
      return (rho, np.zeros(dist[1].shape))
    elif phi >= h:
      if not alpha:
        return (rho, np.zeros(dist[1].shape))
      else:
        grad = np.zeros(dist[1].shape)
        grad[-1] = rho_orig
        return (rho, grad)
    elif glob.boundary == 'linear':
      grad *= .5*((rmx-rm)/h)
      if alpha:
        grad[-1] = rho_orig
      return (rho, grad)
    elif glob.boundary == 'poly':
      grad *= -3.0/4.0*(rmx - rm)*(1/h - phi**2/(h**3))
      if alpha:
        grad[-1] = rho_orig
      return (rho, grad)


# create a idx field for fast access where >= 0 for closest part, -1 for inside and -2 for too far outside
# uses glob
# @param discretization is by default glob.n[0]+1
# @return fields of size (n+1)^dim, to capture the nodes of the density elements.
#         The first field contains index vectors (len(shapes)), the second a vector of rhos (len(shapes))
def create_idx_field(discretization = None):

  shapes = glob.shapes
  Nx = discretization if discretization else glob.n[0]+1
  Ny = discretization if discretization else glob.n[1]+1
  Nz = discretization if discretization else glob.n[2]+1

  h = glob.transition *.55 # is save to add a little buffer

  x_ls = np.linspace(0,glob.dx*glob.n[0],Nx)
  y_ls = np.linspace(0,glob.dx*glob.n[1],Ny)
  z_ls = np.linspace(0,glob.dx*glob.n[2],Nz)

  idx = np.ones((Nx,Ny,Nz,len(shapes)),dtype=int)
  idx_shapes_only = np.ones((Nx,Ny,Nz,len(shapes)),dtype=int)
  grad_field = np.zeros((Nx,Ny,Nz,len(shapes),glob.total()))
  dist  = np.ones((Nx,Ny,Nz,len(shapes)))

  for i, x in enumerate(x_ls):
    for j, y in enumerate(y_ls):
      for k, z in enumerate(z_ls):
        X = [x,y,z]
        for si, s in enumerate(shapes):
          mini, grad = s.dist(X, what='didxgrad')
          d = mini[0]
          dist[i,j,k,si] = d
          grad_field[i,j,k,si,s.base:(s.base+s.num_optvar)] = grad
          idx_shapes_only[i,j,k,si] = mini[1]
          idx[i,j,k,si] = mini[1] if d > -h and d < h else (-1 if d < h else -2)
          # print('cif: i,j,X,d,k',i,j,X,d,k)

  return idx, dist, idx_shapes_only

# setup fe lists in shapes
def create_fe_lists():
  with Pool(glob.num_threads) as pool:
    res = pool.map(setup_fe_list_glob, glob.shapes)
    pool.close()
    pool.join()
  for i, s in enumerate(glob.shapes):
    s.fe_list_intermediate = res[i][0]
    s.fe_list_full = res[i][1]
  #for s in glob.shapes:
   # s.setup_fe_list()

# setup fe lists in shapes
def create_fe_lists_serial():
  for s in glob.shapes:
    s.setup_fe_list()

# assemble fe field list from shapes
def assemble_fe_field(discretization = None):
  start = ti.time()
  create_fe_lists()
  if not glob.silent:
    print('setting up FE lists in parallel:', ti.time()-start, 'seconds')
  #start = ti.time()
  #create_fe_lists_serial()
  #print('setting up FE lists serially:', ti.time()-start, 'seconds')
  start = ti.time()
  glob.fe_list = [[[FE(i,j,k) for i in range(glob.n[2])] for j in range(glob.n[1])] for k in range(glob.n[0])]
  if not glob.silent:
    print(len(glob.fe_list),len(glob.fe_list[0]),len(glob.fe_list[0][0]))
  for snum, s in enumerate(glob.shapes):
    for idx in s.fe_list_full:
      glob.fe_list[idx[0]][idx[1]][idx[2]].shapes_full.append(snum)
    for idx in s.fe_list_intermediate:
      glob.fe_list[idx[0]][idx[1]][idx[2]].shapes_idx_list.append(idx[3])
  if not glob.silent:
    print('assembling global FE field:', ti.time()-start, 'seconds')


# integrate fe field
def integrate_fe_field():
  order = glob.order
  #glob.dist_field = np.ones(glob.n)
  glob.rho = np.zeros(glob.n)
  glob.grad_rho = np.zeros((glob.n[0],glob.n[1],glob.n[2],glob.total()))
  dx = glob.dx
  eta = .5*dx/np.sqrt(3) # distance for Gauss-points
  total_fe = 0
  for i in range(glob.n[0]):
    for j in range(glob.n[1]):
      for k in range(glob.n[2]):
        rho_sm = 0
        grad_rho = np.zeros(glob.total())
        w_s_full = []
        rho_s_full = []
        if glob.fe_list[i][j][k].shapes_full: # at least one shape full
          #glob.dist_field[i,j,k] = -glob.transition
          for snum in glob.fe_list[i][j][k].shapes_full:
            rho_s_full.append(glob.shapes[snum].alpha)
            w_s_full.append(np.exp(glob.p*glob.shapes[snum].alpha))
        rho_s_full = np.array(rho_s_full)
        w_s_full = np.array(w_s_full)
        sum_w_s_full = np.sum(w_s_full)
        num_full = len(glob.fe_list[i][j][k].shapes_full)
        if glob.fe_list[i][j][k].shapes_idx_list: # list not empty, have to integrate
          total_fe += 1
          Xmid = glob.mid_lower + np.array([i,j,k],dtype=int)*dx
          XX = [[Xmid[0]+it[0], Xmid[1]+it[1], Xmid[2]+it[2]] for it in product([-eta, eta],repeat=3)] if order == 2 else [Xmid]
          rho_s = []
          grad_rho_s = []
          s_idx_intermediate = []
          #d_min = 100*glob.transition
          for shape_idx_list in glob.fe_list[i][j][k].shapes_idx_list:
            s_idx_intermediate.append(shape_idx_list[0])
            #d, gradd = glob.shapes[shape_idx_list[0]].dist_by_idx(Xmid, shape_idx_list[1:])
            #if d < d_min:
            #  glob.dist_field[i,j,k] = d
            #  d_min = d
            #rho = boundary([d, gradd], derivative=True)
            rho, grad_ip = zip(*[boundary(glob.shapes[shape_idx_list[0]].dist_by_idx(X, shape_idx_list[1:]),derivative=True,alpha=glob.shapes[shape_idx_list[0]].alpha) for X in XX])
            rho_s.append(rho)
            grad_rho_s.append(grad_ip)
          rho = np.array(rho_s)
          w_s = np.exp(glob.p*rho)
          sum_w_s = np.sum(w_s,axis=0)+sum_w_s_full+len(glob.shapes)-len(w_s)-num_full
          rho_sm = (np.sum(w_s*rho,axis=0)+np.sum(rho_s_full*w_s_full))/sum_w_s
          rho_sm = np.sum(rho_sm)/order**3
          for num, sidx in enumerate(s_idx_intermediate):
            grad_rho[glob.shapes[sidx].base:glob.shapes[sidx].base+glob.shapes[sidx].num_optvar] = np.sum(np.expand_dims(w_s[num]/sum_w_s*(1+glob.p*(rho_s[num] - rho_sm)),axis=1)*grad_rho_s[num],axis=0)/order**3
          for num, sidx in enumerate(glob.fe_list[i][j][k].shapes_full):
            grad_rho[glob.shapes[sidx].base+glob.shapes[sidx].num_optvar-1] = np.sum(w_s_full[num]/sum_w_s)/order**3*(1+glob.p*(rho_s_full[num]-rho_sm))
          #rho_sm = boundary(d_min)
        else:
          rho_sm = np.sum(rho_s_full*w_s_full)/(sum_w_s_full+len(glob.shapes)-num_full) # doesn't work anymore with glob.alpha...
          for num, sidx in enumerate(glob.fe_list[i][j][k].shapes_full):
            grad_rho[glob.shapes[sidx].base+glob.shapes[sidx].num_optvar-1] = w_s_full[num]/(sum_w_s_full+len(glob.shapes)-num_full)*(1+glob.p*(rho_s_full[num]-rho_sm))
        glob.rho[i,j,k] = rho_sm
        glob.grad_rho[i,j,k] = grad_rho

  print(total_fe)

def integrate_fe(fe_num):
  i,j,k = np.unravel_index(fe_num, glob.n)

  order = glob.order
  dx = glob.dx
  eta = .5*dx/np.sqrt(3) # distance for Gauss-points

  rho_sm = 0
  grad_rho = np.zeros(glob.total())
  w_s_full = []
  rho_s_full = []
  if glob.fe_list[i][j][k].shapes_full: # at least one shape full
    #glob.dist_field[i,j,k] = -glob.transition
    for snum in glob.fe_list[i][j][k].shapes_full:
      rho_s_full.append(glob.rhomax*glob.shapes[snum].alpha)
      w_s_full.append(np.exp(glob.rhomax*glob.p*glob.shapes[snum].alpha))
  rho_s_full = np.array(rho_s_full)
  w_s_full = np.array(w_s_full)
  sum_w_s_full = np.sum(w_s_full)
  num_full = len(glob.fe_list[i][j][k].shapes_full)
  if glob.fe_list[i][j][k].shapes_idx_list: # list not empty, have to integrate
    Xmid = glob.mid_lower + np.array([i,j,k],dtype=int)*dx
    XX = [[Xmid[0]+it[0], Xmid[1]+it[1], Xmid[2]+it[2]] for it in product([-eta, eta],repeat=3)] if order == 2 else [Xmid]
    rho_s = []
    grad_rho_s = []
    s_idx_intermediate = []
    #d_min = 100*glob.transition
    for shape_idx_list in glob.fe_list[i][j][k].shapes_idx_list:
      s_idx_intermediate.append(shape_idx_list[0])
      #d, gradd = glob.shapes[shape_idx_list[0]].dist_by_idx(Xmid, shape_idx_list[1:])
      #if d < d_min:
      #  glob.dist_field[i,j,k] = d
      #  d_min = d
      #rho = boundary([d, gradd], derivative=True)
      rho, grad_ip = zip(*[boundary(glob.shapes[shape_idx_list[0]].dist_by_idx(X, shape_idx_list[1:]),derivative=True,alpha=glob.shapes[shape_idx_list[0]].alpha) for X in XX])
      rho_s.append(rho)
      grad_rho_s.append(grad_ip)
    rho = np.array(rho_s)
    w_s = np.exp(glob.p*rho)
    sum_w_s = np.sum(w_s,axis=0)+sum_w_s_full+len(glob.shapes)-len(w_s)-num_full
    rho_sm = (np.sum(w_s*rho,axis=0)+np.sum(rho_s_full*w_s_full))/sum_w_s
    rho_sm = np.sum(rho_sm)/order**3
    for num, sidx in enumerate(s_idx_intermediate):
      grad_rho[glob.shapes[sidx].base:glob.shapes[sidx].base+glob.shapes[sidx].num_optvar] = np.sum(np.expand_dims(w_s[num]/sum_w_s*(1+glob.p*(rho_s[num] - rho_sm)),axis=1)*grad_rho_s[num],axis=0)/order**3
    for num, sidx in enumerate(glob.fe_list[i][j][k].shapes_full):
      grad_rho[glob.shapes[sidx].base+glob.shapes[sidx].num_optvar-1] = np.sum(w_s_full[num]/sum_w_s)/order**3*(1+glob.p*(rho_s_full[num]-rho_sm))
    #rho_sm = boundary(d_min)
  else:
    rho_sm = np.sum(rho_s_full*w_s_full)/(sum_w_s_full+len(glob.shapes)-num_full)
    for num, sidx in enumerate(glob.fe_list[i][j][k].shapes_full):
      grad_rho[glob.shapes[sidx].base+glob.shapes[sidx].num_optvar-1] = w_s_full[num]/(sum_w_s_full+len(glob.shapes)-num_full)*(1+glob.p*(rho_s_full[num]-rho_sm))
  return (rho_sm, grad_rho)

# get constraint sparsity pattern
# did not implement exact sparsity pattern yet, just adding all variables of spaghetti
def cfs_get_sparsity_arc_overlap(opt):
  glob.jac = [] # list of arrays of 0-based variable indices
  cfs_jac = []
  for snum, s in enumerate(glob.shapes):
    if len(s.radii) == 2:
      continue
    base = s.base
    sparse = np.arange(base,base+s.num_optvar)
    for i in range(1,s.m+1):
      cfs_jac.append(np.array(sparse))
      glob.jac.append((snum, i, cfs_jac[-1]))
  return cfs_jac

# get constraint value to prevent arcs moving too close together and overlapping
def cfs_get_constraint_arc_overlap(constraint_num):
  shapenum, a_num, cfs_jac = glob.jac[constraint_num]
  shape = glob.shapes[shapenum]
  return shape.get_constraint_arc_overlap(a_num)

# get constraint gradient
def cfs_get_gradient_arc_overlap(constraint_num):
  shapenum, a_num, cfs_jac = glob.jac[constraint_num]
  shape = glob.shapes[shapenum]
  gradient_full = shape.get_gradient_arc_overlap(a_num)
  return gradient_full[cfs_jac-shape.base]


def get_vector_arc_overlap(var_all):
  const = []
  for snum, s in enumerate(glob.shapes):
    var = var_all[s.base:s.base+s.num_optvar]
    cfs_set_spaghetti(s.id, np.reshape(var[0:s.P.size],s.P.shape), var[s.P.size:-2], var[-2], var[-1])
    if len(s.radii) == 2:
      continue
    for i in range(1,s.m+1):
      const.append(s.get_constraint_arc_overlap(i))
  return np.array(const)

def get_jacobian_arc_overlap(var_opt):
  c_jac = np.zeros((len(get_vector_arc_overlap(var_opt)),len(var_opt)))
  cnum = 0
  for snum, s in enumerate(glob.shapes):
    if len(s.radii) == 2:
      continue
    for i in range(1,s.m+1):
      c_jac[cnum][s.base:s.base+s.num_optvar]=(s.get_gradient_arc_overlap(i))
      cnum += 1
  return c_jac


# do L2-tracking of density
def track(shapes, design_track, saveall = False, track_iter = 30):
  rho = cfs_map_to_design()
  dist = np.sum((design_track.reshape(np.prod(glob.n), order='F')-rho)**2)
  var_all = glob.var_all()
  glob.iter = 0
  print("l2-distance: ", dist)
  #err = sciopt.check_grad(eval_l2,grad_l2,var_all,design_track)
  if glob.gradient_check:
    print("gradient check of l2-tracking: ", sciopt.approx_fprime(var_all, eval_l2, 1.5e-8, [design_track,False])-grad_l2(var_all,[design_track,False]))
  bound = sciopt.Bounds(lb=np.zeros(len(var_all)),ub=2*np.ones(len(var_all)))
  #cfs_get_sparsity_arc_overlap(None)
  n_const = len(get_vector_arc_overlap(var_all))
  if n_const > 0:
    nonlinear_constraints = sciopt.NonlinearConstraint(get_vector_arc_overlap, np.zeros(n_const), np.inf*np.ones(n_const), jac=get_jacobian_arc_overlap)
    res = sciopt.minimize(eval_l2, var_all, args=[design_track, saveall], method="trust-constr", jac=grad_l2, constraints=nonlinear_constraints, bounds=bound, options={"disp": True, "verbose": 2, "maxiter": track_iter})
  else:
    res = sciopt.minimize(eval_l2, var_all, args=[design_track, saveall], method="trust-constr", jac=grad_l2, bounds=bound, options={"disp": True, "verbose": 2, "maxiter": track_iter})
  return res

# function value for L2 tracking optimization
def eval_l2(var_all, arg): # args = [density_track, save_figures]
  for shape in glob.shapes:
    shape_num = shape.id
    var = var_all[shape.base:shape.base+shape.num_optvar]
    cfs_set_spaghetti(shape_num, np.reshape(var[0:shape.P.size],shape.P.shape), var[shape.P.size:-2], var[-2], var[-1])
  des = cfs_map_to_design()
  if arg[1]:
    ot.write_density_file('giffiles/dens' + str(glob.iter).zfill(4) + '.density.xml', np.reshape(des, (glob.n[0],glob.n[1]), 'F'))
  return np.sum((des - arg[0].reshape(np.prod(glob.n), order='F'))**2)

# gradient for L2 tracking optimization
def grad_l2(var_all, arg): # args = [density_track, save_figures]
  for shape in glob.shapes:
    shape_num = shape.id
    var = var_all[shape.base:shape.base+shape.num_optvar]
    cfs_set_spaghetti(shape_num, np.reshape(var[0:shape.P.size],shape.P.shape), var[shape.P.size:-2], var[-2], var[-1])
  if arg[1]:
    try:
      fig = plot_data(1200,shapes,args.detail)
      fig.savefig('giffiles/track' + str(glob.iter).zfill(4) + '.png')
      plt.close(fig)
      glob.iter += 1
    except:
      print("Could not save file in iteration ", glob.iter)
    # read the xml file to memory
    xml = ut.open_xml(args.input)
    for num_var, var in enumerate(var_all):
      # modify the xml file
      ut.replace(xml, '//shapeParamElement[@nr='+str(num_var)+']/@design', str(var))

    # write to a new name
    xml.write('giffiles/track' + '.density.xml')
  drho = 2*(cfs_map_to_design() - arg[0].reshape(np.prod(glob.n), order='F'))
  sens = cfs_get_gradient(drho, "L2-tracking")
  return sens

# write distance values and that stuff
def write_vtk_lists(name):
  from pyevtk.hl import gridToVTK

  x_ls = np.linspace(glob.mid_lower[0],glob.box_upper[0]-0.5*glob.dx,int(glob.n[0]))
  y_ls = np.linspace(glob.mid_lower[1],glob.box_upper[1]-0.5*glob.dx,int(glob.n[1]))
  z_ls = np.linspace(glob.mid_lower[2],glob.box_upper[2]-0.5*glob.dx,int(glob.n[2]))

  pd={"rho": glob.rho}
  #pd["dist"] = glob.dist_field
  print('Writing results to', name+'.vtr')
  gridToVTK(name, x_ls,y_ls,z_ls, pointData=pd)

# write distance values and that stuff
def write_vtk(name,N):
  from pyevtk.hl import gridToVTK

  shapes = glob.shapes

  x_ls = np.linspace(0.5/float(N),1-(0.5/N),N)

  # this is generally not optimized but for understanding and validation
  dist_idx     = np.zeros((N,N,N),dtype=int) # part indices
  dist         = np.zeros((N,N,N))
  rho          = np.zeros((N,N,N))

  # this is using optimized code meant for cfs usage
  #assert len(glob.shapes) == 1
 # field_idx = create_idx_field(N) # has many shapes!

  total = sum([s.num_optvar for s in shapes])
  assert total == shapes[-1].base + shapes[-1].num_optvar

  for j, y in enumerate(x_ls):
    for i, x in enumerate(x_ls):
      for k, z in enumerate(x_ls):
        dist_min = 100*glob.transition
        for s in shapes:
          val, idx = s.dist([x,y,z])
          dist_min = np.min((dist_min,val))
          dist[i,j,k]     = dist_min
          dist_idx[i,j,k] = idx
          rho[i,j,k]      = boundary(dist_min)

  pd={"dist": dist}
  pd["dist_idx"]  = dist_idx
  pd["rho"] = rho

  gridToVTK(name, x_ls,x_ls,x_ls , pointData=pd)
  print("# wrote '" + name + ".vtr'")

# find the mininimum of (value, index) tuples, where only value is compared and value might contain None
def idx_min(a,b):
  if a[0] == None:
    return b
  else:
    if b[0] == None or a[0] < b[0]:
      return a
    else:
      return b

# reads 2D and returns list of Spaghettis
# @param radius if given overwrites the value from the xml header
def read_xml(filename, set = None):

  xml = ot.open_xml(filename)

  shapes = []
  sq = 'last()' if set == None else '@id="' + str(set) + '"'

  glob.n[0] = int(ot.xpath(xml, '//header/mesh/@x'))
  glob.n[1] = int(ot.xpath(xml, '//header/mesh/@y'))
  glob.n[2] = int(ot.xpath(xml, '//header/mesh/@z'))

  # this feature is only written by modern cfs (solar_heater, 12.2022) but we also want to read old files
  pn = xml.xpath('//header/coordinateSystems/domain')
  if len(pn) == 1: # we assume a single coordinate system
    dic = pn[0].attrib
    glob.box_lower = np.array([float(dic['min_x']),float(dic['min_y']),float(dic['min_z'])])
    glob.box_upper = np.array([float(dic['max_x']),float(dic['max_y']),float(dic['max_z'])])

  glob.dx = (glob.box_upper[0]-glob.box_lower[0])/glob.n[0]
  glob.mid_lower = glob.box_lower + 0.5*glob.dx*np.ones((3))

  while True: # exit with break
    idx = len(shapes)
    base = '//set[' + sq + ']/shapeParamElement[@shape="' + str(idx) + '"]'
    test = xml.xpath(base)
    if len(test) == 0:
      break

    # width of noodle is 2*w -> don't confuse with P
    p  = float(ot.xpath(xml, base + '[@type="profile"]/@design'))

    Px = read_xml_list(xml, base, '[@dof="x"]')
    Py = read_xml_list(xml, base, '[@dof="y"]')
    Pz = read_xml_list(xml, base, '[@dof="z"]')
    r = read_xml_list(xml, base, '[@type="radius"]')
    alpha = float(ot.xpath(xml, base + '[@type="alpha"]/@design'))

    base = sum([s.num_optvar for s in shapes])
    noodle = Spaghetti(idx, base, Px, Py, Pz, r, p, alpha)
    shapes.append(noodle)
    print('# read noodle', noodle)

  return shapes

def read_xml_list(xml, base, name):
  values = []
  last = -1
  list = xml.xpath(base + name)
  for entry in list:
    nr = int(entry.get('nr'))
    des = float(entry.get('design'))
    if not nr > last:
      raise('numbering for normal in shape ' + str(idx) + ' seems out of order: nr=' + str(nr))
    values.append(des)
    last = nr
  return values


def visualize(shapes, filename):
  # Create the graphics structure. The renderer renders into the render
  # window. The render window interactor captures mouse events and will
  # perform appropriate camera or actor manipulation depending on the
  # nature of the events.
  ren = vtkRenderer()
  renWin = vtkRenderWindow()
  renWin.AddRenderer(ren)
  iren = vtkRenderWindowInteractor()
  iren.SetRenderWindow(renWin)
  colors = vtkNamedColors()
  #print(colors)
  # Set the background color.
  mycolors = ['Red', 'LawnGreen', 'SteelBlue', 'MediumTurquoise', 'Coral', 'Pink', 'Lavender', 'Ivory']
  bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
  colors.SetColor("BkgColor", *bkg)
  for s in shapes:
    if s.alpha < 0.05:
      continue
    rgb = [0,0,0]
    colors.GetColorRGB(mycolors[s.id % len(mycolors)], rgb)

  # create cylinders
    for num, L in enumerate(s.L):
      L1 = s.L[num][0]
      L2 = s.L[num][1]
      # This creates a polygonal cylinder model with eight circumferential
      # facets.
      cylinder = vtkCylinderSource()
      cylinder.SetResolution(50)
      cylinder.SetRadius(s.profile)
      cylinder.SetCapping(False)
      vec = L2-L1
      height = norm(vec)
      vec = vec/height
      cylinder.SetHeight(height)
      rot = np.cross(vec,np.array((0,1,0)))
      angle = -np.arccos(np.dot(vec, (0,1,0)))*180/np.pi
      center = 0.5*(L1+L2)

      # The mapper is responsible for pushing the geometry into the graphics
      # library. It may also do color mapping, if scalars or other
      # attributes are defined.
      cylinderMapper = vtkPolyDataMapper()
      cylinderMapper.SetInputConnection(cylinder.GetOutputPort())

      # The actor is a grouping mechanism: besides the geometry (mapper), it
      # also has a property, transformation matrix, and/or texture map.
      # Here we set its color and rotate it -22.5 degrees.
      cylinderActor = vtkActor()
      cylinderActor.SetMapper(cylinderMapper)
      cylinderActor.GetProperty().SetColor(rgb)
      cylinderActor.GetProperty().SetOpacity((s.alpha**3)*.9)
      #cylinderActor.RotateX(90)
      cylinderActor.RotateWXYZ(angle, rot[0], rot[1], rot[2])
      cylinderActor.SetPosition(center)

      # Add the actors to the renderer, set the background and size
      ren.AddActor(cylinderActor)

      # Add spherical caps at the ends
      if num == 0:
        vec = s.t[0]
        # create spheres
        sphere = vtkSphereSource()
        #sphere.SetCenter(L1[0],L1[1],L1[2])
        sphere.SetRadius(s.profile)
        sphere.SetThetaResolution(20)
        sphere.SetPhiResolution(20)
        d = np.dot(vec, (0,0,1))
        if d < 0:
          vec = -vec
          sphere.SetEndPhi(90)
        else:
          sphere.SetStartPhi(90)
        rot = np.cross(np.array((0,0,1)),vec)
        angle = np.arccos(np.dot(vec, (0,0,1)))*180/np.pi
        #sphere.SetStartPhi(90)
        sphereMapper = vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())
        sphereActor = vtkActor()
        sphereActor.SetMapper(sphereMapper)
        sphereActor.GetProperty().SetColor(rgb)
        sphereActor.GetProperty().SetOpacity((s.alpha**3)*.9)
        sphereActor.RotateWXYZ(angle, rot[0], rot[1], rot[2])
        sphereActor.SetPosition(s.P[0])
        ren.AddActor(sphereActor)
      if num == s.m-1:
        vec = -s.t[-1]
        # create spheres
        sphere = vtkSphereSource()
        #sphere.SetCenter(L1[0],L1[1],L1[2])
        sphere.SetRadius(s.profile)
        sphere.SetThetaResolution(20)
        sphere.SetPhiResolution(20)
        d = np.dot(vec, (0,0,1))
        if d < 0:
          vec = -vec
          sphere.SetEndPhi(90)
        else:
          sphere.SetStartPhi(90)
        rot = np.cross(np.array((0,0,1)),vec)
        angle = np.arccos(np.dot(vec, (0,0,1)))*180/np.pi
        #sphere.SetStartPhi(90)
        sphereMapper = vtkPolyDataMapper()
        sphereMapper.SetInputConnection(sphere.GetOutputPort())
        sphereActor = vtkActor()
        sphereActor.SetMapper(sphereMapper)
        sphereActor.GetProperty().SetColor(rgb)
        sphereActor.GetProperty().SetOpacity((s.alpha**3)*.9)
        sphereActor.RotateWXYZ(angle, rot[0], rot[1], rot[2])
        sphereActor.SetPosition(s.P[-1])
        ren.AddActor(sphereActor)

    for C, gamma, L1, L2, r in s.arc:
      v1 = L1-C
      v2 = L2-C
      vec = np.cross(v1, v2)
      vec = vec/norm(vec)
      rot = np.cross(vec,np.array((0,0,1)))
      rot = rot/norm(rot)
      angle = np.arccos(np.dot(vec, (0,0,1)))

      # rotate border vectors of torus segment to default torus using Rodrigues' formula
      v1rot = v1*np.cos(angle) + (np.cross(rot, v1))*np.sin(angle) + rot*np.dot(rot, v1)*(1-np.cos(angle))
      v2rot = v2*np.cos(angle) + (np.cross(rot, v2))*np.sin(angle) + rot*np.dot(rot, v2)*(1-np.cos(angle))
      # get angle segment on default vtk torus
      u1 = np.arctan2(v1rot[1],v1rot[0])-.5*np.pi
      u2 = np.arctan2(v2rot[1],v2rot[0])-.5*np.pi
      gamma1 = min(u1,u2)
      gamma2 = max(u1,u2)
      if (gamma2-gamma1 > np.pi): # for edge case of arctan2 domain, angle difference should never be this large
        u1 = gamma2
        u2 = gamma1+2*np.pi

      fn = vtkParametricTorus()
      fn.SetRingRadius(r)
      fn.SetCrossSectionRadius(s.profile)
      # torus is parametrized clockwise, so change sign...
      fn.SetMinimumU(-u2)
      fn.SetMaximumU(-u1)
      fn.SetJoinU(False)

      torus = vtkParametricFunctionSource()
      torus.SetParametricFunction(fn)
      torus.SetUResolution(20)
      torus.SetVResolution(20)
      torus.Update()

      torusMapper = vtkPolyDataMapper()
      torusMapper.SetInputConnection(torus.GetOutputPort())
      torusActor = vtkActor()
      torusActor.GetProperty().SetColor(rgb)
      torusActor.GetProperty().SetOpacity((s.alpha**3)*.9)
      # RotateWXYZ in degree
      angle = -angle*180/np.pi
      torusActor.RotateWXYZ(angle, rot[0], rot[1], rot[2])
      torusActor.SetPosition(C)
      #print(fn)

      torusActor.SetMapper(torusMapper)




      ren.AddActor(torusActor)


    ren.AddActor(sphereActor)
  colors.GetColorRGB("BkgColor", rgb)
  ren.SetBackground(rgb)
  renWin.SetSize(1600, 1200)
  renWin.SetWindowName('CylinderExample')

  axes = vtkAxesActor()
  axes.AxisLabelsOff()

  transform = vtkTransform()
  transform.Translate(0.0, 0.0, 0.0)

  #  The axes are positioned with a user transform
  axes.SetUserTransform(transform)

  widget = vtkOrientationMarkerWidget()
  rgba = [0] * 4
  colors.GetColor('Carrot', rgba)
  widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
  widget.SetOrientationMarker(axes)
  widget.SetInteractor(iren)
  widget.SetViewport(0.0, 0.0, 0.4, 0.4)
  widget.SetEnabled(1)
  widget.InteractiveOn()

  # This allows the interactor to initalize itself. It has to be
  # called before an event loop.
  iren.Initialize()

  # We'll zoom in a little by accessing the camera and invoking a "Zoom"
  # method on it.
  ren.ResetCamera()
  ren.GetActiveCamera().Zoom(.8)
  renWin.Render()

  exporter = vtkOBJExporter()
  exporter.SetRenderWindow( renWin )
  exporter.SetFilePrefix(os.getcwd()+'/'+filename) #create mtl and obj file.
  exporter.Write()
  exporter = vtkX3DExporter()
  exporter.SetRenderWindow( renWin )
  exporter.SetFileName(os.getcwd()+'/'+filename+'.x3d') #create mtl and obj file.
  exporter.Write()

  # Start the event loop.
  print('Starting visualization')
  iren.Start()

      #points = vtkPoints()
      #points.InsertNextPoint(v1rot)
      #points.InsertNextPoint(v2rot)
      #polydata = vtkPolyData()
      #polydata.SetPoints(points)
      # Create anything you want here, we will use a cube for the demo.
      #arrowSource = vtkSphereSource()
      #glyph3D = vtkGlyph3D()
      #glyph3D.SetSourceConnection(arrowSource.GetOutputPort())
      #glyph3D.SetInputData(polydata)
      #glyph3D.Update()
      # Visualize
      #mapper = vtkPolyDataMapper()
      #mapper.SetInputConnection(glyph3D.GetOutputPort())
      #actor = vtkActor()
      #actor.SetMapper(mapper)
      #ren.AddActor(actor)

# reads 2D and returns list of Spaghetti
# @param radius if given overwrites the value from the xml header
def write_xml(filename, shapes, remove_ghosts=-1, remove_alpha=0, padnormals=0):
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <mesh x="' + str(glob.n[0]) + '" y="' + str(glob.n[1]) + '" z="' + str(glob.n[2]) + '"/>\n')
  out.write('    <spaghetti radius="-1"/>\n')
  out.write('  </header>\n')
  out.write('  <set id="spaghetti.py">\n')

  shapeid = 0
  nr = 0
  for shape in shapes:
    if shape.alpha <= remove_ghosts+1e-16:
      continue
    for num in range(shape.P.shape[0]):
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="node" dof="x" tip="' + ('start' if (num == 0) else ('end' if (num == shape.P.shape[0]-1) else "inner")) + '" shape="' + str(shapeid) + '" design="' + str(shape.P[num][0]) + '"/>\n')
      out.write('    <shapeParamElement nr="' + str(nr+1) + '" type="node" dof="y" tip="' + ('start' if (num == 0) else ('end' if (num == shape.P.shape[0]-1) else "inner")) + '" shape="' + str(shapeid) + '" design="' + str(shape.P[num][1]) + '"/>\n')
      out.write('    <shapeParamElement nr="' + str(nr+2) + '" type="node" dof="z" tip="' + ('start' if (num == 0) else ('end' if (num == shape.P.shape[0]-1) else "inner")) + '" shape="' + str(shapeid) + '" design="' + str(shape.P[num][2]) + '"/>\n')
      nr += 3
    for idx, r in enumerate(shape.radii[1:-1]):
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="radius" shape="' + str(shapeid) + '" design="' + str(r) + '"/>\n')
      nr += 1
    for i in range(len(shape.radii)-2, padnormals):
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="normal" shape="' + str(shapeid) + '" design="0.000000001"/>\n')
      nr += 1
    if remove_alpha == 0:
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="profile" shape="' + str(shapeid) + '" design="' + str(shape.profile) + '"/>\n')
      out.write('    <shapeParamElement nr="' + str(nr+1) + '" type="alpha" shape="' + str(shapeid) + '" design="' + str(shape.alpha) + '"/>\n')
      nr += 2
    else:
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="profile" shape="' + str(shapeid) + '" design="' + (str(shape.profile) if remove_alpha == 1 else str(shape.profile*shape.alpha)) + '"/>\n')
      out.write('    <shapeParamElement nr="' + str(nr+1) + '" type="alpha" shape="' + str(shapeid) + '" design="1"/>\n')
      nr += 2
    shapeid += 1


  out.write('  </set>\n')
  out.write(' </cfsErsatzMaterial>\n')
  out.close()



# __name__ is 'spaghetti' if imported or '__main__' if run as commandline tool
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml")
  parser.add_argument("--set", help="set within a .density.file", type=int)
  parser.add_argument("--order", help="number of integration points per dimension and element", choices=[1, 2], default=1, type=int)  
  parser.add_argument("--silent", help="verbosity of print output", action='store_true')
  parser.add_argument('--vtk', help="write vtk file for given name (w/o extenstion)")
  parser.add_argument('--vtk_res', help="resolution for vtk export", type=int, default=200)
  parser.add_argument('--vtk_lists', help="write vtk file for given name (w/o extention)")
  parser.add_argument('--noshow', help="omit visualization", action='store_true')
  parser.add_argument('--gradient_check', help="check internal spaghetti derivatives", action='store_true')
  parser.add_argument('--track', help="track given density file to represent topology using spaghetti shapes")
  parser.add_argument('--track_iter', help="number of optimization iterations for tracking", type=int, default=30)
  parser.add_argument('--remove_ghosts', help="load .density.xml file and remove ghost shapes with alpha <= value", type=float, default=-1)
  parser.add_argument('--remove_alpha', help="remove alpha variable from .density.xml file: 0=off, 1=remove, 2=scale profile with previous alpha", choices=[0,1,2], default=0, type=int)
  parser.add_argument('--saveall', help="save all sets as image with the given format. Might be png, pdf, eps, vtp", action='store_true')
  parser.add_argument('--padradii', help="pad radii with additional radii (zero-valued) so there is a total of x radii", type=int, default=0)


  args = parser.parse_args()

  glob.shapes = read_xml(args.input, args.set)
  glob.order = args.order
  glob.rhomax = 1
  glob.rhomin = 0
  glob.gradient_check = args.gradient_check
  glob.silent = args.silent

  if args.vtk_lists:
    cfs_map_to_design()
    write_vtk_lists(args.vtk_lists)
  if args.vtk:
    write_vtk(args.vtk, args.vtk_res)

  if args.gradient_check:
    eps = np.sqrt(np.finfo(float).eps)
    print('checking gradient of primitives with FD step: ' + str(eps))
    for shape in glob.shapes:
      optvar = shape.optvar()
      which =  []
      for st in ['P', 't', 'gamma', 'C', 'n0_arc']:
        for ii in range(len(shape.t)+1):
          if (st == 't' and ii == len(shape.t)) or (st == 'gamma' and (ii == 0 or ii == len(shape.radii)-1)) or (st == 'n0_arc' and (ii == 0 or ii == shape.m)):
            #print('skip ' + st + str(ii))
            continue
          which.append(st+str(ii))
      for qname in which:
        for xyz in [0,1,2]:
          if qname.startswith('gamma'):
            if xyz == 1 or xyz == 2:
              continue
          fd_grad = sciopt.approx_fprime(optvar, shape.func, eps, [xyz, qname])
          grad = shape.grad(optvar, [xyz, qname])
          if norm(fd_grad-grad) > 1e-6:
            print(qname + ',variable:' + str(xyz))
            print("FD approximation:" + str(fd_grad))
            print("Analytic gradient:" + str(grad))
            print("L2 error:" + str(norm(fd_grad-grad)))
    print('\n\n')

#    n  = glob.n
#    nx = int(n[0])
#    ny = int(n[1])
#    nz = int(n[2])
#    dx = glob.dx
#    print('checking gradient of distances with FD step: ' + str(eps))
#    for shape in glob.shapes:
#      for i in range(nx):
#        for j in range(ny):
#          for k in range(nz):
#            X = [(i+0.5) * dx, (j+0.5)*dx, (k+0.5)*dx]
#           #tt = df(var,i, j,derivative='autograd') # autograd gives no nice vector
#            grad = shape.grad(optvar, [X, 'distance'])
#            fd_grad = sciopt.approx_fprime(optvar, shape.func, eps, [X, 'distance'])
#            #print(grad)
#            #print(fd_grad)
#            if norm(fd_grad-grad) > 1e-6:
#              niceidx = shape.nice_idx(shape.dist(X, what='index'))
#              print('X:', X, 'idx:', niceidx, 'with', shape.m, 'segments')
#              print("FD approximation:" + str(fd_grad))
#              print("Analytic gradient:" + str(grad))
#              print("L2 error:" + str(norm(fd_grad-grad)))
#    print('\n\n')

    assemble_fe_field()
    integrate_fe_field()
    n = glob.n
    count = 0
    opt = np.zeros(glob.total())
    for s in glob.shapes:
      opt[s.base:s.base+s.num_optvar] = s.optvar()
    with Pool(glob.num_threads) as pool:
      for i in range(n[0]):
        for j in range(n[1]):
          for k in range(n[2]):
            Xmid = glob.mid_lower + np.array([i,j,k],dtype=int)*glob.dx
  #          for snum, s in enumerate(glob.shapes):
  #            if s.dist(Xmid)[0] > glob.h + glob.dx:
  #              continue
  #            fd_grad = sciopt.approx_fprime(optvar, s.func, eps, [Xmid, 'boundary'])
  #            if norm(fd_grad-glob.grad_boundary[i,j,k,snum][s.base:s.base+s.num_optvar]) > 1e-4:
  #              count += 1
  #              #niceidx = shape.nice_idx(shape.dist(Xmid, what='index'))
  #              d = s.dist(Xmid)
  #              print('X:', Xmid, 'dist:', d, 'rho:', boundary(d))#, 'idx:', niceidx, 'with', shape.m, 'segments')
  #              for shape_idx_list in glob.fe_list[i][j][k].shapes_idx_list:
  #                dist, gradd = glob.shapes[shape_idx_list[0]].dist_by_idx(Xmid, shape_idx_list[1:])
  #                rho = boundary([dist, gradd], derivative=True)
  #                print(shape_idx_list, 'dist_by_idx:', dist, 'rho:', rho)
  #              print("FD approximation:" + str(fd_grad))
  #              print("Analytic gradient:" + str(glob.grad_boundary[i,j,k,snum][s.base:s.base+s.num_optvar]))
  #              print("L2 error:" + str(norm(fd_grad-glob.grad_boundary[i,j,k,snum][s.base:s.base+s.num_optvar])))
            rho_ref_s = np.array([s.func(opt[s.base:s.base+s.num_optvar], [Xmid, 'boundary']) for s in glob.shapes])
            w_s = np.exp(glob.p*rho_ref_s)
            rho_ref = np.sum(w_s*rho_ref_s)/np.sum(w_s)
            if rho_ref < 1e-16:
              continue
            for optvar in range(glob.total()):
              delta = np.zeros(glob.total())
              delta[optvar] = eps
              #rho_sm = boundary(s.func(s.optvar()+delta, [Xmid, 'distance']))
              rho_s = np.array([s.func(opt[s.base:s.base+s.num_optvar]+delta[s.base:s.base+s.num_optvar], [Xmid, 'boundary']) for s in glob.shapes])
              w_s = np.exp(glob.p*rho_s)
              rho_sm = np.sum(w_s*rho_s)/np.sum(w_s)
              fd = (rho_sm - rho_ref)/eps
              #print(rho_sm, glob.rho[i,j,k], fd, glob.grad_rho[i,j,k][optvar])
              err = abs(fd-glob.grad_rho[i,j,k][optvar])
              if err > 5e-5:
                count += 1
                print('FD error at point', Xmid, 'for variable', optvar, 'is', err, 'with gradients', fd, glob.grad_rho[i,j,k][optvar], 'and values', rho_ref, rho_sm, glob.rho[i,j,k])
                for shape_idx_list in glob.fe_list[i][j][k].shapes_idx_list:
                  print('shape_idx_list:', shape_idx_list, 'with dist', glob.shapes[shape_idx_list[0]].dist(Xmid))
                print('Number of full material shapes:', len(glob.fe_list[i][j][k].shapes_full))
      pool.close()
      pool.join()
    print('number of false gradients:', count)

  if args.remove_alpha or args.remove_ghosts > -.1:
    write_xml(args.input[0:-12] + '_modified.density.xml', glob.shapes, args.remove_ghosts, args.remove_alpha, padnormals=args.padradii)

  if args.track:
    to_track = ot.read_density(args.track, attribute="design", fill=1.0)
    res = track(glob.shapes, to_track, args.saveall, args.track_iter)
    write_xml(args.track[0:-12] + '_tracked.density.xml', glob.shapes, args.remove_ghosts, padnormals=args.padradii)

  if not args.noshow:
    visualize(glob.shapes, args.input[:-12])


