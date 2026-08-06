"""Analytic pill (capsule) volume as a python constraint directly on the feature variables.

This is the linear-in-alpha overlap variant (Norato 2015, eq. 28): the capsule volumes, weighted
by the geometry variable alpha if present (eta_v = 1, while the stiffness sees alpha^q - the
mismatch drives alpha to 0/1)
  2D: v_f = a_f (2 p_f l_f + pi p_f^2)          (stadium area)
  3D: v_f = a_f (pi p_f^2 l_f + 4/3 pi p_f^3)   (capsule volume)
Overlaps are counted multiply (conservative, relative error ~ p/l per joint); no density field is
involved, so the p-norm overshoot and the alpha penalization do not enter. For the unit
square/cube the value is directly the volume fraction.

The code is transparent in the dimension and the presence of alpha via cfs.feature_mapping_layout().
The variable order per pill is [Px Py (Pz) Qx Qy (Qz) p (alpha)] as by
cfs.feature_mapping_get_parameters(); the gradient (vol_grad) and the Hessian (vol_hessian) are
returned in the same full feature variable space - cfs recognizes the space by the size (gradient)
respectively expects it ('hessian' attribute) and maps to the optimization variables without the
density chain.

Usage as (part of the) global python kernel, see the pill_volume, pill_volume_hessian and
pill3d_volume tests:
  <python name="volume" init="vol_init" eval="vol_eval" grad="vol_grad" hessian="vol_hessian" script="kernel"/>
"""

import numpy as np
import cfs


def _layout():
  """(features, vars_per_feature, dim, alpha) from cfs"""
  lo = cfs.feature_mapping_layout()
  return int(lo['features']), int(lo['vars_per_feature']), int(lo['dim']), lo['alpha'] == '1'


def _params():
  """all feature variables (including fixed ones) as a features x vars_per_feature matrix"""
  z = np.zeros(cfs.feature_mapping_num_parameters())
  cfs.feature_mapping_get_parameters(z)
  nf, nv, _, _ = _layout()
  return z.reshape(nf, nv)


def _coef(p, dim):
  """v_f = a (c(p) l + d(p)): the profile factors c, c', c'' and d, d', d''"""
  if dim == 2:
    return 2*p, 2.0, 0.0, np.pi*p*p, 2*np.pi*p, 2*np.pi
  return np.pi*p*p, 2*np.pi*p, 2*np.pi, 4.0/3.0*np.pi*p**3, 4*np.pi*p*p, 8*np.pi*p


def _split(row, dim, alpha):
  """P, Q, p, a of a feature row; a = 1 without alpha"""
  return row[:dim], row[dim:2*dim], row[2*dim], row[2*dim+1] if alpha else 1.0


def vol_init(opt):
  pass


def vol_eval(opt):
  _, _, dim, alpha = _layout()
  V = 0.0
  for row in _params():
    P, Q, p, a = _split(row, dim, alpha)
    c, _, _, d, _, _ = _coef(p, dim)
    V += a * (c * np.linalg.norm(Q - P) + d)
  return float(V)


def vol_grad(opt):
  nf, nv, dim, alpha = _layout()
  z = _params()
  g = np.zeros(z.size)
  for f, row in enumerate(z):
    P, Q, p, a = _split(row, dim, alpha)
    dPQ = Q - P
    l = np.linalg.norm(dPQ)
    c, cp, _, d, dp, _ = _coef(p, dim)
    gf = g[f*nv:(f+1)*nv]  # view into g
    if l > 1e-12:  # kink of l at P=Q; keep away via a lower bound distance constraint
      gf[:2*dim] = a * c * np.concatenate((-dPQ, dPQ)) / l
    gf[2*dim] = a * (cp * l + dp)
    if alpha:
      gf[2*dim+1] = c * l + d
  return g


def _pill_hessian(row, nv, dim, alpha):
  """nv x nv Hessian of v_f = a (c(p) l + d(p)) over [P Q p (a)]"""
  P, Q, p, a = _split(row, dim, alpha)
  dPQ = Q - P
  l = np.linalg.norm(dPQ)
  c, cp, cpp, d, dp, dpp = _coef(p, dim)
  H = np.zeros((nv, nv))
  ip, ia = 2*dim, 2*dim + 1  # index of profile and alpha
  H[ip, ip] = a * (cpp * l + dpp)
  if alpha:
    H[ip, ia] = H[ia, ip] = cp * l + dp  # d2v/da2 = 0, v is linear in alpha
  if l < 1e-12:  # see vol_grad
    return H
  s = np.concatenate((-dPQ, dPQ))     # dl/dr * l
  gl = s / l
  I = np.eye(dim)
  M = np.block([[I, -I], [-I, I]])    # ds/dr
  Hl = M / l - np.outer(s, s) / l**3  # d2l/dr2
  r = slice(0, 2*dim)
  H[r, r] = a * c * Hl
  H[r, ip] = H[ip, r] = a * cp * gl
  if alpha:
    H[r, ia] = H[ia, r] = c * gl
  return H


def vol_hessian(opt):
  """the exact Hessian in the full feature variable space (the 'hessian' attribute callback);
  block diagonal by pill, cfs reduces to the optimization variable space"""
  nf, nv, dim, alpha = _layout()
  z = _params()
  H = np.zeros((z.size, z.size))
  for f, row in enumerate(z):
    b = slice(f*nv, (f+1)*nv)
    H[b, b] = _pill_hessian(row, nv, dim, alpha)
  return H
