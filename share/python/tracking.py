"""Reusable feature-mapping objectives for python optimization functions (embedded python).

Tracking fits the aggregated feature pseudo density rho to a target field rho* by least squares
  J_track  = sum_e (rho_e - rho*_e)^2,   dJ_track/drho_e  = 2 (rho_e - rho*_e)
Reward (paper eq. reward_obj) maximizes the overlap with the target region; it is linear in rho, so
its curvature is 0 (first order only)
  J_reward = -sum_e rho*_e rho_e,        dJ_reward/drho_e = -rho*_e
Both read the target field via track_init; select track_* or reward_* as eval/grad in the xml.

For second-order optimization the per-element objective curvature d^2J/drho_e^2 is constant: 2 for
tracking, 0 for reward (track_curvature / reward_curvature). cfs holds the jacobian D = drho/ds and
forms the objective Hessian term D^T diag(curv) D itself; the geometric d^2rho/ds^2 term is added by
cfs from the boundary function, so a C^2 boundary (quintic/bezier) is required. With reward (curv=0)
the D^T diag(curv) D term vanishes, so the exported Hessian is purely the geometric/aggregation part.
"""

import numpy as np
import cfs
import optimization_tools as ot
# the target density field rho*, read once from opt['target'] in track_init
target = None


def track_init(opt):
  global target
  # read the flat target field in file/element order, matching get_pseudo_density
  target = np.array(ot.read_density_as_vector(opt['target']))
  print('tracking.py: read target', opt['target'], 'with', len(target), 'elements')

def _density():
  ne = cfs.get_num_pseudo_density()
  rho = np.zeros(ne)
  cfs.get_pseudo_density(rho)   # aggregated element density mrho_e
  return rho


def track_eval(opt):
  rho = _density()
  return float(np.sum((rho - target) ** 2))


def track_grad(opt):
  rho = _density()
  return 2.0 * (rho - target)


def reward_eval(opt):
  rho = _density()
  return float(-np.sum(target * rho))


def reward_grad(opt):
  # dJ_reward/drho_e = -rho*_e (constant: reward is linear in rho; cfs chains it through drho/ds)
  return -target


def track_curvature(opt):
  # per-element d^2 J_track/drho_e^2 = 2 (constant), returned like the gradient
  return np.full(cfs.get_num_pseudo_density(), 2.0)


def reward_curvature(opt):
  # reward is linear in rho, so d^2 J_reward/drho_e^2 = 0
  return np.zeros(cfs.get_num_pseudo_density())


# The alpha sum J = sum_f alpha_f is the pruning objective (minimal feature material) for the
# formulation min sum alpha s.t. tracking <= eps. It lives directly on the shape variables
# (design="feature" in the xml): eval/grad work in the full feature variable space
# [Px Py Qx Qy p alpha] per pill as by cfs.feature_mapping_get_parameters(). The objective is linear,
# so its exact Hessian is zero; the curvature (in density space) is zero as well - with it the cfs
# shape Hessian machinery correctly returns an all-zero objective Hessian.
NV_ALPHA = 6  # variables per pill when the geometry variable alpha is active

def _alpha_indices():
  n = cfs.feature_mapping_num_parameters()
  assert n % NV_ALPHA == 0, 'alpha_sum expects %d variables per pill (with alpha)' % NV_ALPHA
  return n, np.arange(NV_ALPHA - 1, n, NV_ALPHA)


def alpha_sum_eval(opt):
  n, idx = _alpha_indices()
  z = np.zeros(n)
  cfs.feature_mapping_get_parameters(z)
  return float(np.sum(z[idx]))


def alpha_sum_grad(opt):
  n, idx = _alpha_indices()
  g = np.zeros(n)
  g[idx] = 1.0
  return g


def alpha_sum_curvature(opt):
  # linear in the design, no density dependence at all
  return np.zeros(cfs.get_num_pseudo_density())
