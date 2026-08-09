"""Generic kernel: turn any observed C++ cfs function into a python objective read on demand.

The wrapped function is declared as <constraint type="..." mode="observation" observeGradient="true"/>
and named via the 'observation' option. Value and per-element density gradient are read at the
current design/state through cfs.eval_function(name) / cfs.eval_function_gradient(name, arr) - not
from the commit-logging cache of get_opt_function_values. The curvature callback returns zeros: for a
state-dependent function (e.g. compliance) the d^2J/drho^2 is dense (not per-element diagonal), so
cfs.evalhessian yields only the exact geometric Hessian terms (aggregation + feature).

For a SELF-ADJOINT wrapped function the missing dense state-curvature block is available exactly via
state_hessian(n): with the shape Jacobian D = d_mrho/d_s, the pseudo loads b_i = sum_e D_ei
(dK_e/drho_e) u (cfs.apply_dk_drho) and the state derivative solves K z_i = b_i against the already
factorized state matrix (cfs.solve_state), the block is 2 B^T Z - N_s extra back-substitutions, no
new FE system. The drivers hessian_scipy.py / hessian_ipopt.py call it automatically in their
'observation' mode; without these apis they fall back to a structured BFGS.

Select via <costFunction type="python"><python init="obj_init" eval="obj_eval" grad="obj_grad"
curvature="obj_curvature" script="kernel"><option key="observation" value="<function name>"/> ... .
"""

import numpy as np
import cfs

class _G:
  pass
glob = _G()


def obj_init(opt):
  glob.name = opt.get('observation', 'compliance')  # name of the observed C++ function
  print('observed_objective.py: wrapping observation', glob.name)

def obj_eval(opt):
  return cfs.eval_function(glob.name)

def obj_grad(opt):
  g = np.zeros(cfs.get_num_pseudo_density())
  cfs.eval_function_gradient(glob.name, g)
  return g

def obj_curvature(opt):
  # zero diagonal curvature -> cfs.evalhessian returns the exact geometric terms only, the dense
  # state block comes exactly from state_hessian() (or is learned by the driver's BFGS fallback)
  return np.zeros(cfs.get_num_pseudo_density())


def state_hessian(n):
  """exact state-curvature Hessian block 2 B^T Z of the wrapped self-adjoint function.

  Evaluate at the current design/state (i.e. right after cfs.evalhessian in the driver's Hessian
  callback). n is the number of optimization variables. Cost: n sparse pseudo-load assemblies and
  n back-substitutions against the factorized state matrix."""
  ne = cfs.get_num_pseudo_density()
  D = np.zeros((ne, n))
  cfs.get_shape_jacobian(D)

  nalg = int(np.asarray(cfs.fe_function_total_equations(0)).ravel()[0])
  B = np.zeros((n, nalg))
  Z = np.zeros((n, nalg))
  for i in range(n):
    cfs.apply_dk_drho(glob.name, D[:, i], B[i])
    cfs.solve_state(glob.name, B[i], Z[i])

  H = 2.0 * (B @ Z.T)
  return 0.5 * (H + H.T)  # b_i^T z_j = z_i^T K z_j is symmetric, clean up round-off
