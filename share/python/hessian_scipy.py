"""Reusable second-order python optimizer for feature mapping (PythonOptimizer interface).

Uses scipy.optimize.minimize (trust-constr, handles provided hessian). Don't use for real work but
prefer cyipopt via hessian_ipopt.py. The scipy variant is for the pipeline testsuite

Options (<option key=".." value=".."/> in the optimizer's python element):
  fd_check : "true" -> central-difference-check the analytic gradient against cfs.evalobj and the
             analytic Hessian (cfs.evalhessian) against a central difference of the gradient at the
             initial design, asserting agreement (1st and 2nd order).
  fd_eps   : finite difference step (default 1e-6)
  fd_tol   : finite difference tolerance (default 1e-4)
"""

import numpy as np
import cfs

# organize global properties in a glob. instance
class _G:
  pass
glob = _G()

# called by cfs
def setup():
  pass

# called by cfs
def init(n, m, maxiter, sim_name, options):
  glob.n = n
  glob.m = m
  glob.maxiter = maxiter
  glob.options = options
  objs, cnstrs = cfs.get_opt_function_values()
  assert len(objs) == 1, "hessian.py expects a single objective"
  print('hessian.py: init n=',n,'m=',m,'objective=',list(objs.keys()),'constraints=',list(cnstrs.keys()))

  # sparse constraint jacobian: the structure is static, so query it once here (rows/cols come back
  # as doubles, cast to int for scipy.sparse). The per-iteration values follow in constr_jac.
  glob.nnz = cfs.get_num_jacobian_nonzeros() if m > 0 else 0
  if glob.nnz > 0:
    rows = np.zeros(glob.nnz)
    cols = np.zeros(glob.nnz)
    cfs.get_constraint_sparsity(rows, cols)
    glob.rows = rows.astype(int)
    glob.cols = cols.astype(int)


def eval(x):
  return cfs.evalobj(x)

def grad(x):
  g = np.zeros(glob.n)
  cfs.evalgradobj(x, g)
  return g

def hess(x):
  H = np.zeros((glob.n, glob.n))
  cfs.evalhessian(x, H)
  return H


def constr(x):
  g = np.zeros(glob.m)
  cfs.evalconstrs(x, g)
  return g

def constr_jac(x):
  # sparse constraint jacobian: cfs fills the packed values, we wrap them with the (rows, cols)
  # structure into a scipy.sparse matrix that trust-constr exploits.
  from scipy.sparse import coo_matrix
  vals = np.zeros(glob.nnz)
  cfs.evalgradconstrs_sparse(x, vals)
  return coo_matrix((vals, (glob.rows, glob.cols)), shape=(glob.m, glob.n))

def constr_hess(x, v):
  # multiplier-contracted constraint Hessian sum_c v_c * Hess(c_c)(x). Without it trust-constr
  # approximates the constraint part of the Lagrangian Hessian by BFGS, which mismatches the exact
  # objective Hessian once a constraint is active and makes the second-order steps drift off a better
  # feasible point (see the native 'distance' Hessian, Condition::CalcHessian).
  H = np.zeros((glob.n, glob.n))
  cfs.evalhessian_constr(x, v, H)
  return H


def commit(x):
  # We must evaluate at x first so cfs' internal design is x
  cfs.evalobj(x)
  cfs.commitIteration()

def fd_check(x):
  eps = float(glob.options.get('fd_eps', 1e-6))
  tol = float(glob.options.get('fd_tol', 1e-4))

  g = grad(x)
  H = hess(x)
  gfd = np.zeros(glob.n)
  Hfd = np.zeros((glob.n, glob.n))
  for j in range(glob.n):
    xp = x.copy(); xp[j] += eps
    xm = x.copy(); xm[j] -= eps
    gfd[j] = (eval(xp) - eval(xm)) / (2 * eps)
    Hfd[:, j] = (grad(xp) - grad(xm)) / (2 * eps)

  gerr = float(np.max(np.abs(g - gfd)))
  herr = float(np.max(np.abs(H - Hfd)))
  print('hessian.py: fd_check max|grad-fd|=%.2e max|hess-fd|=%.2e' % (gerr, herr))
  assert gerr < tol, 'gradient finite difference error %.3e exceeds %.3e' % (gerr, tol)
  assert herr < tol, 'Hessian finite difference error %.3e exceeds %.3e' % (herr, tol)


def solve():
  from scipy.optimize import minimize, NonlinearConstraint, Bounds

  x0 = np.zeros(glob.n)
  cfs.initialdesign(x0)
  # No manual iteration-0 commit: trust-constr fires the callback already at x0 (its iteration 0),
  # after it has evaluated gradient and Hessian there - so that callback commit is the complete
  # iteration-0 record (design variables + the exact objective Hessian). A manual evalobj-only commit
  # here would add an incomplete duplicate (same x0, but no objective Hessian since only the value,
  # not the gradient, was computed yet).

  if str(glob.options.get('fd_check', 'false')).lower() == 'true':
    fd_check(x0)

  xl = np.zeros(glob.n); xu = np.zeros(glob.n)
  gl = np.zeros(glob.m); gu = np.zeros(glob.m)
  cfs.bounds(xl, xu, gl, gu)

  bounds = Bounds(xl, xu, keep_feasible=True)

  cons = []
  if glob.m > 0:
    # native cfs constraints (e.g. the per-feature 'distance' pill length), raw values vs raw bounds.
    # hess=constr_hess provides the exact constraint Hessian so the full Lagrangian Hessian is exact.
    cons.append(NonlinearConstraint(constr, gl, gu, jac=constr_jac, hess=constr_hess))

  # trust-constr tolerances; defaults are tight, the testsuite loosens them via <option> for a quick run
  gtol = float(glob.options.get('gtol', 1e-7))
  xtol = float(glob.options.get('xtol', 1e-9))
  res = minimize(eval, x0, jac=grad, hess=hess, bounds=bounds,
                 constraints=cons, method='trust-constr',
                 options={'maxiter': glob.maxiter, 'gtol': gtol, 'xtol': xtol, 'verbose': 0},
                 callback=lambda xk, *rest: commit(xk))  # trust-constr: callback(xk, state)
  print('hessian.py: solve done J=%.6e %s' % (res.fun, res.message))
  commit(res.x)  # export trust-constr's endpoint (the iteration cap is set in the .xml, see there)
  return 0
