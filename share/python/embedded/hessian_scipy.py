"""Reusable second-order python optimizer for feature mapping (PythonOptimizer interface).

Uses scipy.optimize.minimize (trust-constr, handles provided hessian). Don't use for real work but
prefer cyipopt via hessian_ipopt.py. The scipy variant is for the pipeline testsuite

For a state-INDEPENDENT objective (reward/tracking) cfs.evalhessian is the complete exact Hessian and
used as is. For a state-DEPENDENT objective wrapped as a python function (compliance via
observed_objective.py) cfs.evalhessian is only the geometric part D^T diag(curv) D + H_agg + H_feat;
the dense state-curvature block D^T (2 v^T K^-1 v) D is missing. Given the wrapped function's name via
the 'observation' option it is added EXACTLY via observed_objective.state_hessian(n) (n extra
back-substitutions against the factorized state matrix, self-adjoint functions) - giving the true
Newton Hessian. Only if that is unavailable it is learned in shape space by a structured,
Powell-damped BFGS B on the geometric-corrected secant y_struct = (grad_k+1 - grad_k) - H_geo(x_k+1) s.
Auto-enabled iff that function is state-dependent (cfs.is_function_state_dependent); without the
option the classic exact-Hessian path is unchanged.

Options (<option key=".." value=".."/> in the optimizer's python element):
  observation : name of the state-dependent cfs function the objective wraps -> enables the exact
                state block (BFGS fallback). Omit for a directly exact objective.
  exact_state : "false" -> skip the exact state block and use the structured BFGS (testing/comparison)
  fd_check : "true" -> central-difference-check the analytic gradient against cfs.evalobj and the
             analytic Hessian against a central difference of the gradient at the initial design,
             asserting agreement (1st and 2nd order). Only in the BFGS-fallback mode the objective
             Hessian is quasi-Newton and just the gradient is asserted. Constraint Hessians are
             checked per row against a central difference of the Jacobian; rows without exact
             second order info in cfs (zero contribution, e.g. state-dependent observations) are
             skipped.
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

  # state block handling for a wrapped state-dependent objective (see module docstring): exact via
  # observed_objective.state_hessian when available, else structured BFGS. Auto-enabled iff the
  # named function is state-dependent; without the option this stays a no-op (classic path).
  obs = options.get('observation')
  glob.structured = bool(obs) and bool(cfs.is_function_state_dependent(obs))
  glob.exact = False
  if glob.structured and str(options.get('exact_state', 'true')).lower() == 'true':
    try:
      import observed_objective
      glob.exact = hasattr(observed_objective, 'state_hessian')
    except ImportError:
      pass
  if obs:
    print('hessian.py: wrapped function %r state-dependent=%s -> state block %s'
          % (obs, glob.structured, 'exact' if glob.exact else ('structured BFGS' if glob.structured else 'none')))
  glob.B = None          # learned state block (BFGS fallback)
  glob.prev = None       # (x, grad) of the last gradient evaluation
  glob.pending = None    # secant candidate (s, y) waiting for hess() to correct with H_geo
  glob.collect = False   # gather secant pairs (only during solve, not fd_check)
  glob.scaled = False    # B rescaled from the first accepted pair


def eval(x):
  return cfs.evalobj(x)

def grad(x):
  g = np.zeros(glob.n)
  cfs.evalgradobj(x, g)
  if glob.collect:
    if glob.prev is not None and not np.array_equal(x, glob.prev[0]):
      glob.pending = (x - glob.prev[0], g - glob.prev[1])
    glob.prev = (x.copy(), g.copy())
  return g


def bfgs_update(B, s, y):
  """damped BFGS update (Powell): keeps B positive definite for indefinite secant data"""
  Bs = B @ s
  sBs = float(s @ Bs)
  sy = float(s @ y)
  if sy < 0.2 * sBs:
    theta = 0.8 * sBs / (sBs - sy)
    y = theta * y + (1.0 - theta) * Bs
    sy = float(s @ y)
  if sy <= 1e-12 * max(1.0, float(s @ s)):
    return B  # degenerate pair, skip
  return B - np.outer(Bs, Bs) / sBs + np.outer(y, y) / sy


def hess(x):
  H = np.zeros((glob.n, glob.n))
  cfs.evalhessian(x, H)  # exact geometric terms (objective curvature); complete when not structured

  if not glob.structured:
    return H

  if glob.exact:
    # exact dense state-curvature block 2 B^T Z: n back-substitutions against the factorized
    # state matrix, evaluated at the design cfs.evalhessian just synced
    import observed_objective
    return H + observed_objective.state_hessian(glob.n)

  # learn the missing dense state-curvature block: B on the geometric-corrected secant, see docstring
  if glob.B is None:
    glob.B = max(1e-3, abs(np.trace(H)) / glob.n) * np.eye(glob.n)
  if glob.pending is not None:
    s, y = glob.pending
    glob.pending = None
    y_struct = y - H @ s
    sy = float(s @ y_struct)
    if not glob.scaled and sy > 0.0:
      glob.B = (float(y_struct @ y_struct) / sy) * np.eye(glob.n)  # first-pair scaling
      glob.scaled = True
    glob.B = bfgs_update(glob.B, s, y_struct)
  return H + glob.B


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
  if glob.structured and not glob.exact:
    H = np.zeros((glob.n, glob.n))
    cfs.evalhessian(x, H)  # geometric part only, |Hfd-H| below shows the state block to be learned
  else:
    H = hess(x)  # full objective Hessian (exact everywhere)
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
  if glob.structured and not glob.exact:
    # the objective Hessian is quasi-Newton here; |Hfd-Hgeo| is the missing state block (informative)
    print('hessian.py: structured - |Hfd-Hgeo|=%.2e is the learned state block, not asserted' % herr)
  else:
    # exact everywhere: plain exact Hessian or geometric terms + exact state block. The absolute
    # error scales with |H| (e.g. |H|~1e4 for compliance), so assert relative to it.
    hscale = max(1.0, float(np.max(np.abs(Hfd))))
    assert herr / hscale < tol, 'Hessian finite difference error %.3e (rel %.3e) exceeds %.3e' % (herr, herr / hscale, tol)

  # constraint Hessians: FD the Jacobian once, verify every row cfs has exact second order info for
  # (local sparse like 'distance', curvature based like the native volume, or the 'hessian' attribute
  # of a python function like pill_volume.py). A zero contribution means no info in cfs -> skip.
  if glob.m > 0 and glob.nnz > 0:
    Jfd = np.zeros((glob.m, glob.n, glob.n))
    for j in range(glob.n):
      xp = x.copy(); xp[j] += eps
      xm = x.copy(); xm[j] -= eps
      Jfd[:, :, j] = (constr_jac(xp).toarray() - constr_jac(xm).toarray()) / (2 * eps)
    for r in range(glob.m):
      v = np.zeros(glob.m); v[r] = 1.0
      Hc = constr_hess(x, v)
      if not np.any(Hc):
        continue
      cerr = float(np.max(np.abs(Hc - Jfd[r])))
      cscale = max(1.0, float(np.max(np.abs(Jfd[r]))))
      print('hessian.py: fd_check constraint %d max|hess-fd|=%.2e' % (r, cerr))
      assert cerr / cscale < tol, 'constraint %d Hessian finite difference error %.3e (rel %.3e) exceeds %.3e' % (r, cerr, cerr / cscale, tol)


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

  glob.collect = glob.structured and not glob.exact  # secant pairs only for the BFGS fallback

  # trust-constr tolerances; defaults are tight, the testsuite loosens them via <option> for a quick run
  opts = {'maxiter': glob.maxiter,
          'gtol': float(glob.options.get('gtol', 1e-7)),
          'xtol': float(glob.options.get('xtol', 1e-9)),
          'verbose': int(glob.options.get('verbose', 0))}
  # further numeric trust-constr options pass through verbatim. Notably initial_barrier_parameter
  # (default 0.1): a larger value (e.g. 1.0) keeps the early iterates away from the variable bounds,
  # which can be the difference between a clean KKT point and a spurious xtol stall (see the
  # hessian_alpha pruning test, where the default collapses the trust region on the compensation ridge).
  for k in ('barrier_tol', 'initial_tr_radius', 'initial_barrier_parameter',
            'initial_barrier_tolerance', 'initial_constr_penalty'):
    if k in glob.options:
      opts[k] = float(glob.options[k])
  res = minimize(eval, x0, jac=grad, hess=hess, bounds=bounds,
                 constraints=cons, method='trust-constr',
                 options=opts,
                 callback=lambda xk, *rest: commit(xk))  # trust-constr: callback(xk, state)
  print('hessian.py: solve done J=%.6e nit=%d status=%d %s' % (res.fun, res.nit, res.status, res.message))
  commit(res.x)  # export trust-constr's endpoint (the iteration cap is set in the .xml, see there)
  return 0
