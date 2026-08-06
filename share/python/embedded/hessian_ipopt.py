"""Drop-in cyipopt replacement for hessian.py (same PythonOptimizer interface and cfs.* calls).

Where hessian.py drives scipy trust-constr, this drives Ipopt via cyipopt with the *exact* Lagrangian
Hessian: the objective Hessian cfs.evalhessian plus the multiplier-contracted constraint Hessian
cfs.evalhessian_constr. The native 'distance' constraints come from cfs (sparse Jacobian + exact
Hessian), exactly as in hessian.py.

This is a developer/comparison driver: cyipopt is NOT available on the CI runners, so the testsuite
keeps hessian.py. Swap <python file="hessian.py"> -> <python file="ipopt_hession.py"> to compare. Unlike
hessian.py it does NOT keep the best feasible iterate - it exports Ipopt's actual final solution, so the
comparison shows whether Ipopt genuinely converges to a good point (it should, with the exact Hessian).

For a state-INDEPENDENT objective (reward/tracking) cfs.evalhessian is the complete exact objective
Hessian. For a state-DEPENDENT objective wrapped as a python function (compliance via
observed_objective.py) it is only the geometric part; the dense state-curvature block is missing and,
given the wrapped function's name via the 'observation' option, added EXACTLY via
observed_objective.state_hessian(n) (n back-substitutions against the factorized state matrix,
self-adjoint functions) - the true Newton Hessian Ipopt's 'exact' mode is designed for. Only if that
is unavailable it is learned by a structured Powell-damped BFGS on the geometric-corrected secant
(see hessian_scipy.py). Auto-enabled iff the named function is state-dependent; without the option
the classic exact-Hessian path is unchanged.

Options (<option key=".." value=".."/> in the optimizer's python element):
  observation: name of the state-dependent cfs function the objective wraps -> enables the exact
               state block (BFGS fallback). Omit for a directly exact objective.
  exact_state: "false" -> skip the exact state block and use the structured BFGS (testing/comparison)
  fd_check   : "true" -> central-difference-check the analytic objective gradient and Hessian at x0.
               With 'observation' the Hessian is quasi-Newton, so only the gradient is asserted.
  fd_eps     : finite difference step (default 1e-6)
  fd_tol     : finite difference tolerance (default 1e-4)
  tol        : Ipopt convergence tolerance (default Ipopt's)
  ipopt.<k>  : passed verbatim to nlp.add_option(<k>, value) with int/float/str inferred, e.g.
               <option key="ipopt.mu_strategy" value="adaptive"/>
"""

import numpy as np
import cfs
import cyipopt # usually not available on CI runners but much better than scipy.optimizer interior-point 


# generate a global object 'glob' where we store e.g. the options from init()
class _G:
  pass
glob = _G()


def setup():
  cfs.opt_register_log_property('grad_fd_err', '-1')
  cfs.opt_register_log_property('hess_fd_err', '-1')


def init(n, m, maxiter, sim_name, options):
  glob.n = n
  glob.m = m
  glob.maxiter = maxiter
  glob.options = options
  glob.sim_name = sim_name  # <sim_name>.ipopt receives Ipopt's own output log, see solve()
  objs, cnstrs = cfs.get_opt_function_values()
  assert len(objs) == 1, "ipopt_hession.py expects a single objective"
  print('ipopt_hession.py: init n=%d m=%d objective=%s constraints=%s'
        % (n, m, list(objs.keys()), list(cnstrs.keys())))

  # static sparse constraint jacobian structure (rows/cols come back as doubles, cast to int)
  glob.nnz = cfs.get_num_jacobian_nonzeros() if m > 0 else 0
  if glob.nnz > 0:
    rows = np.zeros(glob.nnz)
    cols = np.zeros(glob.nnz)
    cfs.get_constraint_sparsity(rows, cols)
    glob.jrows = rows.astype(int)
    glob.jcols = cols.astype(int)

  # we provide a dense lower-triangular Hessian (n is small for feature mapping)
  glob.hrows, glob.hcols = np.tril_indices(n)
  glob.last_x = None

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
    print('ipopt_hession.py: wrapped function %r state-dependent=%s -> state block %s'
          % (obs, glob.structured, 'exact' if glob.exact else ('structured BFGS' if glob.structured else 'none')))
  glob.B = None          # learned state block (BFGS fallback)
  glob.prev = None       # (x, grad) of the last gradient evaluation
  glob.pending = None    # secant candidate (s, y) waiting for the Hessian to correct with H_geo
  glob.collect = False   # gather secant pairs (only during solve, not fd_check)
  glob.scaled = False    # B rescaled from the first accepted pair


def _grad(x):
  g = np.zeros(glob.n)
  cfs.evalgradobj(x, g)
  return g

def _hess(x):
  H = np.zeros((glob.n, glob.n))
  cfs.evalhessian(x, H)
  return H


def _bfgs_update(B, s, y):
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

def _learned_block(Hgeo):
  """update and return the learned state block B from the pending geometric-corrected secant"""
  if glob.B is None:
    glob.B = max(1e-3, abs(np.trace(Hgeo)) / glob.n) * np.eye(glob.n)
  if glob.pending is not None:
    s, y = glob.pending
    glob.pending = None
    y_struct = y - Hgeo @ s
    sy = float(s @ y_struct)
    if not glob.scaled and sy > 0.0:
      glob.B = (float(y_struct @ y_struct) / sy) * np.eye(glob.n)  # first-pair scaling
      glob.scaled = True
    glob.B = _bfgs_update(glob.B, s, y_struct)
  return glob.B


def _fd_check(x):
  eps = float(glob.options.get('fd_eps', 1e-6))
  tol = float(glob.options.get('fd_tol', 1e-4))
  g = _grad(x)
  H = _hess(x)  # geometric part; complete when not structured
  if glob.exact:
    import observed_objective
    H = H + observed_objective.state_hessian(glob.n)  # full exact objective Hessian
  gfd = np.zeros(glob.n)
  Hfd = np.zeros((glob.n, glob.n))
  for j in range(glob.n):
    xp = x.copy(); xp[j] += eps
    xm = x.copy(); xm[j] -= eps
    gfd[j] = (cfs.evalobj(xp) - cfs.evalobj(xm)) / (2 * eps)
    Hfd[:, j] = (_grad(xp) - _grad(xm)) / (2 * eps)
  gerr = float(np.max(np.abs(g - gfd)))
  herr = float(np.max(np.abs(H - Hfd)))
  print('ipopt_hession.py: fd_check max|grad-fd|=%.2e max|hess-fd|=%.2e' % (gerr, herr))
  cfs.opt_set_log_property('grad_fd_err', '%.3e' % gerr)
  cfs.opt_set_log_property('hess_fd_err', '%.3e' % herr)
  assert gerr < tol, 'gradient finite difference error exceeds %.3e' % tol
  if glob.structured and not glob.exact:
    # the objective Hessian is quasi-Newton here; |Hfd-Hgeo| is the missing state block (informative)
    print('ipopt_hession.py: structured - |Hfd-Hgeo|=%.2e is the learned state block, not asserted' % herr)
  else:
    # exact everywhere; the absolute error scales with |H| (e.g. ~1e4 for compliance) -> relative
    hscale = max(1.0, float(np.max(np.abs(Hfd))))
    assert herr / hscale < tol, 'Hessian finite difference error %.3e (rel %.3e) exceeds %.3e' % (herr, herr / hscale, tol)


class _Problem:
  """cyipopt callback object, all derivatives exact from cfs."""

  def objective(self, x):
    glob.last_x = np.array(x)  # for the iterate commit in intermediate()
    return cfs.evalobj(x)

  def gradient(self, x):
    g = _grad(x)
    if glob.collect:
      xa = np.asarray(x)
      if glob.prev is not None and not np.array_equal(xa, glob.prev[0]):
        glob.pending = (xa - glob.prev[0], g - glob.prev[1])
      glob.prev = (xa.copy(), g.copy())
    return g

  def constraints(self, x):
    g = np.zeros(glob.m)
    cfs.evalconstrs(x, g)
    return g

  def jacobianstructure(self):
    return (glob.jrows, glob.jcols)

  def jacobian(self, x):
    vals = np.zeros(glob.nnz)
    cfs.evalgradconstrs_sparse(x, vals)
    return vals

  def hessianstructure(self):
    return (glob.hrows, glob.hcols)

  def hessian(self, x, lagrange, obj_factor):
    # Ipopt wants the lower triangle of the Lagrangian Hessian
    #   obj_factor * grad^2 f + sum_i lagrange_i grad^2 c_i
    Hobj = _hess(x)  # geometric part; complete objective Hessian when not structured
    if glob.structured:
      if glob.exact:
        # exact dense state-curvature block 2 B^T Z (n back-substitutions, see observed_objective)
        import observed_objective
        Hobj = Hobj + observed_objective.state_hessian(glob.n)
      else:
        Hobj = Hobj + _learned_block(Hobj)  # learned dense state-curvature block
    H = obj_factor * Hobj
    if glob.m > 0:
      Hc = np.zeros((glob.n, glob.n))
      cfs.evalhessian_constr(x, np.asarray(lagrange, dtype=float), Hc)
      H += Hc
    return H[glob.hrows, glob.hcols]

  def intermediate(self, *args):
    # commit the current Ipopt iterate as a cfs iteration (evaluate first so the design is synced)
    x = glob.last_x
    try:
      x = np.array(glob.nlp.get_current_iterate()['x'])
    except Exception:
      pass
    if x is not None:
      cfs.evalobj(x)
      cfs.commitIteration()


def solve():
  assert cyipopt is not None, "ipopt_hession.py requires cyipopt (not available on the CI runners)"
  n = glob.n
  m = glob.m

  x0 = np.zeros(n)
  cfs.initialdesign(x0)
  cfs.evalobj(x0)
  cfs.commitIteration()  # iteration 0

  if str(glob.options.get('fd_check', 'false')).lower() == 'true':
    _fd_check(x0)

  glob.collect = glob.structured and not glob.exact  # secant pairs only for the BFGS fallback

  xl = np.zeros(n); xu = np.zeros(n)
  gl = np.zeros(m); gu = np.zeros(m)
  cfs.bounds(xl, xu, gl, gu)  # design bounds and (raw) constraint bounds, +-1e19 for one-sided

  nlp = cyipopt.Problem(n=n, m=m, problem_obj=_Problem(), lb=xl, ub=xu, cl=gl, cu=gu)
  glob.nlp = nlp  # so intermediate() can query get_current_iterate()

  nlp.add_option('max_iter', int(glob.maxiter))
  nlp.add_option('hessian_approximation', 'exact')  # use our exact Lagrangian Hessian
  # keep the original bounds strict: by default Ipopt relaxes them by bound_relax_factor=1e-8, so
  # iterates (which we commit and the FE evaluates) can end up slightly outside, e.g. alpha ~ -1e-8.
  # With 0 the interior-point iterates stay strictly within the bounds. Override via
  # <option key="ipopt.bound_relax_factor" value="..."/> if convergence suffers.
  nlp.add_option('bound_relax_factor', 0.0)
  # write Ipopt's own convergence log to <sim_name>.ipopt (iterations, objective, KKT errors, EXIT)
  ipopt_file = glob.sim_name + '.ipopt'
  nlp.add_option('output_file', ipopt_file)
  nlp.add_option('file_print_level', 5)
  if 'tol' in glob.options:
    nlp.add_option('tol', float(glob.options['tol']))
  # generic pass-through: ipopt.<name> -> add_option(<name>, value) with type inference
  for k, v in glob.options.items():
    if k.startswith('ipopt.'):
      name = k[len('ipopt.'):]
      try:
        val = int(v)
      except ValueError:
        try:
          val = float(v)
        except ValueError:
          val = v
      nlp.add_option(name, val)

  x, info = nlp.solve(x0)
  print('ipopt_hession.py: solve done J=%.6e status=%d %s (log: %s)'
        % (info['obj_val'], info['status'], info['status_msg'], ipopt_file))

  cfs.evalobj(x)
  cfs.commitIteration()  # export Ipopt's final solution (no best-feasible fallback, for honest comparison)
  return 0
