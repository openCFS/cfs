"""Drop-in cyipopt replacement for hessian.py (same PythonOptimizer interface and cfs.* calls).

Where hessian.py drives scipy trust-constr, this drives Ipopt via cyipopt with the *exact* Lagrangian
Hessian: the objective Hessian cfs.evalhessian plus the multiplier-contracted constraint Hessian
cfs.evalhessian_constr. The native 'distance' constraints come from cfs (sparse Jacobian + exact
Hessian), exactly as in hessian.py.

This is a developer/comparison driver: cyipopt is NOT available on the CI runners, so the testsuite
keeps hessian.py. Swap <python file="hessian.py"> -> <python file="ipopt_hession.py"> to compare. Unlike
hessian.py it does NOT keep the best feasible iterate - it exports Ipopt's actual final solution, so the
comparison shows whether Ipopt genuinely converges to a good point (it should, with the exact Hessian).

Options (<option key=".." value=".."/> in the optimizer's python element):
  fd_check   : "true" -> central-difference-check the analytic objective gradient and Hessian at x0.
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


def _grad(x):
  g = np.zeros(glob.n)
  cfs.evalgradobj(x, g)
  return g

def _hess(x):
  H = np.zeros((glob.n, glob.n))
  cfs.evalhessian(x, H)
  return H


def _fd_check(x):
  eps = float(glob.options.get('fd_eps', 1e-6))
  tol = float(glob.options.get('fd_tol', 1e-4))
  g = _grad(x)
  H = _hess(x)
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
  assert gerr < tol and herr < tol, 'finite difference error exceeds %.3e' % tol


class _Problem:
  """cyipopt callback object, all derivatives exact from cfs."""

  def objective(self, x):
    glob.last_x = np.array(x)  # for the iterate commit in intermediate()
    return cfs.evalobj(x)

  def gradient(self, x):
    return _grad(x)

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
    H = obj_factor * _hess(x)
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

  xl = np.zeros(n); xu = np.zeros(n)
  gl = np.zeros(m); gu = np.zeros(m)
  cfs.bounds(xl, xu, gl, gu)  # design bounds and (raw) constraint bounds, +-1e19 for one-sided

  nlp = cyipopt.Problem(n=n, m=m, problem_obj=_Problem(), lb=xl, ub=xu, cl=gl, cu=gu)
  glob.nlp = nlp  # so intermediate() can query get_current_iterate()

  nlp.add_option('max_iter', int(glob.maxiter))
  nlp.add_option('hessian_approximation', 'exact')  # use our exact Lagrangian Hessian
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
