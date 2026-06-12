#ifndef OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_
#define OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_

#include <utility>
#include <tuple>
#include <filesystem>
#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField
{
/** Optimizer for shape- and topology optimization! */
class PythonOptimizer : public BaseOptimizer
{
public:
  PythonOptimizer(Optimization* optimization, PtrParamNode pn);

  virtual ~PythonOptimizer();

  /** the static functions are generally optimization oriented and e.g. available for python functions */

  /** returns 4 values: dim, nx, ny, nz */
  static PyObject* GetDims(PyObject* args);

  /** number of design variables */
  static PyObject* GetNumDesign(PyObject* args);

  /** return single plain design variable */
  static PyObject* GetDesignValue(PyObject* args);

  /** set design variables to provided numpy array of proper size + optional acesss */
  static PyObject* GetDesignValues(PyObject* args);

  /** fill provided numpy array (size = number of pseudo densities, see GetNumDesign) with the
   * element pseudo density DesignSpace::data; for feature mapping the aggregated mrho_e, which
   * differs from the design variables (the shape parameters). Optional second argument access. */
  static PyObject* GetPseudoDensity(PyObject* args);

  /** evaluate transfer function.
   * @param args first argument index, optional second argument is index */
  static PyObject* Transfer(PyObject* args, bool derivative);

  /** expects two 1Dim Arrays for design bounds and two 1Dim arrays for constraint bounds. Not normalized, equal is both the same */
  void GetBounds(PyObject* args);

  void GetInitialDesign(PyObject* args);

  /** Evaluate objective and expects 1Dim design as numpy array */
  double EvalObjective(PyObject* args);

  /** automatically makes this a new iteration, could be changed by function to by called by python */
  void EvalGradObjective(PyObject* args);

  /** evaluate the exact total shape Hessian of the objective in design variable space.
   * Expects a 1D design array x (size n) and a 2D array H (n x n) which is filled.
   * For second-order python optimizers (e.g. scipy trust-constr) and Hessian finite difference tests.
   * Only feature mapping with a curvature-providing objective (see Function::CalcCurvature). */
  void EvalHessian(PyObject* args);

  /** Evaluate constraints and expects 1Dim design as numpy array and 1Dim array of size m for values */
  void EvalConstraints(PyObject* args);

  void EvalGradConstraints(PyObject* args);

  /** Sparse constraint jacobian interface (for many local constraints / scipy.sparse). The structure
   * is static and queried once: get the nnz, then fill the (row, col) index arrays; per iteration fill
   * the packed value array. The python side assembles a scipy.sparse matrix from (vals, (rows, cols)).
   * @see ConditionContainer::VirtualView::CalcNumberOfJacobianNonZeros */
  int GetNumberOfJacobianNonZeros();

  /** fills two 1D arrays of size nnz: rows[k] = constraint index, cols[k] = design variable index */
  void GetConstraintSparsity(PyObject* args);

  /** Expects a 1D design array x (size n) and a 1D value array of size nnz, filled with the packed
   * constraint jacobian values in the same order as GetConstraintSparsity */
  void EvalGradConstraintsSparse(PyObject* args);

  /** the multiplier-contracted constraint Hessian sum_c v_c * Hess(c_c)(x) as a dense n x n matrix,
   * for a second-order optimizer's Lagrangian Hessian (e.g. scipy NonlinearConstraint hess=). Expects
   * a 1D design x (size n), a 1D multiplier array v (size m) and a 2D n x n array which is filled.
   * Uses the raw (unnormalized) per-constraint Hessian, see Condition::CalcHessian. */
  void EvalConstraintHessian(PyObject* args);

  /** this should be in PythonKernel but somehow this gives a segfault which indicates that PyArg_ParseTuple is not defined?!
   * @param decref if false make sure to decref the objects via the return array
   * @return you must not use the PyObjects when decref is true */
  static StdVector<PyObject*> ParseArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref);

  double GetSimpExponent();

  /** Returns derivative of compliance w.r.t. material tensor */
  void Get_dfdH(PyObject *args);

  void GetCoreStiffness(PyObject *args);

protected:

  /** Implements virtual function from BaseOptimizer */
  void SolveProblem();
private:

  /** the options from the xml file */
  StdVector<std::pair<std::string, std::string> > options;

  /** given filename. Possibly via file and path which path = cfs:share:python key */
  std::string givenname;

  /** keeps the python script */
  PyObject* module = NULL;

  /** number of design variables */
  unsigned int n = 0;

  /** number of constraints */
  int m = -1;

  /** our python info instead of optimizer info_ */
  PtrParamNode pyinf_ = NULL;

};


} // end of namespace

#endif /* OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_ */
