#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include <assert.h>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/PythonTools.hh"
#include "MatVec/Vector.hh"
#include "Utils/tools.hh"


// declare class specific logging stream
DEFINE_LOG(pyopt, "pyopt")


namespace CoupledField
{
using std::string;
using std::to_string;
using std::make_pair;

/** global pointer to the snopt class to be used by the callback function */
PythonOptimizer* static_pyopt = NULL;

/** cfs.bound(xl,xu,gl,gu) sets bounds for design and constraints in properly sized 1d numpy arrays */
static PyObject* cfs_bounds(PyObject *self, PyObject *args)
{
  static_pyopt->GetBounds(args);
  Py_RETURN_NONE;
}

/** cfs.initialdesign(x) fills 1d numpy array */
static PyObject* cfs_initialdesign(PyObject *self, PyObject *args)
{
  static_pyopt->GetInitialDesign(args);
  Py_RETURN_NONE;
}

/** cfs.evalobj(x) returns a float) */
static PyObject* cfs_evalobj(PyObject *self, PyObject *args)
{
  return PyFloat_FromDouble(static_pyopt->EvalObjective(args));
}

/** cfs.evalgradobj(x,grad) fills numpy 1d arrays of size n */
static PyObject* cfs_evalgradobj(PyObject *self, PyObject *args)
{
  static_pyopt->EvalGradObjective(args);

  Py_RETURN_NONE;
}

/** cfs.cfs_evalconstrs(x,g) fills numpy 1d arrays of size n and m */
static PyObject* cfs_evalconstrs(PyObject *self, PyObject *args)
{
  static_pyopt->EvalConstraints(args);

  Py_RETURN_NONE;
}

/** cfs.cfs_evalgradconstrs(x,grad) fills numpy 1d arrays of size n and m*n */
static PyObject* cfs_evalgradconstrs(PyObject *self, PyObject *args)
{
  static_pyopt->EvalGradConstraints(args);

  Py_RETURN_NONE;
}

/** return true if cfs's stopping criteria is met, including finding the file HALTOPT */
static PyObject* cfs_dostop(PyObject *self, PyObject *args)
{
  bool stop = static_pyopt->optimization->DoStopOptimization();

  if(stop)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

static PyMethodDef cfs_methods[] = {
    {"bounds", cfs_bounds, METH_VARARGS, "Give design and constraints bounds. Expects 1D arrays for xl, xu, gl, gu"},
    {"initialdesign", cfs_initialdesign, METH_VARARGS, "Sets the initial design to the given 1D array"},
    {"evalobj", cfs_evalobj, METH_VARARGS, "Evaluate objective. Expects 1D design array"},
    {"evalconstrs", cfs_evalconstrs, METH_VARARGS, "Evaluate constraints. Expects 1D design array and 1D value array or size m"},
    {"evalgradobj", cfs_evalgradobj, METH_VARARGS, "Evaluate objective gradient. Expects 1D design array and dense gradient array "},
    {"evalgradconstrs", cfs_evalgradconstrs, METH_VARARGS, "Evaluate constrains gradients. Expects 1D design array of size n and dense 1D value array of size m*n"},
    {"dostop", cfs_dostop, METH_VARARGS, "Uses stopping criteria from cfs (including HALTOPT file) and returns True or False"},
    {NULL, NULL, 0, NULL}
};

static PyModuleDef cfs_modules = {
    PyModuleDef_HEAD_INIT, "cfs", NULL, -1, cfs_methods, NULL, NULL, NULL, NULL
};

static PyObject* PyInit_cfs(void)
{
  // https://stackoverflow.com/questions/37943699/crash-when-calling-pyarg-parsetuple-on-a-numpy-array
  import_array();

  return PyModule_Create(&cfs_modules);
}

PythonOptimizer::PythonOptimizer(Optimization* opt, PtrParamNode pn) :
  BaseOptimizer(opt, pn, Optimization::PYTHON_SOLVER)
{
  static_pyopt = this;
  pyinf_ = info_->Get(Optimization::optimizer.ToString(Optimization::PYTHON_SOLVER));

  BaseOptimizer::PostInitScale(1.0);

  n = optimization->GetDesign()->GetNumberOfVariables();
  m =  optimization->constraints.view->GetNumberOfActiveConstraints();

  assert(this_opt_pn_ != NULL);

  givenname = this_opt_pn_->Get("file")->As<string>();

  string version;
  module = InitializePythonModule(givenname.string(), PyInit_cfs, &version);

  pyinf_->Get(ParamNode::HEADER)->Get("givenname")->SetValue(givenname.string());
  pyinf_->Get(ParamNode::HEADER)->Get("python")->SetValue(version);

  // the options are given to the python functions setup() and solve() for their usage
  ParamNodeList lst = this_opt_pn_->GetList("option");
  options = ParseOptions(lst);

  pyinf_->Get(ParamNode::HEADER)->Get("options")->SetValue(lst);

  // call optional setup(n,m,options)
  std::string val;
  PyObject* setup = PyObject_GetAttrString(module, "setup");
  if(setup && PyCallable_Check(setup)) {
    PyObject* arg = PyTuple_New(3);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(n));
    PyTuple_SetItem(arg, 1, PyLong_FromLong(m));
    PyTuple_SetItem(arg, 2, CreatePythonDict(options));
    PyObject* ret = PyObject_CallObject(setup, arg);
    val = ret != NULL ? to_string(PyLong_AsLong(ret)) : "NULL";
    Py_XDECREF(ret);
    Py_XDECREF(arg);
    Py_XDECREF(setup);
  }
  else
    val = "not callable";
  pyinf_->Get(ParamNode::HEADER)->Get("setup()")->SetValue(val);
}

PythonOptimizer::~PythonOptimizer()
{
  Py_XDECREF(module);
  Py_Finalize();
}


void PythonOptimizer::SolveProblem()
{
  assert(module);

  // call mandatory solve(n,m,maxiter, options)
  std::string val;
  PyObject* solve = PyObject_GetAttrString(module, "solve");
  if(!solve || !PyCallable_Check(solve))
    throw Exception("no solve() in python module " + givenname.string());

  PyObject* arg = PyTuple_New(4);
  PyTuple_SetItem(arg, 0, PyLong_FromLong(n));
  PyTuple_SetItem(arg, 1, PyLong_FromLong(m));
  PyTuple_SetItem(arg, 2, PyLong_FromLong(optimization->GetMaxIterations()));
  PyTuple_SetItem(arg, 3, CreatePythonDict(options));
  PyObject* ret = PyObject_CallObject(solve, arg);
  CheckPythonReturn(ret);
  pyinf_->Get(ParamNode::SUMMARY)->Get("solve()")->SetValue(PyLong_AsLong(ret));
  Py_DECREF(ret);
  Py_XDECREF(arg);
  Py_DECREF(solve);
}


StdVector<PyObject*> PythonOptimizer::ParseArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref)
{
  if(PyTuple_GET_SIZE(args) != expect)
    throw Exception("expects " + to_string(expect) + " parameters in python function call, got " + to_string(PyTuple_GET_SIZE(args)));

  StdVector<PyObject*> obj(expect);

  // could clearly be done smarter!
  switch(expect)
  {
  case 1:
    PyArg_ParseTuple(args, "O!", &PyArray_Type, &obj[0]);
    break;

  case 2:
    PyArg_ParseTuple(args, "O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1]);
    break;

  case 3:
    PyArg_ParseTuple(args, "O!|O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1], &PyArray_Type, &obj[2]);
    break;

  case 4:
    PyArg_ParseTuple(args, "O!|O!|O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1], &PyArray_Type, &obj[2],  &PyArray_Type, &obj[3]);

    break;
  default:
    EXCEPTION("not implemented size " << expect);
  }
  data.Resize(expect);

  for(unsigned int i = 0; i < obj.GetSize(); i++)
    data[i].Fill(obj[i], decref);

  return obj;
}

void PythonOptimizer::GetInitialDesign(PyObject* args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 1, data, false);

  Vector<double>& x = data[0];

  optimization->GetDesign()->WriteDesignToExtern(x);

  assert(obj.GetSize() == 1);
  x.Export(obj[0]);
}

void PythonOptimizer::GetBounds(PyObject *args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 4, data, false);

  Vector<double>& xl = data[0];
  Vector<double>& xu = data[1];
  Vector<double>& gl = data[2];
  Vector<double>& gu = data[3];

  assert(xl.GetSize() == xu.GetSize());
  assert(gl.GetSize() == gu.GetSize());

  BaseOptimizer::GetBounds(xl.GetSize(), xl.GetPointer(), xu.GetPointer(), gl.GetSize(), gl.GetPointer(), gu.GetPointer());

  assert(obj.GetSize() == data.GetSize());
  for(unsigned int i = 0; i < obj.GetSize(); i++)
    data[i].Export(obj[i]); // we must to decref obj!
}

double PythonOptimizer::EvalObjective(PyObject *args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 1, data, false);

  const Vector<double>& x = data[0];

  return BaseOptimizer::EvalObjective(x.GetSize(), x.GetPointer(), true);
}

void PythonOptimizer::EvalGradObjective(PyObject *args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 2, data, false);

  const Vector<double>& x = data[0];
  Vector<double>& grad = data[0];

  assert(x.GetSize() == n && grad.GetSize() == n);

  // we use stdgrad just as a wrapper for grad as EvalGrad*() uses StdVector because of the Window.
  StdVector<double> stdgrad;
  stdgrad.Assign(grad.GetPointer(), grad.GetSize(), true); // copies data pointer from grad.

  BaseOptimizer::EvalGradObjective(x.GetSize(), x.GetPointer(), true, stdgrad); // grad also gets the new data

  assert(close(grad.Avg(),Average(stdgrad.GetPointer(), stdgrad.GetSize()), 1e-10));

  stdgrad.Assign(NULL, 0, false); // make stdgrad forget about data in grad to not delete it with stdgrad destrutor

  grad.Export(obj[1]);

  // could be moved to an (optional) python callback
  optimization->CommitIteration(); // return paramnode ignored
}


void PythonOptimizer::EvalConstraints(PyObject *args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 2, data, false);

  const Vector<double>& x = data[0];
  Vector<double>& gval = data[1];

  // one could add "normalize" as option
  BaseOptimizer::EvalConstraints(x.GetSize(), x.GetPointer(), gval.GetSize(), true, gval.GetPointer(), false); // scale=true, normalize=false

  gval.Export(obj[1]);
}

void PythonOptimizer::EvalGradConstraints(PyObject *args)
{
  // see EvalGradObjective() for comments

  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 2, data, false);

  const Vector<double>& x = data[0];
  Vector<double>& grad = data[0];

  assert(x.GetSize() == n && grad.GetSize() == m*n);

  StdVector<double> stdgrad;
  stdgrad.Assign(grad.GetPointer(), grad.GetSize(), true); // copies data pointer from grad

  BaseOptimizer::EvalGradConstraints(x.GetSize(), x.GetPointer(), m, m*n, true, false, stdgrad); // scale=true, normalize=false

  stdgrad.Assign(NULL, 0, false);

  grad.Export(obj[1]);
}




} // namespace
