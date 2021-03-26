#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include <assert.h>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"
#include "Optimization/Design/DesignElement.hh"
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
static PyObject* cfs_getDims(PyObject *self, PyObject *args)
{
  return static_pyopt->GetDims(args);
}

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

/** cfs.cfs_commitIteration() commits iteration to cfs */
static PyObject* cfs_commitIteration(PyObject *self, PyObject *args)
{
  static_pyopt->optimization->CommitIteration(); // return paramnode ignored;

  Py_RETURN_NONE;
}

/** cfs.evalgradobj(x,grad) fills numpy 1d arrays of size n */
static PyObject* cfs_evalgradobj(PyObject *self, PyObject *args)
{
  static_pyopt->EvalGradObjective(args);
//  cfs_commitIteration(self,args);

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

static PyObject* cfs_getSimpExponent(PyObject *self, PyObject *args)
{
  return PyFloat_FromDouble(static_pyopt->GetSimpExponent());
}

/** returns derivative of compliance w.r.t stiffness tensor entries of original (core) material */
static PyObject* cfs_get_dfdH(PyObject *self, PyObject *args)
{
  static_pyopt->Get_dfdH(args);

  Py_RETURN_NONE;
}

/** cfs.cfs_getOrgStiffness(stiffness) returns stiffness tensor of original (core) material */
static PyObject* cfs_getOrgStiffness(PyObject *self, PyObject *args)
{
  static_pyopt->GetCoreStiffness(args);

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
    {"getDims", cfs_getDims, METH_VARARGS, "Returns info on mesh dimensions: dim, nx, ny, nz."},
    {"bounds", cfs_bounds, METH_VARARGS, "Give design and constraints bounds. Expects 1D arrays for xl, xu, gl, gu"},
    {"initialdesign", cfs_initialdesign, METH_VARARGS, "Sets the initial design to the given 1D array"},
    {"evalobj", cfs_evalobj, METH_VARARGS, "Evaluate objective. Expects 1D design array"},
    {"evalconstrs", cfs_evalconstrs, METH_VARARGS, "Evaluate constraints. Expects 1D design array and 1D value array or size m"},
    {"evalgradobj", cfs_evalgradobj, METH_VARARGS, "Evaluate objective gradient w. r. t design variables. Expects 1D design array and dense gradient array "},
    {"evalgradconstrs", cfs_evalgradconstrs, METH_VARARGS, "Evaluate constrains gradients. Expects 1D design array of size n and dense 1D value array of size m*n"},
    {"getSimpExponent", cfs_getSimpExponent, METH_VARARGS, "Returns power-law exponent."},
    {"get_dfdH", cfs_get_dfdH, METH_VARARGS, "Returns derivative of objective w.r.t. material tensor. Expects 3D array of size nelems x 3 x 3"},
    {"getOrgStiffness", cfs_getOrgStiffness, METH_VARARGS, "Returns stiffness tensor of original (core) material."},
    {"commitIteration", cfs_commitIteration, METH_VARARGS, "Commit current iteration to cfs."},
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

  string version;
  module = InitializePythonModule(this_opt_pn_->Get("file")->As<string>(), this_opt_pn_->Get("path")->As<string>(), PyInit_cfs, &givenname,&version);

  pyinf_->Get(ParamNode::HEADER)->Get("file")->SetValue(givenname);
  pyinf_->Get(ParamNode::HEADER)->Get("version")->SetValue(version);

  // the options are given to the python functions setup() and init() for their usage
  ParamNodeList lst = this_opt_pn_->GetList("option");
  options = ParseOptions(lst);

  pyinf_->Get(ParamNode::HEADER)->Get("options")->SetValue(lst);

  // call optional setup()
  std::string val;
  PyObject* setup = PyObject_GetAttrString(module, "setup");
  if(setup && PyCallable_Check(setup)) {
    PyObject* ret = PyObject_CallObject(setup, NULL);
    val = ret != NULL ? to_string(PyLong_AsLong(ret)) : "NULL";
    Py_XDECREF(ret);
    Py_XDECREF(setup);
    if (val == "NULL")
      WARN("Failed to execute setup() in python!");
  }
  else
    val = "not callable";
  pyinf_->Get(ParamNode::HEADER)->Get("setup()")->SetValue(val);

  // call init(n,m,iters,options)
  PyObject* init = PyObject_GetAttrString(module, "init");
  if(init && PyCallable_Check(init)) {
    PyObject* arg = PyTuple_New(5);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(n));
    PyTuple_SetItem(arg, 1, PyLong_FromLong(m));
    PyTuple_SetItem(arg, 2, PyLong_FromLong(optimization->GetMaxIterations()));
    PyTuple_SetItem(arg, 3, PyUnicode_FromString(progOpts->GetSimName().c_str()));
    PyTuple_SetItem(arg, 4, CreatePythonDict(options));
    PyObject* ret = PyObject_CallObject(init, arg);
    val = ret != NULL ? to_string(PyLong_AsLong(ret)) : "NULL";
    Py_XDECREF(ret);
    Py_XDECREF(arg);
    Py_XDECREF(init);
    if (val == "NULL")
      throw Exception("Failed to execute init() in python!");
  }
  else
    val = "not callable";

  pyinf_->Get(ParamNode::HEADER)->Get("init()")->SetValue(val);
  optimizer_timer_->Stop();
}

PythonOptimizer::~PythonOptimizer()
{
  Py_XDECREF(module);
  Py_Finalize();
}


void PythonOptimizer::SolveProblem()
{
  assert(module);

  // call mandatory solve()
  std::string val;
  PyObject* solve = PyObject_GetAttrString(module, "solve");
  if(!solve || !PyCallable_Check(solve))
    throw Exception("no solve() in python module " + givenname);

  PyObject* ret = PyObject_CallObject(solve, NULL);
  if(!ret) {
    PyErr_Print();
    throw Exception("error calling solve() in python module " + givenname);
  }
  pyinf_->Get(ParamNode::SUMMARY)->Get("solve()")->SetValue(PyLong_AsLong(ret));
  Py_DECREF(ret);
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
  for(unsigned int i = 0; i < obj.GetSize(); i++) {
    // can only perform this when parsed np arrays are 1D
    if (PyArray_NDIM((PyArrayObject*)obj[i]) == 1)
      data[i].Fill(obj[i], decref);
  }

  return obj;
}

void PythonOptimizer::GetInitialDesign(PyObject* args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 1, data, false);

  Vector<double>& x = data[0];

  optimization->GetDesign()->WriteDesignToExtern(x);

  LOG_DBG3(pyopt) << "Initial design: " << x.ToString();

  assert(obj.GetSize() == 1);
  x.Export(obj[0]);
}

PyObject* PythonOptimizer::GetDims(PyObject* args)
{
  unsigned int dim = domain->GetGrid()->GetDim();
  // assume rectilinear mesh with one region only
  assert(optimization->GetDesign()->GetRegionIds().GetSize() == 1);
  StdVector<unsigned int> bounds = domain->GetGrid()->GetBoundaries(optimization->GetDesign()->GetRegionIds().First());
  assert(bounds.GetSize() == 3);

  PyObject* res = PyTuple_New(4);
  PyTuple_SetItem(res, 0, PyLong_FromLong(dim));
  PyTuple_SetItem(res, 1, PyLong_FromLong(bounds[0]));
  PyTuple_SetItem(res, 2, PyLong_FromLong(bounds[1]));
  PyTuple_SetItem(res, 3, PyLong_FromLong(bounds[2]));

  return res;
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

//  optimization->GetDesign()->DisableTransferFunctions();

  double val = BaseOptimizer::EvalObjective(x.GetSize(), x.GetPointer(), true);
//  optimization->GetDesign()->EnableTransferFunctions();
  return val;
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

  //std::cout << "grad.Avg():" << grad.Avg() << " Average(stdgrad.GetPointer(), stdgrad.GetSize()):" << Average(stdgrad.GetPointer(), stdgrad.GetSize()) << std::endl;
  //assert(close(grad.Avg(),Average(stdgrad.GetPointer(), stdgrad.GetSize()), 1e-10));

  stdgrad.Assign(NULL, 0, false); // make stdgrad forget about data in grad to not delete it with stdgrad destructor

  LOG_DBG3(pyopt) << "objective gradient:" << grad.ToString() << std::endl;

  grad.Export(obj[1]);
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

inline double PythonOptimizer::GetSimpExponent() {
  return optimization->GetDesign()->GetTransferFunction(DesignElement::DENSITY, App::MECH)->GetParam();
}

void PythonOptimizer::Get_dfdH(PyObject *args)
{
  // we get numpy 3d array from python and interprete it as a list of 2d numpy matrices
  PyArrayObject *list;
  PyArg_ParseTuple(args, "O!", &PyArray_Type, &list);
  if(!list)
    PyErr_Print();

  // number of matrices inside list
  unsigned int n_mat = PyArray_DIM(list,0);
  // assume all matrices have same dimensions
  assert(PyArray_DIM(list,1) == PyArray_DIM(list,2));

  LOG_DBG3(pyopt) << "GTD: Got " << n_mat << " matrices with dim=(" << PyArray_DIM(list,1) << "," << PyArray_DIM(list,2) << ") from python" << std::endl;

  unsigned int rows = (domain->GetGrid()->GetDim() == 3) ? 6 : 3;
  // get objective function
  assert(optimization->objectives.Has(Function::COMPLIANCE));
  Function* f =  optimization->objectives.Get(Function::COMPLIANCE);

  assert(!f->elements.IsEmpty());

  // for element e: df/dD = -B_e * u * u^T * B_e^T
  // need to find Bmat for all elements: basically a copy of MagSIMP:CalcMagFluxDensity, see also BDBInt::CalcElementMatrix
  // get the form first
  Context* context = Optimization::context;

  DesignSpace* space = optimization->GetDesign();
  // GetRegionId() is only valid if we have 1 region
  assert(space->GetRegionIds().GetSize() == 1);
  BaseBDBInt* bdb = dynamic_cast<BaseBDBInt*>(context->pde->GetAssemble()->GetBiLinForm("LinElastInt", space->GetRegionId(), context->pde)->GetIntegrator());
  assert(bdb != NULL);

  // annoying entity iterator got hold the elem
  ElemList el(domain->GetGrid());

  ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
  StateContainer& forward = em->GetForwardStates();

  unsigned int nelems = space->GetNumberOfElements();
  assert(nelems == n_mat);
  LOG_DBG3(pyopt) << "GdfDH: f=" << f->ToString() << " nelems in design space=" << nelems;
  // we have n_elems B matrices
  StdVector<Matrix<double>> res;
  res.Resize(nelems);
  for (unsigned int e = 0; e < nelems; e++) {
    res[e].Resize(rows,rows);
    res[e].Init();
  }

  StdVector<LocPoint> intPoints; // Get integration Points
  LocPointMapped lp;
  StdVector<double> weights;
  IntegOrder order;
  IntScheme::IntegMethod method = IntScheme::UNDEFINED;

  // do we have multiple excitations?
  MultipleExcitation* me = optimization->me;
  bool mult_excite = me->IsEnabled();

  // the multiple excitation case is a special case - for all other cases this is executed once
  for(unsigned int e = 0; e < (mult_excite ? me->excitations.GetSize() : 1); e++)
  {
    Excitation* excite = mult_excite ? &(me->excitations[e]) : f->GetExcitation();
    assert(excite != NULL);
    assert(excite->index < (int) me->excitations.GetSize());
    // the stored element solution vector
    StdVector<SingleVector*>& sol = forward.Get(excite)->elem[App::MECH];

    for(unsigned int e = 0; e < nelems; e++)
    {
      // elemNum is 1-based
      const Elem* elem = domain->GetGrid()->GetElem(e+1);

      Vector<double>& u_e = dynamic_cast<Vector<double>&>(*sol[e]);
      // precompute u*u^T - outer product gives 8 x 8 matrix
      // need to convert u to a matrix and let matrix class do the multiplication
      Matrix<double> umat;
      umat.Assign(u_e, u_e.GetSize(), 1, true);
      // transposed u_e
      Matrix<double> umatT;
      umat.Transpose(umatT);

      umat *= umatT;

      LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " u_e=" << u_e.ToString(2);
      LOG_DBG3(pyopt) << " mat u_e^T=" << umatT.ToString(2,true);
      LOG_DBG3(pyopt) << " u_e*u_e^T=" << umat.ToString(2,true);

      // prepare to get the curl operator
      el.SetElement(elem);

      BaseFE* ptFe = bdb->GetFeSpace1()->GetFe(el.GetIterator(), method, order );

      bdb->GetIntScheme()->GetIntPoints(Elem::GetShapeType(elem->type), method, order, intPoints, weights );
      assert(method != IntScheme::UNDEFINED);
      assert(!intPoints.IsEmpty());
      // Get shape map from grid
      shared_ptr<ElemShapeMap> esm = domain->GetGrid()->GetElemShapeMap(elem);
      // for intermediate steps
      Matrix<double> bMatT;
      // this is the result
      Matrix<double> bMat;
      for(unsigned int ip = 0; ip < intPoints.GetSize(); ip++)
      {
        // Calculate for each integration point the LocPointMapped
        lp.Set(intPoints[ip], esm, weights[ip]);

        // Call the CalcBMat()-method
        assert(bdb->GetBOp());
        bdb->GetBOp()->CalcOpMat(bMat, lp, ptFe);
        bMat.Transpose(bMatT);

        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " bMat=" << bMat.ToString(2,true);
        // multiplications for derivative of compliance w. r. t. material tensor: B*u*u^T*B^T
        Matrix<double> uuTBT;
        uuTBT.Resize(u_e.GetSize(),rows);
        // u*u^T*B^T
        umat.Mult(bMatT,uuTBT);
        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " u*u^T*B^T=" << uuTBT.ToString(2, true);
        // B*u*u^T*B^T
        bMat *= uuTBT;

        res[e] -= bMat * weights[ip] * lp.jacDet;

  //      std::cout << "elem" <<  elem->elemNum << " ip=" << ip << " weight=" << weights[ip] << " jac det=" << lp.jacDet << std::endl;

        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " B*u*u^T*B^T=" << bMat.ToString(2, true);
        LOG_DBG3(pyopt) << " tmp=" << res[e].ToString(2,true);
      } // integration points
      LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " dfdH=" << res[e].ToString(2,true);
    } // elements
  } // excitation

  for (unsigned int e = 0; e < n_mat; e++) {
    Matrix<double>& dfdH_e = res[e];
    assert(dfdH_e.GetNumRows() == rows);
    assert(dfdH_e.GetNumCols() == rows);
    for (unsigned int i = 0; i < rows; i++)
      for (unsigned int j = 0; j < rows; j++){
        (*((double*) PyArray_GETPTR3(list,e,i,j))) = dfdH_e[i][j];
      }
    // elemNum is 1-based
    const Elem* elem = domain->GetGrid()->GetElem(e+1);
    LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " dfdH=" << dfdH_e.ToString(2,true);
  }
}

void PythonOptimizer::GetCoreStiffness(PyObject *args)
{
  StdVector<Vector<double> > data;
  StdVector<PyObject*> obj = ParseArrays(args, 1, data, false);

  PyArrayObject* arr = (PyArrayObject*) obj[0];
  // assume all matrices have same dimensions
  unsigned int rows = PyArray_DIM(arr,0);
  unsigned int cols = PyArray_DIM(arr,1);
  LOG_DBG2(pyopt) << "numpy array rows=" <<  rows << " cols=" << cols;
  assert(rows == cols);
  assert(rows == 6 || rows == 3 );

  MaterialTensor<double> tens(VOIGT);
  DesignMaterial* dm = Optimization::context->dm;
  assert(dm != NULL);
  // take first element in domain
  const Elem* elem = domain->GetGrid()->GetElem(1);

  // we want pure stiffness tensor without density rho etc -> provide external identity transfer function to design material
  dm->GetTensor(tens, DesignElement::ALL_DESIGNS, Optimization::context->stt, elem, DesignElement::NO_DERIVATIVE, VOIGT, true);
  Matrix<double>& material = tens.GetMatrix(VOIGT);

  // fill numpy array
  for (unsigned int r = 0; r < rows; r++)
    for (unsigned int c = 0; c < cols; c++) {
      *((double*) PyArray_GETPTR2(arr,r,c) ) = material[r][c] ;
    }

}

} // namespace
