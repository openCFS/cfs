#define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
//#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <assert.h>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "Driver/Assemble.hh"
#include "Forms/BiLinForms/BDBInt.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/TransferFunction.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "MatVec/Vector.hh"
#include "Utils/PythonKernel.hh"
#include "Utils/tools.hh"

// declare class specific logging stream
DEFINE_LOG(pyopt, "pyopt")


namespace CoupledField
{
using std::string;
using std::to_string;
using std::make_pair;

PythonOptimizer::PythonOptimizer(Optimization* opt, PtrParamNode pn) :
  BaseOptimizer(opt, pn, Optimization::PYTHON_SOLVER)
{
  pyinf_ = info_->Get(Optimization::optimizer.ToString(Optimization::PYTHON_SOLVER));
  PtrParamNode pnh = pyinf_->Get(ParamNode::HEADER);

  BaseOptimizer::PostInitScale(1.0);

  n = optimization->GetDesign()->GetNumberOfVariables();
  m =  optimization->constraints.view->GetNumberOfActiveConstraints();

  assert(this_opt_pn_ != NULL);

  python->Register(this);

  PythonKernel::LoadStatus stat = python->LoadPythonModule(this_opt_pn_);
  module = stat.module;
  givenname = stat.full_file;

  // this stores the python file in the h5 output file
  std::ifstream fstream(this_opt_pn_->Get("file")->As<string>());
  std::stringstream buffer;
  buffer << fstream.rdbuf();
  domain->GetSimState()->GetOutputWriter()->WriteStringToUserData("PythonFile", buffer.str());

  pnh->Get("file")->SetValue(givenname);
  if(progOpts->DoDetailedInfo())
    pnh->Get("syspath")->SetValue(stat.sys_path.ToString(TS_PLAIN, ":"));

  // the options are given to the python functions setup() and init() for their usage
  ParamNodeList lst = this_opt_pn_->GetList("option");
  options = PythonKernel::ParseOptions(lst);
  pnh->Get("options")->SetValue(lst);

  // https://stackoverflow.com/questions/37943699/crash-when-calling-pyarg-parsetuple-on-a-numpy-array
  import_array1(); // having this in the first call magically prevents segfaults with PyArg_ParseTuple

  // call optional setup()
  std::string val = "not callable";
  PyObject* setup = PyObject_GetAttrString(module, "setup");
  if(setup && PyCallable_Check(setup)) {
    PyObject* ret = PyObject_CallObject(setup, NULL);
    PythonKernel::CheckPythonReturn(ret, "setup"); // setup is optional, nut if there is function it shall work!
    val = PythonKernel::ToString(ret);
    Py_XDECREF(ret);
    Py_XDECREF(setup);
  }
  pnh->Get("setup")->SetValue(val);

  // call mandatory init(n ,m, iters, name, options)
  PyObject* init = PyObject_GetAttrString(module, "init");
  if(init && PyCallable_Check(init)) {
    PyObject* arg = PyTuple_New(5);
    PyTuple_SetItem(arg, 0, PyLong_FromLong(n));
    PyTuple_SetItem(arg, 1, PyLong_FromLong(m));
    PyTuple_SetItem(arg, 2, PyLong_FromLong(optimization->GetMaxIterations()));
    PyTuple_SetItem(arg, 3, PyUnicode_FromString(progOpts->GetSimName().c_str()));
    PyTuple_SetItem(arg, 4, PythonKernel::CreatePythonDict(options));
    PyObject* ret = PyObject_CallObject(init, arg);
    PythonKernel::CheckPythonReturn(ret, "init");
    pnh->Get("init")->SetValue(PythonKernel::ToString(ret));
    Py_XDECREF(ret);
    Py_XDECREF(arg);
    Py_XDECREF(init);
  }
  optimizer_timer_->Stop();
}

PythonOptimizer::~PythonOptimizer()
{
  python->Register(this, false);
  Py_XDECREF(module);
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
  if(!ret)
    throw Exception("error calling solve() in python module " + givenname + ": " + PythonKernel::PyErr());

  pyinf_->Get(ParamNode::SUMMARY)->Get("solve()")->SetValue(PyLong_AsLong(ret));
  Py_DECREF(ret);
  Py_DECREF(solve);
}


StdVector<PyObject*> PythonOptimizer::ParseArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref)
{
  if(PyTuple_Size(args) != expect)
    throw Exception("expects " + to_string(expect) + " parameters in python function call, got " + to_string(PyTuple_Size(args)));

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
  Vector<double>& grad = data[1];

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
  Vector<double>& grad = data[1];

  LOG_DBG(pyopt) << "EGC: #x=" << x.GetSize() << " #g=" << grad.GetSize() << " n=" << n << " m=" << m;

  assert(x.GetSize() == n && grad.GetSize() == m*n);

  StdVector<double> stdgrad;
  stdgrad.Assign(grad.GetPointer(), grad.GetSize(), true); // copies data pointer from grad

  BaseOptimizer::EvalGradConstraints(x.GetSize(), x.GetPointer(), m, m*n, true, false, stdgrad); // scale=true, normalize=false

  stdgrad.Assign(NULL, 0, false);

  grad.Export(obj[1]);
}

double PythonOptimizer::GetSimpExponent() {
  return optimization->GetDesign()->GetTransferFunction(DesignElement::DENSITY, App::MECH)->GetParam();
}

void PythonOptimizer::Get_dfdH(PyObject *args)
{
  // we get numpy 3d array from python and interprete it as a list of 2d numpy matrices
  PyArrayObject *list;
  PyArg_ParseTuple(args, "O!", &PyArray_Type, &list);
  if(!list)
    throw Exception("no list given to Get_dfdH: " + PythonKernel::PyErr());

  // number of matrices inside list
  unsigned int n_mat = PyArray_DIM(list,0);
  // assume all matrices have same dimensions
  assert(PyArray_DIM(list,1) == PyArray_DIM(list,2));

  LOG_DBG3(pyopt) << "GTD: Got " << n_mat << " matrices with dim=(" << PyArray_DIM(list,1) << "," << PyArray_DIM(list,2) << ") from python" << std::endl;

  unsigned int rows = (domain->GetDim() == 3) ? 6 : 3;
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
  for(unsigned int ex = 0; ex < (mult_excite ? me->excitations.GetSize() : 1); ex++)
  {
    Excitation* excite = mult_excite ? &(me->excitations[ex]) : f->GetExcitation();
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

      LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " u_e=" << u_e.ToString();
      LOG_DBG3(pyopt) << " mat u_e^T=" << umatT.ToString();
      LOG_DBG3(pyopt) << " u_e*u_e^T=" << umat.ToString();

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

        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " bMat=" << bMat.ToString();
        // multiplications for derivative of compliance w. r. t. material tensor: B*u*u^T*B^T
        Matrix<double> uuTBT;
        uuTBT.Resize(u_e.GetSize(),rows);
        // u*u^T*B^T
        umat.Mult(bMatT,uuTBT);
        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " u*u^T*B^T=" << uuTBT.ToString();
        // B*u*u^T*B^T
        bMat *= uuTBT;

        res[e] -= bMat * weights[ip] * lp.jacDet;

  //      std::cout << "elem" <<  elem->elemNum << " ip=" << ip << " weight=" << weights[ip] << " jac det=" << lp.jacDet << std::endl;

        LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " ip=" << ip << " B*u*u^T*B^T=" << bMat.ToString();
        LOG_DBG3(pyopt) << " tmp=" << res[e].ToString();
      } // integration points
      LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " dfdH=" << res[e].ToString();
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
    LOG_DBG3(pyopt) << "GdfDH: elem=" << elem->elemNum << " dfdH=" << dfdH_e.ToString();
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

PyObject* PythonOptimizer::GetDims(PyObject* args)
{
  unsigned int dim = domain->GetDim();
  // assume rectilinear mesh with one region only
  const auto& regids = domain->GetOptimization()->GetDesign()->GetRegionIds();


  assert(regids.GetSize() == 1);
  StdVector<unsigned int> bounds = domain->GetGrid()->GetRegularDiscretization(regids.First());
  assert(bounds.GetSize() == 3);

  PyObject* res = PyTuple_New(4);
  PyTuple_SetItem(res, 0, PyLong_FromLong(dim));
  PyTuple_SetItem(res, 1, PyLong_FromLong(bounds[0]));
  PyTuple_SetItem(res, 2, PyLong_FromLong(bounds[1]));
  PyTuple_SetItem(res, 3, PyLong_FromLong(bounds[2]));

  return res;
}

/** number of design variables */
PyObject* PythonOptimizer::GetNumDesign(PyObject* args)
{
  return PyLong_FromLong(domain->GetOptimization()->GetDesign()->data.GetSize());
}


/** arg can be enum value or string */
Function::Access ParseAccess(PyObject* arg)
{
  Function::Access ret = Function::NO_ACCESS;
  if(PyLong_Check(arg))
    ret = (Function::Access) PyLong_AsLong(arg);
  else
  {
    const char* c_str = PyUnicode_AsUTF8(arg);
    ret = Function::access.Parse(c_str);
  }
  return ret;
}

double GetDesign(const BaseDesignElement* bde, Function::Access ac)
{
  switch(ac)
  {
  case Function::DEFAULT:
  case Function::NO_ACCESS:
  case Function::PLAIN:
    return bde->GetDesign(DesignElement::PLAIN);

  case Function::PHYSICAL:
  case Function::FILTERED:
  {
    const DesignElement* de = dynamic_cast<const DesignElement*>(bde);
    if(de == nullptr)
      return bde->GetDesign(); // slack can only be plain
    else
      return ac == Function::PHYSICAL ? de->GetPhysicalDesign() : de->GetDesign(DesignElement::SMART);
  }
  }
  assert(false);
  return -1;
}


/** return single plain design value */
PyObject* PythonOptimizer::GetDesignValue(PyObject* args)
{
  DesignSpace* space = domain->GetOptimization()->GetDesign();

  if(PyTuple_Size(args) < 1 || PyTuple_Size(args) > 2)
    throw("arguments for opt_get_design_value() are 0-based index and optionally access");

  unsigned int idx = PyLong_AsLong(PyTuple_GetItem(args,0));

  Function::Access ac = PyTuple_Size(args) == 2 ? ParseAccess(PyTuple_GetItem(args,1)) : Function::PLAIN;

  if(idx < 0 || idx >= space->GetNumberOfVariables())
    EXCEPTION("invalid design index " << idx);

  return PyFloat_FromDouble(GetDesign(space->GetDesignElement(idx), ac));
}


/** set design variables to provided numpy array of proper size (1. arg) with optional access (2.arg).
 * access can be int or string
 * @see BaseDesignElement::ValueSpecifier and BaseDesignElement::Access */
PyObject* PythonOptimizer::GetDesignValues(PyObject* args)
{
  DesignSpace* space = domain->GetOptimization()->GetDesign();

  if(PyTuple_Size(args) < 1 || PyTuple_Size(args) > 2)
    throw("arguments for opt_get_design_values() are numpy-return, access (optional)");

  Function::Access ac = PyTuple_Size(args) == 2 ? ParseAccess(PyTuple_GetItem(args,1)) : Function::PLAIN;

  LOG_DBG(pyopt) << "GDVs access=" << ac << "=" << Function::access.ToString(ac) << " n=" << space->GetNumberOfVariables();

  Vector<double> py(space->GetNumberOfVariables());

  for(unsigned int i = 0, n = space->GetNumberOfVariables(); i < n; i++)
    py[i] = GetDesign(space->GetDesignElement(i), ac);

  py.Export(PyTuple_GetItem(args,0)); // PyTuple_GetItem = Borrowed reference

  Py_RETURN_NONE;
}

PyObject* PythonOptimizer::Transfer(PyObject* args, bool derivative)
{
  DesignSpace* space = domain->GetOptimization()->GetDesign();
  assert(space);

  if(PyTuple_Size(args) < 1 || PyTuple_Size(args) > 2)
    throw "arguments for opt_(d_)_transfer() are value, index (optional)";

  int idx =  PyTuple_Size(args) == 2 ? PyLong_AsLong(PyTuple_GetItem(args,1)) : 0;

  if(idx < 0 || idx > (int) space->transfer.GetSize()-1)
    throw "invalid index for opt_(d_)_transfer() " + to_string(idx) + " range is 0 ... " + to_string(space->transfer.GetSize()-1);

  double val = PyFloat_AsDouble(PyTuple_GetItem(args,0));

  if(derivative)
    return PyFloat_FromDouble(space->transfer[idx].Derivative(val));
  else
    return PyFloat_FromDouble(space->transfer[idx].Transform(val));
}

} // namespace
