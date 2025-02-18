/** this file contains the embedded python specific implementations for Function, ErsatzMaterial, ... */

#include "General/Exception.hh"
#include "Optimization/Function.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Optimizer/OptimalityCondition.hh"
#include "Optimization/Optimizer/MMA.hh"
#include "Optimization/Optimizer/DumasMMA.hh"
#include <MMASolver.h>
#include <GCMMASolver.h>
#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

EXTERN_LOG(func) // opt_func
EXTERN_LOG(em)
EXTERN_LOG(designSpace)

namespace CoupledField
{

void Optimization::PythonStopOptimization(PyObject* args)
{
  int ok = -1;
  char* msg = nullptr;

  PyArg_ParseTuple(args, "ps", &ok, &msg);

  this->user_break_converged = ok;
  this->user_break_reason = std::string(msg);
}

PyObject* Optimization::PythonFunctionValues() const
{
  // from Optimization::LogFileLine()
  StdVector<std::pair<std::string, std::string> > objs;
  for(Function* f : objectives.data)
    objs.Push_back(std::make_pair(f->ToString(), std::to_string(f->GetValue())));

  StdVector<std::pair<std::string, std::string> > cnstrs;
  for(Function* g : constraints.all)
    cnstrs.Push_back(std::make_pair(g->ToString(), std::to_string(g->GetValue())));

  PyObject* ret = PyTuple_New(2);
  PyTuple_SetItem(ret, 0, PythonKernel::CreatePythonDict(objs));
  PyTuple_SetItem(ret, 1, PythonKernel::CreatePythonDict(cnstrs));
  return ret;
}

PyObject* Optimization::PythonFunctionProperties(PyObject* args)
{
  char* ns = nullptr;
  PyArg_ParseTuple(args, "s", &ns);
  Function* f = GetFunction(string(ns));
  StdVector<std::pair<string, string> > map;
  f->DescribeProperties(map);
  return PythonKernel::CreatePythonDict(map);
}


PyObject* Optimization::PythonGetOptimizerProperties() const
{
  StdVector<std::pair<std::string, std::string> > map;

  map.Push_back(std::make_pair("optimizer", optimizer.ToString(optimizer_)));
  map.Push_back(std::make_pair("max_iter", std::to_string(maxIterations)));
  baseOptimizer_->DescribeProperties(map);
  return PythonKernel::CreatePythonDict(map);
}

std::pair<string, string> ParseStringString(PyObject* args)
{
  char* fc = nullptr;
  char* sc = nullptr;

  PyArg_ParseTuple(args, "ss", &fc, &sc);

  return std::make_pair(string(fc),string(sc));
}


void OptimalityCondition::PythonSetProperty(PyObject* args)
{
  auto ss = ParseStringString(args);

  if(ss.first == "damping")
    oc_damping_ = std::stod(ss.second);
  else if(ss.first == "move_limit")
    move_limit_ = std::stod(ss.second);
  else
    throw "Unknown property " + ss.first + " for 'ocm'";
}


void MMA::PythonSetProperty(PyObject* args)
{
  auto ss = ParseStringString(args);

  if(ss.first == "move_limit")
    move_limit = std::stod(ss.second);
  else
    throw "Unknown property " + ss.first + " for 'MMA'";
}

void DumasMMA::PythonSetProperty(PyObject* args)
{
  auto ss = ParseStringString(args);

  if(ss.first == "move_limit")
    move_limit = std::stod(ss.second);
  else if(ss.first == "asymdec")
    asymdec = std::stod(ss.second);
  else if(ss.first == "asyminc")
    asyminc = std::stod(ss.second);
  else
    throw "Unknown property " + ss.first + " for 'DumasMMA'";

  // no harm, when we do this also for move_limit, the functions just sets the internal variables without further action
  if(mma)
    mma->SetAsymptotes(asyminit, asymdec, asyminc);
  else
    gcmma->SetAsymptotes(asyminit, asymdec, asyminc);
}

void Function::InitPythonFunction(PtrParamNode pn, DesignSpace* design)
{
  assert(type_ == PYTHON_FUNCTION || type_ == LOCAL_PYTHON_FUNCTION);

  if(!pn->Has("python")) // this is the sub-element
    throw Exception("python function requires 'python' element");

  PtrParamNode ppn = pn->Get("python");

  // such that a later Function::ToInfo() call has the properties
  preInfo_->Get("python")->SetValue(ppn);

  py_name_= ppn->Get("name")->As<string>();

  // we have as kernel either kernel or design or it is possibly not set
  PyObject* kernel = ppn->Get("script")->As<string>() == "kernel" ? python->GetKernel() : design->GetPythonModule();
  if(kernel == nullptr)
    throw Exception("No python file set. Possibly set it via the global python element");

  // the options are flattened and contain the function and the python element. Going deeper is ugly
  auto list = pn->ToStringList(2);
  // possibly there are options
  auto optlist = PythonKernel::ParseOptions(pn->Get("python")->GetList("option"));
  list.Append(optlist);

  PyObject* opt = PythonKernel::CreatePythonDict(list);
  assert(opt != NULL);
  py_arg_ = PyTuple_New(1); // keep this for all python function calls
  Py_INCREF(opt); // to compensate SetItem()
  PyTuple_SetItem(py_arg_, 0, opt); // PyTuple_SetItem() steals the reference


  // init is optional
  if(ppn->Has("init"))
  {
    PyObject* func = PyObject_GetAttrString(kernel, ppn->Get("init")->As<string>().c_str());
    if(func)
    {
      PyObject* call = PyObject_CallObject(func, py_arg_);
      PythonKernel::CheckPythonReturn(call);

      Py_XDECREF(call);
      Py_XDECREF(func);
    }
  }

  // create eval and grad pointers
  py_eval_ = PyObject_GetAttrString(kernel, ppn->Get("eval")->As<string>().c_str());
  PythonKernel::CheckPythonFunction(py_eval_, ppn->Get("eval")->As<string>().c_str());

  py_grad_ = PyObject_GetAttrString(kernel, ppn->Get("grad")->As<string>().c_str());
  PythonKernel::CheckPythonFunction(py_grad_, ppn->Get("grad")->As<string>().c_str());

  // we check this in SetLocalPythonVirtualElementMap()
  if(ppn->Has("sparsity"))
    py_sparsity_ = PyObject_GetAttrString(kernel, ppn->Get("sparsity")->As<string>().c_str());
}


void Function::DeletePythonFunction()
{
  Py_XDECREF(py_sparsity_);
  Py_XDECREF(py_grad_);
  Py_XDECREF(py_eval_);
  Py_XDECREF(py_arg_);

}

void Function::SetLocalPythonVirtualElementMap(StdVector<Function::Local::Identifier>& virtual_elem_map, DesignSpace* space)
{
  if(!py_sparsity_)
    throw Exception("localPython requires python/sparsity to be given");

  PythonKernel::CheckPythonFunction(py_sparsity_);
  PyObject* call = PyObject_CallObject(py_sparsity_, py_arg_);
  PythonKernel::CheckPythonReturn(call);

  StdVector<Vector<int> > vec;
  PythonKernel::ConvertPythonList(vec, call); // does not decref the numpy arrays

  virtual_elem_map.Reserve(vec.GetSize());
  StdVector<BaseDesignElement*> buddies;

  for(const Vector<int>& item : vec)
  {
    LOG_DBG2(func) << "SLPVEM: item=" << item.ToString();
    assert(item.GetSize() > 1);
    BaseDesignElement* de = space->GetDesignElement(item[0]);

    buddies.Resize(0); // keeps capacity
    for(unsigned int i = 1; i < item.GetSize(); i++)
      buddies.Push_back(space->GetDesignElement(item[i]));

    virtual_elem_map.Push_back(Local::Identifier(de, buddies));
  }


  Py_XDECREF(call);
}


PyObject* Function::CallPythonFunction(bool eval)
{
  PyObject* func = eval ? py_eval_ : py_grad_;
  PyObject* call = PyObject_CallObject(func, py_arg_);
  PythonKernel::CheckPythonReturn(call);
  return call;
}


double ErsatzMaterial::CalcPython(Excitation& excite, Function* f, bool derivative)
{
  PyObject* pyret = f->CallPythonFunction(!derivative);
  double ret = -1;
  if(!derivative)
  {
    ret = PyFloat_AsDouble(pyret);
  }
  else
  {
    Vector<double> grad(pyret, false);
    if(grad.GetSize() != f->elements.GetSize())
      EXCEPTION("got gradient of size " << grad.GetSize() << " excpected " << f->elements.GetSize());

    for(unsigned int i = 0; i < grad.GetSize(); i++)
    {
      DesignElement* de = f->elements[i];
      de->AddGradient(f, grad[i]);
      LOG_DBG2(em) << "CP: i=" << i << " pygrad=" << grad[i];
    }
  }

  Py_XDECREF(pyret);
  return ret;
}

PyObject* CallLocalPythonFunction(const Function* func, PyObject* py_func)
{
  int idx = dynamic_cast<const LocalCondition*>(func)->GetCurrentRelativePosition();

  PyObject* arg = PyTuple_New(1);
  PyTuple_SetItem(arg, 0, PyLong_FromLong(idx)); // PyTuple_SetItem() steals the reference
  PyObject* call = PyObject_CallObject(py_func, arg);
  PythonKernel::CheckPythonReturn(call);
  Py_XDECREF(arg);
  return call;
}

/** implemented in PythonFunction.cc */
double Function::Local::Identifier::CalcLocalPythonFunc(const Function::Local* local) const
{
  PyObject* call = CallLocalPythonFunction(local->func_, local->func_->py_eval_);

  double ret = PyFloat_AsDouble(call);
  Py_XDECREF(call);
  return ret;
}

/** we return the "full" local python function gradient at once. Implemented in PythonFunction.cc */
void Function::Local::Identifier::CalcLocalPythonGrad(Vector<double>& grad, const Function::Local* local) const
{
  PyObject* call = CallLocalPythonFunction(local->func_, local->func_->py_grad_);
  grad.Fill(call, true); // decref
}

PyObject* DesignSpace::PythonGetFilterProperties(PyObject* args) const
{
  if(filter.GetSize() == 0)
    throw("it seams no filter is set for optimization");

  if(PyTuple_Size(args) > 2)
    throw("max one optional argument for get_opt_filter_values() / DesignSpace::PythonGetFilterProperties()");

  unsigned int idx = PyTuple_Size(args) == 1 ? PyLong_AsLong(PyTuple_GetItem(args,0)) : 0;

  if(idx >= filter.GetSize())
    throw("largest argument for get_opt_filter_values(idx) is " + std::to_string(filter.GetSize()-1) + " given is " + std::to_string(idx));

  const GlobalFilter& gf = filter[idx];

  StdVector<std::pair<std::string, std::string> > map;
  // values which can be set
  map.Push_back(std::make_pair("value", std::to_string(gf.value)));
  map.Push_back(std::make_pair("beta", std::to_string(gf.beta)));
  map.Push_back(std::make_pair("eta", std::to_string(gf.eta)));
  map.Push_back(std::make_pair("non_lin_scale", std::to_string(gf.non_lin_scale)));
  map.Push_back(std::make_pair("non_lin_offset", std::to_string(gf.non_lin_offset)));
  // types info
  map.Push_back(std::make_pair("type", Filter::type.ToString(gf.type)));
  map.Push_back(std::make_pair("filterSpace", Filter::filterSpace.ToString(gf.filterspace)));
  map.Push_back(std::make_pair("sensitivity", Filter::sensitivity.ToString(gf.sensitivity)));
  map.Push_back(std::make_pair("density", Filter::density.ToString(gf.density)));
  // design and statistics
  map.Push_back(std::make_pair("design", BaseDesignElement::type.ToString(gf.design)));
  map.Push_back(std::make_pair("region", domain->GetGrid()->regionData[gf.region].name));
  map.Push_back(std::make_pair("elements", std::to_string(gf.elements)));
  map.Push_back(std::make_pair("avg_radius", std::to_string(gf.avg_radius)));
  map.Push_back(std::make_pair("avg_neigbor", std::to_string(gf.avg_neigbor)));
  // little dirty helper :)
  map.Push_back(std::make_pair("total_filters", std::to_string(filter.GetSize())));

  return PythonKernel::CreatePythonDict(map);
}

void DesignSpace::PythonSetFilterProperties(PyObject* args)
{
  int idx = -1;
  double beta = -1.0;
  double eta = -1.0;
  double scale = -1.0;
  double offset = -1e42;

  LOG_DBG(designSpace) << "PSFP args=" << PyTuple_Size(args);
  if(PyTuple_Size(args) < 2 || PyTuple_Size(args) > 5)
      throw("set_opt_filter_values(idx, beta, [eta, [scale, [offset]]]) has 2 ... 5 arguments, keywords are not allowed");
  int ok = PyArg_ParseTuple(args,"id|ddd",&idx, &beta, &eta, &scale, &offset);
  LOG_DBG(designSpace) << "PSFP args=" << PyTuple_Size(args) << " ok=" << ok << " beta=" << beta << " eta=" << eta << " scale=" << scale << " offset=" << offset;
  PythonKernel::CheckPythonReturn(ok, "set_opt_filter_values");

  if(idx >= (int) filter.GetSize())
    throw("largest idx for set_opt_filter_values() is " + std::to_string(filter.GetSize()-1) + " given is " + std::to_string(idx));

  GlobalFilter& gf = filter[idx];

  gf.beta = beta;
  if(eta != -1.0)
    gf.eta = eta;
  if(scale != -1.0 && offset != -1e42)
  {
    gf.non_lin_scale = scale;
    gf.non_lin_offset = offset;
  }
  else
  {
    // from DesignStructure::SetFilter()
    DesignRegion* dr = GetRegion(gf.region, gf.design);
    gf.SetNonLinCorrection(&data[dr->base]);
    LOG_DBG(designSpace) << "PSFP ref=" << dr->base << " SetNonLinCorrection -> scale=" << gf.non_lin_scale << " offset=" << gf.non_lin_offset;
  }
}



} // end of namespace
