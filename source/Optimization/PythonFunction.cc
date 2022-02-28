/** this file contains the embedded python specific implementations for Function and ErsatzMaterial */

#include "General/Exception.hh"
#include "Optimization/Function.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

EXTERN_LOG(func)
EXTERN_LOG(em)

namespace CoupledField
{

void Function::InitPythonFunction(PtrParamNode pn, DesignSpace* design)
{
  assert(type_ == PYTHON_FUNCTION || type_ == LOCAL_PYTHON_FUNCTION);

  if(!pn->Has("python")) // this is the sub-element
    throw Exception("python function requires 'python' element");

  PtrParamNode ppn = pn->Get("python");

  // such that a later Function::ToInfo() call has the properties
  preInfo_->Get("python")->SetValue(ppn);

  py_name_= ppn->Get("name")->As<string>();



  // we have as kernel either kernel or design
  PyObject* kernel = ppn->Get("script")->As<string>() == "kernel" ? python->GetKernel() : design->GetPythonModule();

  // the options are flattened and contain everything
  StdVector<std::pair<std::string, std::string> > list;
  pn->ToStringList(list);

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



} // end of namespace
