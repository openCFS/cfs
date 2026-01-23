#include <fstream>
#include <forward_list>
#include <mutex>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "Domain/CoefFunction/CoefFunctionPython.hh"
#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::string;

namespace CoupledField
{

// Mutex to serialize Python API calls from multiple threads
// Using mutex instead of GIL to avoid deadlock when main thread holds GIL
// and worker threads try to acquire it during OpenMP parallel regions
static std::mutex python_call_mutex;

CoefFunctionPython::CoefFunctionPython(PtrParamNode pn, unsigned int dim) : CoefFunction()
{
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = dim == 1 ? SCALAR : VECTOR;
  dependType_ = SPACE;


  if(python->GetKernel() == NULL)
    throw Exception("To use python coef function, define a python kernel module in document root python element.");

  dim_       = dim;
  init_      = pn->Get("init")->As<string>();
  function_  = pn->Get("function")->As<string>();
  normalize_ = pn->Get("normalize")->As<bool>();
  by_coord_  = pn->Get("by_coord")->As<bool>();

  auto opt = PythonKernel::ParseOptions(pn->GetList("option"));
  if(opt.GetSize() > 0 && init_ == "")
    throw Exception("Options cannot be used as init function is not given");

  if(init_ != "")
  {
    PyObject* func = PyObject_GetAttrString(python->GetKernel(), init_.c_str());
    PythonKernel::CheckPythonFunction(func, init_.c_str());

    PyObject* arg = PyTuple_New(1);
    PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonDict(opt)); // PyTuple_SetItem() steals the reference

    PyObject* call = PyObject_CallObject(func, arg);
    PythonKernel::CheckPythonReturn(call);

    Py_XDECREF(call);
    Py_XDECREF(arg);
    Py_XDECREF(func);
  }

  eval_ = PyObject_GetAttrString(python->GetKernel(), function_.c_str());
  PythonKernel::CheckPythonFunction(eval_, function_.c_str());

}

CoefFunctionPython::~CoefFunctionPython()
{
  Py_XDECREF(eval_);
}


PyObject* CoefFunctionPython::CallFunction(const LocPointMapped& lpm)
{
  PyObject* arg = PyTuple_New(1);

  if(by_coord_) {
    Vector<double> pointCoord;
    lpm.GetGlobal(pointCoord); // either from shapeMap or from lpm.lp.number
    PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonList(pointCoord)); // PyTuple_SetItem() steals the reference
  }
  else {
    assert(lpm.lp.number != LocPoint::NOT_SET);
    unsigned long number = lpm.lp.number;
    PyTuple_SetItem(arg, 0, PyLong_FromUnsignedLong(number));
    if(number == 1)
    {
      //std::cout << "CoefFunctionPython::CallFunction(" << number << ")\n";
      assert(arg != nullptr); // to set breakpoint
    }
  }

  assert(eval_ != NULL);
  PyObject* ret = PyObject_CallObject(eval_, arg);
  PythonKernel::CheckPythonReturn(ret);

  Py_XDECREF(arg);

  return ret; // dont't forget to decref ret
}

void CoefFunctionPython::GetVector(Vector<double>& vec, const LocPointMapped& lpm)
{
  // Serialize Python API calls using mutex instead of GIL
  // Using mutex avoids deadlock when main thread holds GIL and enters OpenMP region
  std::lock_guard<std::mutex> lock(python_call_mutex);

  PyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  PythonKernel::ConvertPythonList<double>(vec, ret);
  Py_XDECREF(ret);
}


void CoefFunctionPython::GetScalar(double& val, const LocPointMapped& lpm)
{
  // Serialize Python API calls using mutex instead of GIL
  // Using mutex avoids deadlock when main thread holds GIL and enters OpenMP region
  std::lock_guard<std::mutex> lock(python_call_mutex);

  PyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  val = PyFloat_AsDouble(ret);
  Py_XDECREF(ret);
}


} // end of namespace
