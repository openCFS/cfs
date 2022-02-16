#include <fstream>
#include <forward_list>
#include <boost/tokenizer.hpp>
#include <iostream>

#include "Domain/CoefFunction/CoefFunctionPython.hh"
#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::string;

namespace CoupledField
{

CoefFunctionPython::CoefFunctionPython(PtrParamNode pn, unsigned int dim) : CoefFunction()
{
  isAnalytic_ = false;
  isComplex_ = false;
  supportDerivative_ = false;
  dimType_ = dim == 1 ? SCALAR : VECTOR;
  dependType_ = SPACE;


  if(python->GetKernel() == NULL)
    throw Exception("To use python coef function, define a python kernel module in document root python element.");

  dim_      = dim;
  init_     = pn->Get("init")->As<string>();
  function_ = pn->Get("function")->As<string>();

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
  Vector<double> pointCoord;
  lpm.GetGlobal(pointCoord); // either from shapeMap or from lpm.lp.number

  PyObject* arg = PyTuple_New(1);
  PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonList(pointCoord)); // PyTuple_SetItem() steals the reference

  assert(eval_ != NULL);
  PyObject* ret = PyObject_CallObject(eval_, arg);
  PythonKernel::CheckPythonReturn(ret);

  Py_XDECREF(arg);

  return ret; // dont't forget to decref ret
}

void CoefFunctionPython::GetVector(Vector<double>& vec, const LocPointMapped& lpm)
{
  PyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  PythonKernel::ConvertPythonList<double>(vec, ret);
  Py_XDECREF(ret);
}


void CoefFunctionPython::GetScalar(double& val, const LocPointMapped& lpm)
{
  PyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  val = PyFloat_AsDouble(ret);
  Py_XDECREF(ret);
}


} // end of namespace
