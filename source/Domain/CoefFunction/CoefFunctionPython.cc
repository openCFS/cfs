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
    pyObject* func = pyObject_GetAttrString(python->GetKernel(), init_.c_str());
    PythonKernel::CheckPythonFunction(func, init_.c_str());

    pyObject* arg = PyTuple_New(1);
    PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonDict(opt)); // PyTuple_SetItem() steals the reference

    pyObject* call = pyObject_CallObject(func, arg);
    PythonKernel::CheckPythonReturn(call);

    Py_XDECREF(call);
    Py_XDECREF(arg);
    Py_XDECREF(func);
  }

  eval_ = pyObject_GetAttrString(python->GetKernel(), function_.c_str());
  PythonKernel::CheckPythonFunction(eval_, function_.c_str());

}

CoefFunctionPython::~CoefFunctionPython()
{
  Py_XDECREF(eval_);
}


pyObject* CoefFunctionPython::CallFunction(const LocPointMapped& lpm)
{
  pyObject* arg = PyTuple_New(1);

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
  pyObject* ret = pyObject_CallObject(eval_, arg);
  PythonKernel::CheckPythonReturn(ret);

  Py_XDECREF(arg);

  return ret; // dont't forget to decref ret
}

void CoefFunctionPython::GetVector(Vector<double>& vec, const LocPointMapped& lpm)
{
  pyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  PythonKernel::ConvertPythonList<double>(vec, ret);
  Py_XDECREF(ret);
}


void CoefFunctionPython::GetScalar(double& val, const LocPointMapped& lpm)
{
  pyObject* ret = CallFunction(lpm);
  assert(ret != NULL);

  val = PyFloat_AsDouble(ret);
  Py_XDECREF(ret);
}


} // end of namespace
