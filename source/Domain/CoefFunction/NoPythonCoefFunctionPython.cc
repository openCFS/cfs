// this implements dummpy functions for when not compiled with USE_EMBEDDED_PYTHON

#include "Domain/CoefFunction/CoefFunctionPython.hh"
using std::string;

namespace CoupledField
{

CoefFunctionPython::CoefFunctionPython(PtrParamNode pn, unsigned int dim) : CoefFunction()
{
  throw Exception("CoefFunctionPython requires USE_EMBEDDED_PYTHON");
}

CoefFunctionPython::~CoefFunctionPython()
{
}


pyObject* CoefFunctionPython::CallFunction(const LocPointMapped& lpm)
{
  return NULL;
}

void CoefFunctionPython::GetVector(Vector<double>& vec, const LocPointMapped& lpm)
{
}


void CoefFunctionPython::GetScalar(double& val, const LocPointMapped& lpm)
{
}


} // end of namespace
