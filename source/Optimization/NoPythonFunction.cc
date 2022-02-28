/** this file implements dummy replacements for PythonFunction.cc if compiled w/o USE_PYTHON_EXCEPTION */

#include "Optimization/Function.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Design/DesignSpace.hh"

namespace CoupledField
{

void Function::InitPythonFunction(PtrParamNode pn, DesignSpace* design)
{
  assert(false);
}


void Function::DeletePythonFunction()
{
  assert(false);
}

void Function::SetLocalPythonVirtualElementMap(StdVector<Function::Local::Identifier>& virtual_elem_map, DesignSpace* space)
{
  assert(false);
}


PyObject* Function::CallPythonFunction(bool eval)
{
  assert(false);
  return NULL;
}

double ErsatzMaterial::CalcPython(Excitation& excite, Function* f, bool derivative)
{
  assert(false);
  return -1;
}


/** implemented in PythonFunction.cc */
double Function::Local::Identifier::CalcLocalPythonFunc(const Function::Local* local) const
{
  assert(false);
  return -1;
}
/** we return the "full" local python function gradient at once. Implemented in PythonFunction.cc */
void Function::Local::Identifier::CalcLocalPythonGrad(Vector<double>& grad, const Function::Local* local) const
{
  assert(false);
}


} // end of namespace
