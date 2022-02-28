#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/SimInOut/python/SimInputPython.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"

/** this file contains methods implementation for the python api provided to the python
 * instance controlled by PythonKernel */

namespace CoupledField
{

PyObject* PythonKernel::mesher_set_nodes(PyObject *self, PyObject *args)
{
  if(!CheckMesher())
    return NULL;

  mesher_->SetNodes(args);
  Py_RETURN_NONE;
}

PyObject* PythonKernel::mesher_set_regions(PyObject *self, PyObject *args)
{
  if(!CheckMesher())
    return NULL;

  mesher_->SetRegions(args);
  Py_RETURN_NONE;
}

PyObject* PythonKernel::mesher_add_elements(PyObject *self, PyObject *args)
{
  if(!CheckMesher())
    return NULL;

  mesher_->AddElements(args);
  Py_RETURN_NONE;
}

PyObject* PythonKernel::mesher_add_named_nodes(PyObject *self, PyObject *args)
{
  if(!CheckMesher())
    return NULL;

  mesher_->AddNamedNodes(args);
  Py_RETURN_NONE;
}

PyObject* PythonKernel::mesher_add_named_elements(PyObject *self, PyObject *args)
{
  if(!CheckMesher())
    return NULL;

  mesher_->AddNamedElements(args);
  Py_RETURN_NONE;
}

/** cfs.bound(xl,xu,gl,gu) sets bounds for design and constraints in properly sized 1d numpy arrays */
PyObject* opt_getDims(PyObject *self, PyObject *args)
{
  if(!PythonKernel::CheckOpt())
    return NULL;

  return PythonOptimizer::GetDims(args);
}

PyObject* get_opt_design_size(PyObject *self, PyObject *args)
{
  if(!PythonKernel::CheckOpt())
    return NULL;

  return PythonOptimizer::GetNumDesign(args);
}

PyObject* get_opt_design_value(PyObject *self, PyObject *args)
{
  if(!PythonKernel::CheckOpt())
    return NULL;

  return PythonOptimizer::GetDesignValue(args);
}

PyObject* get_opt_design_values(PyObject *self, PyObject *args)
{
  if(!PythonKernel::CheckOpt())
    return NULL;

  return PythonOptimizer::GetDesignValues(args);
}


/** cfs.bound(xl,xu,gl,gu) sets bounds for design and constraints in properly sized 1d numpy arrays */
PyObject* opt_bounds(PyObject *self, PyObject *args)

{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->GetBounds(args);
  Py_RETURN_NONE;
}

/** cfs.initialdesign(x) fills 1d numpy array */
PyObject* opt_initialdesign(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->GetInitialDesign(args);
  Py_RETURN_NONE;
}

/** cfs.evalobj(x) returns a float) */
PyObject* opt_evalobj(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  return PyFloat_FromDouble(python->GetPyOpt()->EvalObjective(args));
}

/** cfs.cfs_commitIteration() commits iteration to cfs */
PyObject* opt_commitIteration(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->CommitIteration(); // return paramnode ignored;

  Py_RETURN_NONE;
}

/** cfs.evalgradobj(x,grad) fills numpy 1d arrays of size n */
PyObject* opt_evalgradobj(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->EvalGradObjective(args);
  Py_RETURN_NONE;
}

/** cfs.cfs_evalconstrs(x,g) fills numpy 1d arrays of size n and m */
PyObject* opt_evalconstrs(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->EvalConstraints(args);
  Py_RETURN_NONE;
}

/** cfs.cfs_evalgradconstrs(x,grad) fills numpy 1d arrays of size n and m*n */
PyObject* opt_evalgradconstrs(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->EvalGradConstraints(args);
  Py_RETURN_NONE;
}

PyObject* opt_getSimpExponent(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  return PyFloat_FromDouble(python->GetPyOpt()->GetSimpExponent());
}

/** returns derivative of compliance w.r.t stiffness tensor entries of original (core) material */
PyObject* opt_get_dfdH(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->Get_dfdH(args);
  Py_RETURN_NONE;
}

/** cfs.cfs_getOrgStiffness(stiffness) returns stiffness tensor of original (core) material */
PyObject* opt_getOrgStiffness(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  python->GetPyOpt()->GetCoreStiffness(args);
  Py_RETURN_NONE;
}

/** return true if cfs's stopping criteria is met, including finding the file HALTOPT */
PyObject* opt_dostop(PyObject *self, PyObject *args)
{
  if(!python->CheckPyOpt())
    return NULL;

  bool stop = domain->GetOptimization()->DoStopOptimization();

  if(stop)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

// pde stuff

/** return list with single pde names */
PyObject* single_pde_names(PyObject *self, PyObject *args)
{
  assert(domain != NULL);

  const StdVector<SinglePDE*>& pdes = domain->GetSinglePDEs();
  StdVector<string> vec(pdes.GetSize());
  for(unsigned int i = 0; i < pdes.GetSize(); i++)
    vec[i] = pdes[i]->GetName();

  return PythonKernel::CreatePythonList(vec);
}

PyMethodDef PythonKernel::cfs_methods[] =
{
  /* python mesh reader */
  {"set_nodes", PythonKernel::mesher_set_nodes, METH_VARARGS, "Set node coordinates via numpy array with 2 or 3 columns"},
  {"set_regions", PythonKernel::mesher_set_regions, METH_VARARGS, "set region names by list of strings"},
  {"add_elements", PythonKernel::mesher_add_elements, METH_VARARGS, "parameters total_number, fe_type, int numpy array of 1-based node ids"},
  {"add_named_nodes", PythonKernel::mesher_add_named_nodes, METH_VARARGS, "string and int numpy array of 1-based node ids"},
  {"add_naned_elements", PythonKernel::mesher_add_named_elements, METH_VARARGS, "string and int numpy array of 1-based element ids"},

  /* general optimization */
  {"getDims", opt_getDims, METH_VARARGS, "Returns info on optimization design domain dimensions: dim, nx, ny, nz."},
  {"get_opt_design_size", get_opt_design_size, METH_VARARGS, "Returns number of design variables"},
  {"get_opt_design_value", get_opt_design_value, METH_VARARGS, "Give single DesignSpace::GetDesignValue() for 0-based index with optional access (default plain)"},
  {"get_opt_design_values", get_opt_design_values, METH_VARARGS, "Call DesignElement::GetValues() with attributes numpy array and optional access"},

  /** python optimizer */
  {"bounds", opt_bounds, METH_VARARGS, "Give design and constraints bounds. Expects 1D arrays for xl, xu, gl, gu"},
  {"initialdesign", opt_initialdesign, METH_VARARGS, "Sets the initial design to the given 1D array"},
  {"evalobj", opt_evalobj, METH_VARARGS, "Evaluate objective. Expects 1D design array"},
  {"evalconstrs", opt_evalconstrs, METH_VARARGS, "Evaluate constraints. Expects 1D design array and 1D value array or size m"},
  {"evalgradobj", opt_evalgradobj, METH_VARARGS, "Evaluate objective gradient w. r. t design variables. Expects 1D design array and dense gradient array "},
  {"evalgradconstrs", opt_evalgradconstrs, METH_VARARGS, "Evaluate constrains gradients. Expects 1D design array of size n and dense 1D value array of size m*n"},
  {"getSimpExponent", opt_getSimpExponent, METH_VARARGS, "Returns power-law exponent."},
  {"get_dfdH", opt_get_dfdH, METH_VARARGS, "Returns derivative of objective w.r.t. material tensor. Expects 3D array of size nelems x 3 x 3"},
  {"getOrgStiffness", opt_getOrgStiffness, METH_VARARGS, "Returns stiffness tensor of original (core) material."},
  {"commitIteration", opt_commitIteration, METH_VARARGS, "Commit current iteration to cfs."},
  {"dostop", opt_dostop, METH_VARARGS, "Uses stopping criteria from cfs (including HALTOPT file) and returns True or False"},

  /** pde stuff */
  {"single_pde_names", single_pde_names, METH_VARARGS, "return a list with the names of the single pdes"},

  {NULL, NULL, 0, NULL}
};

PyModuleDef PythonKernel::cfs_modules = {
    PyModuleDef_HEAD_INIT, "cfs", NULL, -1, cfs_methods, NULL, NULL, NULL, NULL
};

PyObject* PythonKernel::SetModulFunctions(void)
{
  // https://stackoverflow.com/questions/37943699/crash-when-calling-pyarg-parsetuple-on-a-numpy-array
  //import_array();

  return PyModule_Create(&PythonKernel::cfs_modules);
}


} // end of namespace
