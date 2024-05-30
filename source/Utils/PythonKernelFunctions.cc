#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include "Utils/PythonKernel.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/SimInOut/python/SimInputPython.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"
#include "Optimization/Design/FeaturedDesign.hh"
#include "Driver/Assemble.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "OLAS/algsys/AlgebraicSys.hh"

DEFINE_LOG(pkf, "PythonKernelFunctions")

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
  return PythonKernel::CheckOpt() ? PythonOptimizer::GetDims(args) : NULL;
}

PyObject* get_opt_design_size(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? PythonOptimizer::GetNumDesign(args) : NULL;
}

PyObject* get_opt_design_value(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? PythonOptimizer::GetDesignValue(args) : NULL;
}

PyObject* get_opt_design_values(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? PythonOptimizer::GetDesignValues(args) : NULL;
}

PyObject* get_opt_iteraton(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? PyLong_FromLong(domain->GetOptimization()->GetCurrentIteration()) : NULL;
}

PyObject* get_opt_stopping_rules(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  return PythonKernel::CreatePythonDict(domain->GetOptimization()->GetStoppingRules());
}

PyObject* get_opt_function_values(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? domain->GetOptimization()->PythonFunctionValues() : NULL;
}

PyObject* get_opt_function_properties(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? domain->GetOptimization()->PythonFunctionProperties(args) : NULL;
}

PyObject* get_opt_filter_values(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? domain->GetOptimization()->GetDesign()->PythonGetFilterProperties(args) : NULL;
}

PyObject* set_opt_filter_values(PyObject *self, PyObject *args)
{
  if(!PythonKernel::CheckOpt())
    return NULL;

  domain->GetOptimization()->GetDesign()->PythonSetFilterProperties(args);
  Py_RETURN_NONE;
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
  return python->CheckPyOpt() ? PyFloat_FromDouble(python->GetPyOpt()->EvalObjective(args)) : NULL;
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
  return python->CheckPyOpt() ? PyFloat_FromDouble(python->GetPyOpt()->GetSimpExponent()) : NULL;
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

/** make optimization stop in the next iteration.
 * param args tuple of bool and string with info if converged and the message */
PyObject* opt_stop(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  domain->GetOptimization()->PythonStopOptimization(args);
  Py_RETURN_NONE;
}

/** return true if cfs's stopping criteria is met, including finding the file HALTOPT */
PyObject* opt_dostop(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  bool stop = domain->GetOptimization()->DoStopOptimization();

  if(stop)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

PyObject* opt_register_log_property(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  char* ns = nullptr;
  char* vs = nullptr;

  PyArg_ParseTuple(args, "ss", &ns, &vs);
  domain->GetOptimization()->RegisterAuxLogValue(std::string(ns), std::string(vs));
  Py_RETURN_NONE;
}

PyObject* opt_set_log_property(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  char* ns = nullptr;
  char* vs = nullptr;

  PyArg_ParseTuple(args, "ss", &ns, &vs);
  domain->GetOptimization()->SetAuxLogValue(std::string(ns), std::string(vs));
  Py_RETURN_NONE;
}

PyObject* optimizer_get_properties(PyObject *self, PyObject *args)
{
  return PythonKernel::CheckOpt() ? domain->GetOptimization()->PythonGetOptimizerProperties() : NULL;
}

PyObject* optimizer_set_property(PyObject *self, PyObject *args)
{
  if(!python->CheckOpt())
    return NULL;

  domain->GetOptimization()->GetOptimizerInstance()->PythonSetProperty(args);
  Py_RETURN_NONE;
}

/** feature mapping stuff */

/** the total number of feature mapping parameters - not only their optimization variables */
PyObject* feature_mapping_num_parameters(PyObject *self, PyObject *args)
{
  return PyLong_FromLong(PythonKernel::GetFeaturedDesign()->GetNumberOfFeatureMappingVariables());
}

/** fills the given array with feature mapping parameters. Including fixed ones and no aux or design */
PyObject* feature_mapping_get_parameters(PyObject *self, PyObject *args)
{
  FeaturedDesign* fd = PythonKernel::GetFeaturedDesign(); // throws

  PyArrayObject* pyarray = nullptr;

  /** temporary data - not too big as feature mapping */
  Vector<double> data(fd->GetNumberOfFeatureMappingVariables());
  for(unsigned int i = 0; i < data.GetSize(); i++)
    data[i] = fd->GetFeaturedDesignElement(i)->GetPlainDesignValue();

  int ret = PyArg_ParseTuple(args, "O", &pyarray);
  PythonKernel::CheckPythonReturn(ret);
  data.Export((PyObject*) pyarray);
  Py_RETURN_NONE;
}



// grid stuff


/** return array with region names. The index is the 0-based region id  */
PyObject* grid_get_regions(PyObject *self, PyObject *args)
{
  StdVector<std::string> vec;
  domain->GetGrid()->GetRegionNames(vec);

  return PythonKernel::CreatePythonList(vec);
}


/** return number of nodes in mesh (for all regions) */
PyObject* grid_num_nodes(PyObject *self, PyObject *args)
{
  int regid = -1;
  int ret = PyArg_ParseTuple(args, "|i", &regid);
  PythonKernel::CheckPythonReturn(ret);

  RegionIdType reg = regid < 0 ? ALL_REGIONS : regid;

  return PyLong_FromLong(domain->GetGrid()->GetNumNodes(reg));
}

/** return region nodes for given region id to proper numpy array if dtype=int.
 * Use grid_get_regions() and grid_num_nodes() */
PyObject* grid_get_nodes(PyObject *self, PyObject *args)
{
  int regid;
  PyArrayObject* pyarray = nullptr;
  int ret = PyArg_ParseTuple(args, "iO", &regid, &pyarray);
  PythonKernel::CheckPythonReturn(ret);

  StdVector<unsigned int> nodes;
  domain->GetGrid()->GetNodesByRegion(nodes, regid);
  Vector<unsigned int> vec(nodes.GetSize(), nodes.GetPointer(), false); // wrap empty Vector around StdVector

  vec.Export((PyObject*) pyarray);
  Py_RETURN_NONE;
}



/** fill numpy array with coordinates of given node number (1-based).
 * arguments are node number and numpy array */
PyObject* grid_node_coord(PyObject *self, PyObject *args)
{
  int num;
  PyArrayObject* pycoord = nullptr;

  int ret = PyArg_ParseTuple(args, "iO", &num, &pycoord);
  PythonKernel::CheckPythonReturn(ret);

  Vector<double> vec; // how sad to have this temporary objects :(
  domain->GetGrid()->GetNodeCoordinate(vec, num); // not updated
  vec.Export((PyObject*) pycoord);
  Py_RETURN_NONE;
}

/** fill numpy array with all nodal coordinates */
PyObject* grid_all_nodes_coord(PyObject *self, PyObject *args)
{
  PyArrayObject* numpy = nullptr;
  int ret = PyArg_ParseTuple(args, "O", &numpy);
  PythonKernel::CheckPythonReturn(ret);

  unsigned int rows = PyArray_DIM(numpy,0);
  unsigned int cols = PyArray_DIM(numpy,1);

  GridCFS* grid = dynamic_cast<GridCFS*>(domain->GetGrid());
  const StdVector<Vector<double> >& data = grid->GetNodeCoordinates();

  if(data.GetSize() != rows || data[0].GetSize() != cols)
    EXCEPTION("numpy array is expected to be size (" << data.GetSize() << ", " << data[0].GetSize() << ") but is (" << rows << "," << cols <<")");

  for(unsigned int r = 0; r < rows; r++)
    for(unsigned int c = 0; c < cols; c++)
      *((double*) PyArray_GETPTR2(numpy,r,c)) = data[r][c];

  Py_RETURN_NONE;
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

/** return list of tuples (solution type, solution name) by singe pde idx */
PyObject* single_pde_solutions(PyObject *self, PyObject *args)
{
  const StdVector<SinglePDE*>& pdes = domain->GetSinglePDEs();
  int pde;
  int ret = PyArg_ParseTuple(args, "i", &pde);
  PythonKernel::CheckPythonReturn(ret);
  if(pde < 0 || pde >= (int) pdes.GetSize())
    throw Exception("invalid pde idx for single_pde_solutions(), needs to be below " + std::to_string(domain->GetSinglePDEs().GetSize()));

  PyObject* list = PyList_New(pdes[pde]->GetFeFunctions().size());
  int cnt = 0;
  for(const auto& it : pdes[pde]->GetFeFunctions())
  {
    PyObject* tuple = PyTuple_New(2);
    PyTuple_SetItem(tuple, 0, PyLong_FromLong(it.first));
    PyTuple_SetItem(tuple, 1, PyUnicode_FromString(SolutionTypeEnum.ToString(it.first).c_str()));
    PyList_SetItem(list, cnt,tuple);
    cnt++;
  }

  return list;
}

/** return equation mapping for single pde's solution.
 * Attributes are number of single pde, numpy array of node numbers (dtype=int) and numpy array of indices or size times dofs (dtype=int).
 * Enhance for attribute with solution type of the pde has more solutions. */
PyObject* single_pde_mapping(PyObject *self, PyObject *args)
{
  int pde;
  PyArrayObject* pynodenr = nullptr;
  PyArrayObject* pyeqns = nullptr;

  int ret = PyArg_ParseTuple(args, "iOO", &pde, &pynodenr, &pyeqns);
  PythonKernel::CheckPythonReturn(ret);

  const StdVector<SinglePDE*>& pdes = domain->GetSinglePDEs();
  if(pdes[pde]->GetFeFunctions().size() != 1)
    throw Exception("enhance single_pde_mapping() for solution type. Available for " + pdes[pde]->GetName()
                    + " are " + std::to_string(pdes[pde]->GetFeFunctions().size()));

  BaseFeFunction* fe = pdes[pde]->GetFeFunctions().begin()->second.get();
  FeSpace* space = fe->GetFeSpace().get();
  const FeSpace::SingleEqnMap& map = space->GetNodeMap();
  Vector<int> nodenr(map.eqns.size());
  int dofs = map.eqns.begin()->second.GetSize();
  Vector<int> eqns(nodenr.GetSize() * dofs);
  unsigned int cnt=0;
  for(const auto& itr : map.eqns)
  {
    nodenr[cnt] = itr.first;
    for(unsigned int d=0; d < itr.second.GetSize(); d++)
      eqns[cnt*dofs + d] = itr.second[d];
    cnt++;
  }

  nodenr.Export((PyObject*) pynodenr);
  eqns.Export((PyObject*) pyeqns);

  Py_RETURN_NONE;
}

/** return list with FeFunction names names */
PyObject* fe_function_names(PyObject *self, PyObject *args)
{
  assert(domain != NULL);

  StdVector<std::string> vec = domain->GetBasePDE()->GetAssemble()->GetAlgSys()->GetFeFunctionsNames();
  return PythonKernel::CreatePythonList(vec);
}

/** return total number of equations. Don't call to early!!! */
PyObject* fe_function_total_equations(PyObject *self, PyObject *args)
{
  const StdVector<unsigned int>& eqn = domain->GetBasePDE()->GetAssemble()->GetAlgSys()->GetFeFunctionsTotalEquations();
  StdVector<int> vec(eqn.GetSize());
  std::copy(eqn.begin(), eqn.end(), vec.begin());

  return PythonKernel::CreatePythonList(vec);
}

/** copy solution vector by given FeFunction
 * Argumensts are idx of fe function, bool for setIDBC, numpy array with proper size to be filled. */
PyObject* fe_function_solution(PyObject *self, PyObject *args)
{
  int fe;
  bool idbc;
  PyArrayObject *array = nullptr;
  int ret = PyArg_ParseTuple(args, "ipO", &fe, &idbc, &array);
  PythonKernel::CheckPythonReturn(ret);
  PythonKernel::CheckPythonReturn((PyObject*) array);

  // up to now we assume only real -> extend in case
  Vector<double> sol;
  domain->GetBasePDE()->GetAssemble()->GetAlgSys()->GetSolutionVal(sol, fe, idbc);
  sol.Export((PyObject*) array);
  Py_RETURN_NONE;
}

/** Nodal fe-function solution for given pde and solution type as vector of nodal solution
 * The first parameter is the 1-based node number, the pde-idx is usually 0, the solution type can be obtained by single_pde_solutions() */
PyObject* fe_function_nodal_solution(PyObject *self, PyObject *args)
{
  int pde_idx;
  int sol_type;
  int node;
  int ret = PyArg_ParseTuple(args, "iii", &node, &pde_idx, &sol_type); // to be extendend for default types!
  PythonKernel::CheckPythonReturn(ret);

  NodeList nodeList(domain->GetGrid());
  StdVector<UInt> nodeId(1);
  nodeId[0] = node;
  nodeList.SetNodes(nodeId);
  Vector<double> stateSol; // we get one scalar

  if(pde_idx < 0 || pde_idx >= (int) domain->GetSinglePDEs().GetSize())
    EXCEPTION("invalid pde-idx " << pde_idx << " given. See single_pde_names()");
  SinglePDE* pde =  domain->GetSinglePDEs()[pde_idx];

  if(pde->GetFeFunctions().count((SolutionType) sol_type) == 0)
    EXCEPTION("invalid solution type " << sol_type << " given. See single_pde_solutions()");
  BaseFeFunction* fe = pde->GetFeFunction((SolutionType) sol_type).get();

  fe->GetEntitySolution(stateSol, nodeList.GetIterator());

  return PythonKernel::CreatePythonList(stateSol);
}



// general stuff

/** return the problem name / simulation name of the current cfs run - the name of the output files */
PyObject* get_simulation_name(PyObject *self, PyObject *args)
{
  std::string name = progOpts->GetSimName();
  return PyUnicode_FromString(name.c_str());
}

PyMethodDef PythonKernel::cfs_methods[] =
{
  /** general cfs stuff */
  {"get_simulation_name", get_simulation_name, METH_VARARGS, "return the simulation=problem name"},

  /** grid stuff */
  {"grid_get_regions", grid_get_regions, METH_VARARGS, "return list of region names where the index is the 0-based region-id"},
  {"grid_num_nodes", grid_num_nodes, METH_VARARGS, "return number of nodes in region (if 0-based region-id is given) or for all regions"},
  {"grid_get_nodes", grid_get_nodes, METH_VARARGS, "fill 0-based region-id the numpy array of type int with 1-based node numbers"},
  {"grid_all_nodes_coord", grid_all_nodes_coord, METH_VARARGS, "return all nodes to a numpy array (num_nodes, dim). Note that node numbers are 1-based and here stored 0-based"},
  {"grid_node_coord", grid_node_coord, METH_VARARGS, "arguments are node number (1-based) and numpy array with proper size to take coordinates. Without deformation"},

  /** pde stuff */
  {"single_pde_names", single_pde_names, METH_VARARGS, "return a list with the names of the single pdes"},
  {"single_pde_solutions", single_pde_solutions, METH_VARARGS, "return a list of tuples (solution type, solution name) of the given single pde idx (usually 0)"},
  {"single_pde_mapping", single_pde_mapping, METH_VARARGS, "return mapping for indexed single pde (1. parameter). Provide numpy array (dtype=int) for nodes and one for indices which has size times dof"},
  {"fe_function_names", fe_function_names, METH_VARARGS, "return a list with the FeFunction names at their 0-based idx"},
  {"fe_function_total_equations", fe_function_total_equations, METH_VARARGS, "return total number of equations. Empty if called too early"},
  {"fe_function_solution", fe_function_solution, METH_VARARGS, "give real FeFunction solution to numpy array of proper size. Arguments: fe_idx (int), idbc (bool), sol (numpy array). See fe_function_total_equations()"},
  {"fe_function_nodal_solution", fe_function_nodal_solution, METH_VARARGS, "Give nodal solution vector by node (1-based), pde-idx (0-based), solution type (-> single_pde_solutions()) "},

  /* python mesh reader */
  {"set_nodes", PythonKernel::mesher_set_nodes, METH_VARARGS, "Set node coordinates via numpy array with 2 or 3 columns"},
  {"set_regions", PythonKernel::mesher_set_regions, METH_VARARGS, "set region names by list of strings"},
  {"add_elements", PythonKernel::mesher_add_elements, METH_VARARGS, "parameters total_number, fe_type, int numpy array of 1-based node ids"},
  {"add_named_nodes", PythonKernel::mesher_add_named_nodes, METH_VARARGS, "string and int numpy array of 1-based node ids"},
  {"add_naned_elements", PythonKernel::mesher_add_named_elements, METH_VARARGS, "string and int numpy array of 1-based element ids"},

  /* general optimization */
  {"getDims", opt_getDims, METH_VARARGS, "Returns info on optimization design domain dimensions: dim, nx, ny, nz."},
  {"opt_stop", opt_stop, METH_VARARGS, "make optimization to stop after current iteration give bool if converged and string with message"},
  {"opt_register_log_property", opt_register_log_property, METH_VARARGS, "register a property for file and info.xml iteration log. Give name:string and initial:string. Update via opt_set_log_property()"},
  {"opt_set_log_property", opt_set_log_property, METH_VARARGS, "update file and info.xml iteration propery via name:string and value:string. register first!"},
  {"optimizer_get_properties", optimizer_get_properties, METH_VARARGS, "get properties / infos about the actual optimizer used in optimization as string/string dict"},
  {"optimizer_set_property", optimizer_set_property, METH_VARARGS, "set name:string and value:string for the current optimizer, only a few implement this."},
  {"get_opt_design_size", get_opt_design_size, METH_VARARGS, "Returns number of design variables"},
  {"get_opt_design_value", get_opt_design_value, METH_VARARGS, "Give single DesignSpace::GetDesignValue() for 0-based index with optional access (default plain)"},
  {"get_opt_design_values", get_opt_design_values, METH_VARARGS, "Call DesignElement::GetValues() with attributes numpy array and optional access"},
  {"get_opt_stopping_rules", get_opt_stopping_rules, METH_VARARGS, "return the stoppping rules"},
  {"get_opt_iteraton", get_opt_iteraton, METH_VARARGS, "Return Optimization->GetCurrentIteration()"},
  {"get_opt_function_values", get_opt_function_values, METH_VARARGS, "return two string/string maps with name an value for objective functions and constraints as in .info.xml"},
  {"get_opt_function_properties", get_opt_function_properties, METH_VARARGS, "return string/string maps with with property of objective/constraint/observation"},
  {"get_opt_filter_values", get_opt_filter_values, METH_VARARGS, "return string/string map for given filter index with filter properties. Check total_filter property!"},
  {"set_opt_filter_values", set_opt_filter_values, METH_VARARGS, "set optimization filter idx, beta, eta, scale, offset where idx and beta is required. No keywords allowed!"},

  /* feature mapping optimization design */
  {"feature_mapping_num_parameters", feature_mapping_num_parameters, METH_VARARGS, "Return the number of feature mapping parameters including fixed variables"},
  {"feature_mapping_get_parameters", feature_mapping_get_parameters, METH_VARARGS, "Fills given numpy array with feature mapping parameters"},

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
