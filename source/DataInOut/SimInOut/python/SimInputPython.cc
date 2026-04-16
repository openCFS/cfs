#include <string>
#include <iostream>
#include <algorithm>

#include <def_use_embedded_python.hh>
#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#endif


#define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
//#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include "SimInputPython.hh"
#include "DataInOut/ProgramOptions.hh"
#include "MatVec/Vector.hh"
#include "Utils/ToolsFull.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

// declare class specific logging stream
DEFINE_LOG(pymesh, "pymesh")

using std::string;
using std::to_string;

namespace CoupledField
{

SimInputPython::SimInputPython(std::string fileName, PtrParamNode inputNode, PtrParamNode infoNode) :
  SimInput(fileName, inputNode, infoNode)
{
  info_ = infoNode->Get(ParamNode::HEADER)->Get("domain/python");
  capabilities_.insert( SimInput::MESH);

  python->Register(this);
  // InitModule() does the real work
}

SimInputPython::~SimInputPython()
{
  python->Register(this, true); // remove
}


void SimInputPython::InitModule()
{
  PythonKernel::LoadStatus stat = python->LoadPythonModule(myParam_->Get("fileName")->As<string>(), myParam_->Get("path")->As<string>());

  module = stat.module;
  givenname = stat.full_file;

  info_->Get("file")->SetValue(givenname);
  info_->Get("function")->SetValue("set_cfs_mesh");
  if(progOpts->DoDetailedInfo())
    info_->Get("syspath")->SetValue(stat.sys_path.ToString(TS_PLAIN, ":"));

  // the options are given to the python function set_cfs_mesh()
  ParamNodeList lst = myParam_->GetList("option");
  options = PythonKernel::ParseOptions(lst);

  info_->Get("options")->SetValue(lst);
  // to be continued in ReadMesh()

  import_array1();
}


void SimInputPython::ReadMesh(Grid* grid)
{
  assert(grid != NULL);
  this->grid = grid;
  // call set_cfs_mesh(options)
  // this python function shall call set_nodes, set_regions, add_elements and optionally add_named_nodes and add_named_elements

  assert(grid->GetNumNodes() == 0);
  assert(grid->regionData.IsEmpty());

  PyObject* init = PyObject_GetAttrString(module, "set_cfs_mesh");
  if(init && PyCallable_Check(init)) {
    PyObject* arg = PyTuple_New(1);
    PyTuple_SetItem(arg, 0, PythonKernel::CreatePythonDict(options));
    LOG_DBG(pymesh) << "RM: call set_cfs_mesh()";
    PyObject* ret = PyObject_CallObject(init, arg);
    LOG_DBG(pymesh) << "RM: called set_cfs_mesh()";
    if(!ret)
      throw Exception("set_cfs_mesh() aborted in module " + givenname + ": " + PythonKernel::PyErr());
    Py_XDECREF(ret);
    Py_XDECREF(arg);
    Py_XDECREF(init);
  }
  else
    throw Exception("could not call set_cfs_mesh() in python module " + givenname);

  // to expensive too keep for cfs life time
  Py_XDECREF(module);
  module = NULL;
  // The interpreter stays running!

  if(grid->GetNumNodes() == 0 || total_elements_ == 0 || grid->regionData.IsEmpty())
    throw Exception("set_cfs_mesh() did not call (all of) set_nodes, set_regions and add_elements");

  if(grid->GetNumElems() != total_elements_)
    throw Exception("set_cfs_mesh() called add_elements with total_elements " + to_string(total_elements_) + " but ordered elements in mesh are " + to_string(grid->GetNumElems()));
}


void SimInputPython::SetNodes(PyObject* args)
{
  assert(grid != NULL);

  if(grid->GetNumNodes() > 0)
    throw Exception("set_nodes() called before");

  PyArrayObject *array = NULL;

  PyArg_ParseTuple(args, "O!", &PyArray_Type, &array);
  PythonKernel::CheckPythonReturn((PyObject*) array);

  if(!PyArray_Check(array))
    throw Exception("did not get numpy array with set_nodes()");

  unsigned int rows = PyArray_DIM(array,0);
  unsigned int cols = PyArray_DIM(array,1);

  if(!(rows > 0 && (cols == 2 || cols == 3)))
    EXCEPTION("numpy array for set_nodes has " << rows << " rows and " << cols << " cols, shall have 2 or 3 cols and positive rows");

  LOG_DBG(pymesh) << "SN: rows=" << rows << " cols=" << cols;

  grid->AddNodes(rows);

  Vector<double> loc(3);
  loc.Init(0.0); // third dimension keeps 0.0 in the 2d case
  for(unsigned int i = 0; i < rows; i++)
  {
    loc[0] = *((double*) PyArray_GETPTR2(array,i,0));
    loc[1] = *((double*) PyArray_GETPTR2(array,i,1));
    if(cols > 2)
      loc[2] = *((double*) PyArray_GETPTR2(array,i,2));
    LOG_DBG3(pymesh) << "SN: i=" << i << "->" << loc.ToString();
    grid->SetNodeCoordinate(i+1, loc);
  }

  // we must not decref array!
  Py_XDECREF(args);
}

void SimInputPython::SetRegions(PyObject* args)
{
  if(grid->regionData.GetSize() > 0)
    throw Exception("set_regions() called before");

  // https://stackoverflow.com/questions/3253563/pass-list-as-argument-to-python-c-module
  PyObject* list;
  PyArg_ParseTuple(args, "O!", &PyList_Type, &list);
  PythonKernel::CheckPythonReturn((PyObject*) list);

  unsigned int numLines = PyList_Size(list);
  if(numLines < 1)
    throw Exception("cfs.set_regions() not called with a valid list");

  StdVector<string>       reg_names(numLines);
  StdVector<RegionIdType> reg_ids(numLines); // = int

  for(unsigned int i = 0; i < numLines; i++) {
    PyObject* obj = PyList_GetItem(list, i);
    const char* line = PyUnicode_AsUTF8(obj);
    reg_names[i] = string(line);
    reg_ids[i] = i;
  }
  LOG_DBG(pymesh) << "SR: regs=" << reg_names.ToString();
  grid->AddRegions(reg_names, reg_ids);

  // we must not decref list as this would remove the item from python!
  Py_XDECREF(args);
}

void SimInputPython::AddElements(PyObject* args)
{
  // we get total, cfs_type, array of ints

  int total;
  int type;
  PyArrayObject *array = NULL;
  int ret = PyArg_ParseTuple(args, "iiO!", &total, &type, &PyArray_Type, &array);
  PythonKernel::CheckPythonReturn(ret);
  PythonKernel::CheckPythonReturn((PyObject*) array);
  LOG_DBG(pymesh) << "AE: total=" << total << " type=" << type << "=" << Elem::feType.ToString((Elem::FEType) type);


  if(total_elements_ == 0) {
    grid->AddElems(total);
    total_elements_ = total;
  }
  if(total_elements_ != (unsigned int) total)
    EXCEPTION("called cfs.add_elements() again with inconsistent total_elements: " << total_elements_ << " vs. " << total);

  // the mesh is never 1d. 1d elements we therefore count as 2d
  // note we can call this many times with 1d/2d elements last.
  bool this_is_3d = type >= Elem::ET_HEXA8;
  dim_ = std::max((int) dim_, this_is_3d ? 3 : 2);

  unsigned int rows = PyArray_DIM(array,0);
  unsigned int cols = PyArray_DIM(array,1);
  Vector<unsigned int> connect(cols -2); // skip e-num and reg-id
  for(unsigned int r = 0; r < rows; r++)
  {
    int num =  *((int*) PyArray_GETPTR2(array,r,0));
    int rid =  *((int*) PyArray_GETPTR2(array,r,1));
    for(unsigned int c = 2; c < cols; c++)
      connect[c-2] = *((int*) PyArray_GETPTR2(array,r,c));
    // SetElemData(UInt ielem, Elem::FEType type, RegionIdType region, const UInt* connect)
    LOG_DBG3(pymesh) << "AE: r=" << r << " num=" << num << " rid=" << rid << " connect=" << connect.ToString();
    grid->SetElemData(num, (Elem::FEType) type, rid, connect.GetPointer());
  }

  Py_XDECREF(args);
}

void SimInputPython::AddNamedNodes(PyObject* args)
{
  AddNamedNodesElements(args, true);
}


void SimInputPython::AddNamedElements(PyObject* args)
{
  AddNamedNodesElements(args, false);
}

void SimInputPython::AddNamedNodesElements(PyObject* args, bool nodes)
{
  char* s;
  PyArrayObject *array = NULL;
  int ret = PyArg_ParseTuple(args, "sO!", &s, &PyArray_Type, &array);
  PythonKernel::CheckPythonReturn(ret);
  PythonKernel::CheckPythonReturn((PyObject*) array);
  LOG_DBG(pymesh) << "ANNE: action=" << (nodes ? "nodes" : "elements") << " s=" << s << " rows=" << PyArray_DIM(array,0);

  // copy numpy content without decrec
  Vector<unsigned int> vec((PyObject*) array, false); // shall be np.array(nodes, dtype=uintc)
  LOG_DBG3(pymesh) << "ANNE: vec=" << vec.ToString();
  // simple wrapper
  StdVector<unsigned int> stv(vec.GetPointer(), vec.GetSize()); 

  if(nodes)
    grid->AddNamedNodes(string(s), stv);
  else
    grid->AddNamedElems(string(s), stv);

  Py_XDECREF(args);
}


unsigned int SimInputPython::GetDim()
{
  return dim_;
}

} // end of namespace

