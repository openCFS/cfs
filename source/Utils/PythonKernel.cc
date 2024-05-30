#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "Utils/PythonKernel.hh"
#include "DataInOut/ProgramOptions.hh"
#include "General/Exception.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/SimInOut/python/SimInputPython.hh"
#include "Optimization/Optimizer/PythonOptimizer.hh"


// declare class specific logging stream
DEFINE_LOG(pykernel, "pykernel")


namespace CoupledField
{
using std::string;
using std::to_string;
using std::pair;

// the external definition of PythonKernel* python is in tools.cc as we don't compile and link this file w/o USE_EMBEDDED_PYTHON

// set static attributes
SimInputPython*  PythonKernel::mesher_ = NULL;
PythonOptimizer* PythonKernel::pyopt_  = NULL;


Enum<PythonKernel::Hook> PythonKernel::hook;

/** actual construtor */
void PythonKernel::Init(PtrParamNode pn, PtrParamNode info)
{
  info_ = info->Get(ParamNode::HEADER); // do info->Get("python") only when we have something to show
  python = this; // set the global pointer to this object

  hook.SetName("PythonKernel::Hook");
  hook.Add(ASSEMBLE_RHS, "assemble_rhs");
  hook.Add(OPT_EVAL_FUNC, "opt_eval_func");
  hook.Add(OPT_EVAL_GRAD, "opt_eval_grad");
  hook.Add(OPT_POST_INIT, "opt_post_init");
  hook.Add(OPT_POST_ITER, "opt_post_iter");
  hook.Add(POST_DOMAIN_INIT, "post_domain_init");
  hook.Add(POST_GRID, "post_grid");
  hook.Add(POST_SOLVE_PROBLEM, "post_solve_problem");

  if(pn)
  {
    PtrParamNode pi = info_->Get("python");
    LoadStatus status = LoadPythonModule(pn);
    kernel_ = status.module;
    pi->Get("file")->SetValue(status.full_file);
    if(progOpts->DoDetailedInfo())
      pi->Get("syspath")->SetValue(status.sys_path.ToString(TS_PLAIN, ":"));


    ParamNodeList lst = pn->GetList("callback");
    for(auto callback : lst) {
      string name = callback->Get("hook")->As<string>();
      Hook key = hook.Parse(name);
      if(hooks_.find(key) != hooks_.end())
        throw Exception("python callback hook " + name + " used multiple times");

      HookParam hp;
      hp.function = callback->Get("function")->As<string>();;
      hp.options = ParseOptions(callback->GetList("option"));

      hooks_[key] = hp;

      PtrParamNode cbpn = pi->Get("callback", ParamNode::APPEND);
      cbpn->Get("hook")->SetValue(name);
      cbpn->Get("function")->SetValue(hp.function);
      if(!hp.options.IsEmpty())
        cbpn->Get("options")->SetValue(callback->GetList("option"));
    }
  }
}

/** actual destructor */
void PythonKernel::Release()
{
  if(initialized_)
  {
    Py_XDECREF(kernel_); // can be NULL
    kernel_ = NULL;
    Py_Finalize();
  }
  python = NULL; // reset the global object pointer
}

void PythonKernel::InitInterpreter()
{
  assert(!initialized_); // run once!

  PyImport_AppendInittab("cfs", SetModulFunctions);
  Py_Initialize();

  PyObject* version = PySys_GetObject("version");
  if(!version)
    throw Exception("embedded Python not working: " + PyErr());

  const char* c_str = PyUnicode_AsUTF8(version);
  LOG_DBG(pykernel) << "II: org version=" << c_str;
  assert(c_str);
  string str(c_str);
  boost::replace_all(str, "\n"," ");
  info_->Get("python/version")->SetValue(str);
  Py_XDECREF(version);

  initialized_ = true;
}

PythonKernel::LoadStatus PythonKernel::LoadPythonModule(const string& file, const string opt_path)
{
  if(!initialized_)
    InitInterpreter(); // we assume no critical threadding issue

  LoadStatus ret;

  // our share/python
  boost::filesystem::path share = progOpts->GetSchemaPath().parent_path().append("python");

  boost::filesystem::path givenname(file); // default
  if(opt_path != "")
  {
    if(opt_path == "cfs:share:python")
      givenname = share.append(file);
    else
      givenname = boost::filesystem::path(opt_path).append(file);
  }
  ret.full_file = givenname.string();

  if(!boost::filesystem::exists(givenname))
    throw Exception("cannot find python file '" + givenname.string() + "' for python optimizer");
  boost::filesystem::path absolute = boost::filesystem::absolute(givenname.parent_path());

  // add to PYTHONPATH to make it more realistic to load the module :)
  PyObject* sysPath = PySys_GetObject((char*) "path"); // must not decref after append - note the difference to the attribute syspath!
  LOG_DBG(pykernel) << "IPM: original sysPath=" << ConvertPythonList<string>(sysPath).ToString();
  // first add share
  PyList_Insert(sysPath, 0, PyUnicode_FromString(share.string().c_str()));
  // then add as second our module (might be doubled)
  PyList_Insert(sysPath, 1, PyUnicode_FromString(absolute.string().c_str()));
  LOG_DBG(pykernel) << "IPM: final sysPath=" << ConvertPythonList<string>(sysPath).ToString();
  // return the new sysPath
  ConvertPythonList<string>(ret.sys_path, sysPath); // again copy constructor magic :)

  // now load the module
  boost::filesystem::path name = boost::filesystem::change_extension(givenname.filename(), "");
  ret.module = PyImport_ImportModule(name.string().c_str());

  if(!ret.module)
  {
    throw Exception("Cannot load '" + givenname.string() + "' as embbeded python module: " + PyErr());
  }

  return ret;
}

PythonKernel::LoadStatus PythonKernel::LoadPythonModule(PtrParamNode pn)
{
  if(pn->Has("path"))
    return LoadPythonModule(pn->Get("file")->As<string>(), pn->Get("path")->As<string>());
  else
    return LoadPythonModule(pn->Get("file")->As<string>());
}



void PythonKernel::CallHook(Hook key)
{
  const auto& iter = python->hooks_.find(key);

  if(iter == python->hooks_.end())
    return;

  const HookParam& param = iter->second;

  assert(param.function != "");
  assert(python->kernel_ != NULL);

  PyObject* func = PyObject_GetAttrString(python->kernel_, param.function.c_str());
  CheckPythonFunction(func, param.function.c_str());

  PyObject* arg = PyTuple_New(1);
  PyTuple_SetItem(arg, 0, CreatePythonDict(param.options)); // PyTuple_SetItem() steals the reference

  PyObject* call = PyObject_CallObject(func, arg);
  CheckPythonReturn(call);

  Py_XDECREF(call);
  Py_XDECREF(arg);
  Py_XDECREF(func);
}

bool PythonKernel::CheckMesher()
{
  if(PythonKernel::mesher_)
    return true;

  PyErr_SetString(PyExc_ValueError, "the python mesh input did not (yet) register itself");
  return false;
}

bool PythonKernel::CheckPyOpt()
{
  if(PythonKernel::pyopt_)
    return true;

  PyErr_SetString(PyExc_ValueError, "the python optimizer did not (yet) register itself");
  return false;
}


bool PythonKernel::CheckOpt()
{
  if(domain->GetOptimization())
    return true;

  PyErr_SetString(PyExc_ValueError, "no optimization problem given");
  return false;
}

/** try to case the optimization DesignSpace to FeatureMapping or throw an error */
FeaturedDesign* PythonKernel::GetFeaturedDesign()
{
  FeaturedDesign* fd = dynamic_cast<FeaturedDesign*>(domain->GetDesign());
  if(fd == nullptr)
    throw Exception("no FeaturedDesign defined in openCFS problem");
  return fd;
}



template<class TYPE>
void PythonKernel::ConvertPythonList(Container<TYPE>& vec, PyObject* list)
{
  if (!PyList_Check(list))
    throw Exception("provided object is no list");

  int size = PyList_Size(list);

  vec.Resize(size);

  for(int i = 0; i < size; i++)
  {
    PyObject* item = PyList_GetItem(list, i); // borrowed reference
    vec[i] = Convert<TYPE>(item);
  }
}


template<>
std::string PythonKernel::Convert<std::string>(PyObject* item)
{
  const char* c = PyUnicode_AsUTF8(item);
  if(!c)
    EXCEPTION("cannot convert to string: " << ToString(item));
  return string(c);
}

template<>
double PythonKernel::Convert<double>(PyObject* item)
{
  return PyFloat_AsDouble(item);
}

template<>
Vector<int> PythonKernel::Convert<Vector<int> >(PyObject* item)
{
  Vector<int> vec(item, false); // don't decref
  return vec;
}


template<class TYPE>
PyObject* PythonKernel::CreatePythonList(const Container<TYPE>& vec)
{
  PyObject* list = PyList_New(vec.GetSize());
  for(unsigned int i = 0; i < vec.GetSize(); i++)
    PyList_SetItem(list, i, Create<TYPE>(vec[i]));
  return list;
}


template<>
PyObject* PythonKernel::Create<std::string>(const std::string& val)
{
  return PyUnicode_FromString(val.c_str());
}

template<>
PyObject* PythonKernel::Create<double>(const double& val)
{
  return PyFloat_FromDouble(val);
}

template<>
PyObject* PythonKernel::Create<int>(const int& val)
{
  return PyLong_FromLong(val);
}


PyObject* PythonKernel::CreatePythonDict(const StdVector<pair<string, string> > options)
{
  // test dictionary
  PyObject* dict = PyDict_New();
  for(auto opt : options)
    PyDict_SetItemString(dict, opt.first.c_str(), PyUnicode_FromString(opt.second.c_str()));

  return dict;
}

StdVector<pair<string, string> > PythonKernel::ParseOptions(ParamNodeList lst)
{
  StdVector<pair<string, string> > options;
  options.Reserve(lst.GetSize());
  for(auto opt : lst)
    options.Push_back(std::make_pair(opt->Get("key")->As<string>(),opt->Get("value")->As<string>()));

  return options; // believe in smart return value optimization
}



StdVector<PyObject*> PythonKernel::ParseNumpyArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref)
{
  if(PyTuple_Size(args) != expect)
    throw Exception("expects " + to_string(expect) + " parameters in python function call, got " + to_string(PyTuple_Size(args)));

  StdVector<PyObject*> obj(expect);

  // could clearly be done smarter!
  switch(expect)
  {
  case 1:
    PyArg_ParseTuple(args, "O!", &PyArray_Type, &obj[0]);
    break;

  case 2:
    PyArg_ParseTuple(args, "O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1]);
    break;

  case 3:
    PyArg_ParseTuple(args, "O!|O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1], &PyArray_Type, &obj[2]);
    break;

  case 4:
    PyArg_ParseTuple(args, "O!|O!|O!|O!", &PyArray_Type, &obj[0], &PyArray_Type, &obj[1], &PyArray_Type, &obj[2],  &PyArray_Type, &obj[3]);

    break;
  default:
    EXCEPTION("not implemented size " << expect);
  }
  data.Resize(expect);

  for(unsigned int i = 0; i < obj.GetSize(); i++)
    data[i].Fill(obj[i], decref);

  return obj;
}


void PythonKernel::CheckPythonReturn(PyObject* pyobject, const char* name)
{
  if(!pyobject)
    // allow for stacktrace printing on console
    throw Exception("Python returns error " + string(name != NULL ? name : "") + ": " + PyErr());
}

void PythonKernel::CheckPythonReturn(int ret, const char* name)
{
  if(ret == 0)
    throw Exception("Python returns error " + string(name != NULL ? name : "") + ": " + PyErr());
}


void PythonKernel::CheckPythonFunction(PyObject* pyobject, const char* name)
{
  if(!(pyobject && PyCallable_Check(pyobject)))
    throw Exception("Cannot call Python function " + string(name ? name : ""));
}


string PythonKernel::PyErr(bool call_PyErr_Print)
{
  // prints something like
  // ---
  // Traceback (most recent call last):
  //  File "/Users/fwein/code/cfs/debug/../source/unittests/embeddedpython.py", line 41, in get_vec
  //    print('embeddedpython.py get_vec ->', v,n)
  // NameError: name 'v' is not defined
  if(call_PyErr_Print)
    PyErr_Print();

  // https://docs.python.org/3/library/sys.html
  PyObject* ltv = PySys_GetObject((char*) "last_value"); // borrowed reference, don't decref

  return ToString(ltv);
}

std::string PythonKernel::ToString(PyObject* obj)
{
  if(!obj)
    return string();

  PyObject* so = PyObject_Repr(obj); // new reference!
  assert(so);
  const char* sc = PyUnicode_AsUTF8(so);
  assert(sc != NULL);
  string ret(sc);
  Py_DECREF(so);
  return ret;
}


template<class T>
Matrix<T> PythonKernel::Numpy2DArrayToMatrix(PyObject* numpy, Matrix<T>& out, bool decref)
{
  if(!PyArray_Check(numpy))
    throw Exception("provided PyObject seems not to be a numpy array");

  PyArrayObject* arr = (PyArrayObject*) numpy;

  unsigned int rows = PyArray_DIM(arr,0);
  unsigned int cols = PyArray_DIM(arr,1);

  out.Resize(rows,cols);

  for(unsigned int r = 0; r < rows; r++)
    for(unsigned int c = 0; c < cols; c++)
      out[r][c] = *((T*) PyArray_GETPTR2(arr,r,c) );

  if(decref)
     Py_DECREF(numpy);

  return out;
}

template<class T>
void PythonKernel::MatrixToNumpyArray(const Matrix<T>& in, PyObject* numpy)
{
  if(!PyArray_Check(numpy))
    throw Exception("provided PyObject seems not to be a numpy array");

  PyArrayObject* arr = (PyArrayObject*) numpy;

  unsigned int rows = PyArray_DIM(arr,0);
  unsigned int cols = PyArray_DIM(arr,1);

  if(in.GetNumRows() != rows || in.GetNumCols() != cols)
    EXCEPTION("numpy array is expected to be size (" << in.GetNumRows() << ", " << in.GetNumCols() << ") but is (" << rows << "," << cols <<")");

  for(unsigned int r = 0; r < rows; r++)
    for(unsigned int c = 0; c < cols; c++)
      *((T*) PyArray_GETPTR2(arr,r,c)) = in[r][c];
}

// the pyton function stuff is in PythonKernelFunction.cc

template Matrix<double> PythonKernel::Numpy2DArrayToMatrix<double>(PyObject*, Matrix<double>&, bool);
template Matrix<int> PythonKernel::Numpy2DArrayToMatrix<int>(PyObject*, Matrix<int>&, bool);

template void PythonKernel::MatrixToNumpyArray<double>(const Matrix<double>&, PyObject*);
template void PythonKernel::MatrixToNumpyArray<int>(const Matrix<int>&, PyObject*);

template PyObject* PythonKernel::CreatePythonList<std::string>(const Container<std::string>&);
template PyObject* PythonKernel::CreatePythonList<double>(const Container<double>&);
template PyObject* PythonKernel::CreatePythonList<int>(const Container<int>&);

template void PythonKernel::ConvertPythonList<std::string>(Container<std::string>&, PyObject*);
template void PythonKernel::ConvertPythonList<double>(Container<double>&, PyObject*);
template void PythonKernel::ConvertPythonList<Vector<int> >(Container<Vector<int> >&, PyObject*);


} // end of namespace
