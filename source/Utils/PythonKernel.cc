#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
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
namespace fs = boost::filesystem;

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
  if(Py_IsInitialized())
  {
    Py_XDECREF(kernel_); // can be NULL
    kernel_ = nullptr;
    Py_Finalize(); // note that we may not re-initialize the interpreter as numpy would not load
  }
  python = nullptr; // reset the global object pointer
}

void PythonKernel::InitInterpreter()
{
  assert(!Py_IsInitialized()); // run once!

  char* pp = std::getenv("PYTHONPATH");
  char* ph = std::getenv("PYTHONHOME");
  LOG_DBG(pykernel) << "II: PYTHONPATH=" << (pp ? string(pp) : "NULL") << " PYTHONHOME=" << (ph ? string(ph) : "NULL");
  
  PyImport_AppendInittab("cfs", &SetModulFunctions);

  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  // here config could be modified

  LOG_DBG(pykernel) << "II: home=" << (config.home ? "set" : "NULL");
 
  PyStatus status = Py_InitializeFromConfig(&config);
  assert(!PyStatus_IsError(status));
  if(PyStatus_Exception(status)) // checks the status struct if exist reason was exception
  {
    PyConfig_Clear(&config);
    throw Exception("error initializing python interpreter: " + string(status.err_msg));
  }
  PyConfig_Clear(&config);

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
}

PythonKernel::LoadStatus PythonKernel::LoadPythonModule(const string& file, const string opt_path, bool preload_numpy)
{
  LOG_DBG(pykernel) << "LPM: f=" << file << "op=" << opt_path << " i=" << Py_IsInitialized();

  if(!Py_IsInitialized())
    InitInterpreter(); // we assume no critical threadding issue

  assert(Py_IsInitialized());

  LoadStatus ls;

  // file (e.g. make_mesh.py) is loaded as module (more complicated than running a file)
  // for this the extension will be removed and the module will be imported by it's name in all sys.path

  // we add the path to sys.path. I empty this is cwd, otherwise the special c:s:p or a given path

  // we check the existence of the file manually to get a better error message than C-API
  // to this end, we construct the absolute path as it might contain cf:share:python

  fs::path share = progOpts->GetSchemaPath().parent_path().append("python");  // /home/fwein/code/cfs/share/xml -> share/python
  fs::path givenname(file); // default
  if(opt_path != "")
    givenname = opt_path == "cfs:share:python" ? share / file : fs::path(opt_path) / file;
  ls.full_file = givenname.string();
  if(!boost::filesystem::exists(givenname))
    throw Exception("cannot find python file '" + givenname.string() + "' for python optimizer");
  LOG_DBG(pykernel) << "LPM: gn=" << givenname.string();

  // determine the path to be added to sys.path
  fs::path fs_path;
  if(opt_path == "" || opt_path == ".")
    fs_path = fs::current_path();
  else
    fs_path = opt_path == "cfs:share:python" ? share : fs::path(opt_path);
  LOG_DBG(pykernel) << "LPM: fs_path=" << fs_path.string() << " share=" << share.string();

  // up to now we did not change anything
  // with Python 3.12 vs 3.11 we have issues loading numpy (segfault in many test cases) and
  // this seams to to help. Ususally we need numpy anyway.
  #ifndef NDEBUG
    bool verb = true;
  #else
    bool verb = progOpts->DoDetailedInfo();
  #endif
  // PyRun_SimpleString("import sys;print('python sys.path', sys.path, flush=True)");
  // note that you need for virtual python what you see manual via sysconfig.get_path('platlib')
  // e.g. /Users/fwein/python3.14/lib/python3.14/site-packages
  // however, sysconfig.get_path('platlib') in embedded python shows the non-virtual base.
  // best is to set the path first in PYTHONPATH when you activate your virtual python environment
  if(preload_numpy)
  {
    if(verb) PyRun_SimpleString("print('++ pre-load numpy ... ',flush=True, end='')");
    PyRun_SimpleString("import numpy as np");
    if(verb) PyRun_SimpleString("print(np.__version__,flush=True)");
  }
  ls.preloaded_numpy = preload_numpy;

  // we insert first to sys.path such that we "overwrite" files in the search path
  // note that sys.path contains also PYTHONPATH
  PyObject* sysPath = PySys_GetObject((char*) "path"); // must not decref after append - note the difference to the attribute syspath!
  LOG_DBG(pykernel) << "LPM: original sysPath=" << ConvertPythonList<string>(sysPath).ToString();
  // insert path
  PyList_Insert(sysPath, 0, PyUnicode_FromString(fs_path.string().c_str()));
  // add share/python to have the cfs stuff available
  if(opt_path != "cfs:share:python")
    PyList_Append(sysPath, PyUnicode_FromString(share.string().c_str()));
  LOG_DBG(pykernel) << "LPM: final sysPath=" << ConvertPythonList<string>(sysPath).ToString();
  // return the new sysPath
  ConvertPythonList<string>(ls.sys_path, sysPath); // again copy constructor magic :)
  
  // PyRun_SimpleString("import sys;print('python sys.path', sys.path, flush=True)");

  // now load the module which is the python file without extension
  fs::path name = fs::change_extension(givenname.filename(), "");
  LOG_DBG(pykernel) << "LPM: name=" << name.string();

  ls.module = PyImport_ImportModule(name.string().c_str());

  if(!ls.module)
  {
    throw Exception("Cannot load '" + givenname.string() + "' as embbeded python module: " + PyErr());
  }

  return ls;
}

PythonKernel::LoadStatus PythonKernel::LoadPythonModule(PtrParamNode pn)
{
  LOG_DBG(pykernel) << "LPM: pn";
  bool pre = pn->Has("preload_numpy") ? pn->Get("preload_numpy")->As<bool>() : true;
  if(pn->Has("path"))
    return LoadPythonModule(pn->Get("file")->As<string>(), pn->Get("path")->As<string>(), pre);
  else
    return LoadPythonModule(pn->Get("file")->As<string>(), "", pre);
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
