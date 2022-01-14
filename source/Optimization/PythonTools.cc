#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

#include "Optimization/PythonTools.hh"
#include "DataInOut/ProgramOptions.hh"
#include "General/Exception.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


// declare class specific logging stream
DEFINE_LOG(pytools, "pytools")


namespace CoupledField
{
using std::string;
using std::to_string;
using std::pair;


PyObject* InitializePythonModule(const string& file, const string& opt_path, PyObject* (*init_cfs)(), string* file_out, string* version_out, StdVector<std::string>* syspath)
{
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
  if(file_out)
    *file_out = givenname.string();

  if(!boost::filesystem::exists(givenname))
    throw Exception("cannot find python file '" + givenname.string() + "' for python optimizer");
  boost::filesystem::path absolute = boost::filesystem::absolute(givenname.parent_path());

  // needs to be done before Py_Initialize
  if(init_cfs != NULL)
    PyImport_AppendInittab("cfs", init_cfs);

  Py_Initialize();

  PyObject* version = PySys_GetObject("version");
  if(!version)
    throw Exception("embedded Python not working: " + PyErr());

  const char* c_str = PyUnicode_AsUTF8(version);
  LOG_DBG(pytools) << "IPM: org version=" << c_str;
  assert(c_str);
  string str(c_str);
  boost::replace_all(str, "\n"," ");
  if(version_out)
    *version_out = str;

  Py_XDECREF(version);

  // add to PYTHONPATH to make it more realistic to load the module :)
  PyObject* sysPath = PySys_GetObject((char*) "path"); // must not decref after append - note the difference to the attribute syspath!
  LOG_DBG(pytools) << "IPM: original sysPath=" << Convert(sysPath).ToString();
  // first add share
  PyList_Insert(sysPath, 0, PyUnicode_FromString(share.string().c_str()));
  // then add as second our module (might be doubled)
  PyList_Insert(sysPath, 1, PyUnicode_FromString(absolute.string().c_str()));
  LOG_DBG(pytools) << "IPM: final sysPath=" << Convert(sysPath).ToString();
  // return the new sysPath
  if(syspath)
    *syspath = Convert(sysPath); // again copy constructor magic :)


  // now load the module
  boost::filesystem::path name = boost::filesystem::change_extension(givenname.filename(), "");
  PyObject* module = PyImport_ImportModule(name.string().c_str());

  if(!module)
  {
    throw Exception("Cannot load '" + givenname.string() + "' as embbeded python module: " + PyErr());
  }

  return module;
}


StdVector<std::string> Convert(PyObject* list)
{
  if (!PyList_Check(list))
    throw Exception("provided object is no list");

  int size = PyList_Size(list);

  StdVector<string> res(size);
  for(int i = 0; i < size; i++)
  {
    PyObject* item = PyList_GetItem(list, i); // borrowed reference
    const char* c = PyUnicode_AsUTF8(item);
    if(!c)
      EXCEPTION("element with index " << i << " of python list is no string");
    res[i] = string(c);
  }

  return res;
}

StdVector<PyObject*> ParseNumpyArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref)
{
  if(PyTuple_GET_SIZE(args) != expect)
    throw Exception("expects " + to_string(expect) + " parameters in python function call, got " + to_string(PyTuple_GET_SIZE(args)));

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


StdVector<pair<string, string> > ParseOptions(ParamNodeList lst)
{
  StdVector<pair<string, string> > options;
  options.Reserve(lst.GetSize());
  for(auto opt : lst)
    options.Push_back(std::make_pair(opt->Get("key")->As<string>(),opt->Get("value")->As<string>()));

  return options; // believe in smart return value optimization
}

PyObject* CreatePythonDict(const StdVector<pair<string, string> > options)
{
  // test dictionary
  PyObject* dict = PyDict_New();
  for(auto opt : options)
    PyDict_SetItemString(dict, opt.first.c_str(), PyUnicode_FromString(opt.second.c_str()));

  return dict;
}

void CheckPythonReturn(PyObject* pyobject, const char* name)
{
  if(!pyobject)
    // allow for stacktrace printing on console
    throw Exception("Python returns error " + string(name != NULL ? name : "") + ": " + PyErr());
}

void CheckPythonReturn(int ret, const char* name)
{
  if(ret == 0)
    throw Exception("Python returns error " + string(name != NULL ? name : "") + ": " + PyErr());
}


void CheckPythonFunction(PyObject* pyobject, const char* name)
{
  if(!(pyobject && PyCallable_Check(pyobject)))
    throw Exception("Cannot call Python function " + string(name ? name : ""));
}


string PyErr(bool call_PyErr_Print)
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

std::string ToString(PyObject* obj)
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
Matrix<T> Numpy2DArrayToMatrix(PyObject* numpy, Matrix<T>& out, bool decref)
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
void MatrixToNumpyArray(const Matrix<T>& in, PyObject* numpy)
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

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template Matrix<double> Numpy2DArrayToMatrix<double>(PyObject*, Matrix<double>&, bool);
  template Matrix<int> Numpy2DArrayToMatrix<int>(PyObject*, Matrix<int>&, bool);

  template void MatrixToNumpyArray<double>(const Matrix<double>&, PyObject*);
  template void MatrixToNumpyArray<int>(const Matrix<int>&, PyObject*);
#endif


} // end of namespace
