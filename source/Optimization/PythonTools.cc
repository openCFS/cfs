#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/core/include/numpy/arrayobject.h>

#include "Optimization/PythonTools.hh"
#include "DataInOut/ProgramOptions.hh"
#include "General/Exception.hh"


#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace CoupledField
{
using std::string;
using std::to_string;
using std::pair;


PyObject* InitializePythonModule(const string& file, const string& opt_path, PyObject* (*init_cfs)(), string* file_out, string* version_out)
{
  boost::filesystem::path givenname(file); // default
  if(opt_path != "")
  {
    if(opt_path == "cfs:share:python")
      givenname = progOpts->GetSchemaPath().parent_path().append("python").append(file);
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
  if(!version) {
    PyErr_Print();
    throw Exception("embedded Python not working");
  }
  const char* c_str = PyUnicode_AsUTF8(version);
  assert(c_str);
  string str(c_str);
  boost::replace_all(str, "\n"," ");
  if(version_out)
    *version_out = str;

  Py_XDECREF(version);

  // add to PYTHONPATH to make it more realistic to load the module :)
  PyObject* sysPath = PySys_GetObject((char*) "path"); // must not decref after append
  PyList_Append(sysPath, PyUnicode_FromString(absolute.string().c_str()));


  // now load the module
  boost::filesystem::path name = boost::filesystem::change_extension(givenname.filename(), "");
  PyObject* module = PyImport_ImportModule(name.string().c_str());

  if(!module)
  {
    PyErr_Print();
    throw Exception("Cannot load '" + givenname.string() + "' as embbeded python module.");
  }

  return module;
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

void CheckPythonReturn(PyObject* pyobject)
{
  if(!pyobject)
  {
    // prints something like
    // ---
    // Traceback (most recent call last):
    //  File "/Users/fwein/code/cfs/debug/../source/unittests/embeddedpython.py", line 41, in get_vec
    //    print('embeddedpython.py get_vec ->', v,n)
    // NameError: name 'v' is not defined
    // ---
    // no other way to get the error message is known
    PyErr_Print();
    // allow for stacktrace
    throw Exception("Python error within a Python function from embedded Python");
  }
}

void CheckPythonFunction(PyObject* pyobject, const char* name)
{
  if(!(pyobject && PyCallable_Check(pyobject)))
    throw Exception("Cannot call Python function '" + string(name ? name : "not-given") + "'");
}

} // end of namespace
