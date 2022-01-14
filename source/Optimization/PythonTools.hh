#ifndef PYTHONTOOLS_HH_
#define PYTHONTOOLS_HH_

#include <def_use_embedded_python.hh>
#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include <string>
#include <utility>

/** Helper to use embedded Python as done for PythonOptimizer and SpaghettiDesign */
namespace CoupledField
{
  /** Initialize the embedded Python interpreter via the modules form cfs_init to be available in
   * Python via "import cfs".
   * @param file the python module to be loaded. Can include a full path if attribute path is empty
   * @param path optional path (leave empty if not given). "cfs:share:python" is a special key which uses ProgramOption->GetSchemaPath() + "../python"
   * @param init_cfs static function like PythonOptimizer::PyInit_cfs to contain the C-Functions to be exported. Can be NULL
   * @param file_out the actually used file
   * @param version_out gets the the python version set
   * @return module which needs to be destroyed as in PythonOptimizer::~PythonOptimizer() */
  PyObject* InitializePythonModule(const std::string& file, const std::string& path, PyObject* (*init_cfs)(), std::string* file_out = NULL, std::string* version_out = NULL, StdVector<std::string>* syspath = NULL);

  /** DOES NOT WORK! SHALL REPLACE PythonOptimizer::ParseArrays() BUT SEGAULTS ?!
   * Helper which processes a PyTupleObject which needs to consist only of 1dim numpy arrays
   * @param decref if false make sure to decref the objects via the return array
   * @return you must not uses the PyObjects when decref is true */
  StdVector<PyObject*> ParseNumpyArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref);

  /** helper which from all elements "key" and "value" */
  StdVector<std::pair<std::string, std::string> > ParseOptions(ParamNodeList lst);

  /** helper which creates a python dictionary for a string resource.
   * @see ParseOptions() */
  PyObject* CreatePythonDict(const StdVector<std::pair<std::string, std::string> > options);

  void CheckPythonReturn(int ret, const char* name = NULL);
  /** convenience function which checks the return value and if it fails calls PyErr_Print() and throws an exception
   * @param pyobject e.g. what you get from PyObject_CallObject */
  void CheckPythonReturn(PyObject* pyobject, const char* name = NULL); // { CheckPythonReturn(pyobject == NULL ? 0 : 1, name); }

  /** convenience function which checks if the python function is callable
   * @param pyobject e.g. what you get from PyObject_GetAttrString
   * @param name optional function name for error message */
  void CheckPythonFunction(PyObject* pyobject, const char* name = NULL);

  /** if an exception had been occurred and PyErr_Print() would print an stacktrace
   * @param call PyErr_Print()
   * @return error message as string (type + value) */
  std::string PyErr(bool call_PyErr_Print = true);

  /** returns the string representation of a PyObject()
   * If obj is NULL an empty string is returned */
  std::string ToString(PyObject* obj);

  /** convert a python list of strings */
  StdVector<std::string> Convert(PyObject* strlist);


  /** Write 2D numpy array data to a Matrix.
   * @param numpy PyArrayObject*
   * @param out is properly set
   * @param decref shall numpy be decrefed? */
  template<class T>
  Matrix<T> Numpy2DArrayToMatrix(PyObject* numpy, Matrix<T>& out, bool decref);

  template<class T>
  Matrix<T> Numpy2DArrayToMatrix(PyObject* numpy, bool decref)
  {
    Matrix<T> out;
    return Numpy2DArrayToMatrix(numpy, out, decref);
  }

  /** Write Matrix data to numpy array which needs already to be of proper size */
  template<class T>
  void MatrixToNumpyArray(const Matrix<T>& in, PyObject* numpy);

} // end of namespace




#endif /* PYTHONTOOLS_HH_ */
