#ifndef PYTHONKERNEL_HH_
#define PYTHONKERNEL_HH_

#include <def_use_embedded_python.hh>
#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#endif

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Utils/Container.hh"
#include "Optimization/Design/FeaturedDesign.hh"
#include <string>
#include <utility>


namespace CoupledField
{
class SimInputPython;
class PythonOptimizer;
class FeaturedDesign;
class PythonKernel;

/** The PythonKernel always exists as a global instance to be accessed via this pointer. Defined in tools.cc */
extern PythonKernel* python;

  /** This class exists always in exactly one instance.
   * It handles the embedded python interpreter (USE_EMBEDDED_PYTHON) which the system python lib and its
   * installed modules (e.g. numpy is necessary). The python interpreter gets cfs_modules to be
   * imported via "import cfs" from python. The implementations are static C functions collected in
   * PythonKerneFunctions.cc. This needs to happen very early, as e.g. the python mesh reader needs this functions.
   * Therefore during runtime python reader, python optimizer, ... register themselves here to make
   * the functions viable. The PythonKernel optionally can hold a python module (script) for auxiliary
   * usage (e.g. hooks). This is the kernel script.
   *
   * There are three principle approaches to use embedded python:
   * 1) have a full featured script (like python mesh reader). The use LoadPythonModule() and do the
   *    rest yourself.
   * 2) use the script defined in the root level python element (kernel script via GetKernel()) as in CoefFunctionPython.
   *    The called python function might call cfs functions. This is a lightweight implementation.
   * 3) use hooks to call a function from the kernel script at defined points. E.g. to trigger
   *    something like mapped equations output, ...
   *
   * Note that independent from the interpreter is importing python scripts (e.g. python mesh reader).
   * It needs just the interpreter being initialized, but requires no handle. This scripts are also called
   * modules, don't confuse them with the modules for the interpreter!
   *
   * Anytime USE_EMBEDDED_PYTHON is enabled, the interpreter is called on startup (when not called with --noPython)
   * as we cannot know so early if we want it later.
   *
   * See the corresponding test case (search for <python> and embeddedpython.cc. To google, search for "python c api" */
class PythonKernel
{
public:
  typedef enum { ASSEMBLE_RHS = 0, OPT_EVAL_FUNC, OPT_EVAL_GRAD, OPT_POST_INIT, OPT_POST_ITER, POST_DOMAIN_INIT, POST_GRID, POST_SOLVE_PROBLEM} Hook;

  static Enum<Hook> hook;

  /** make sure to instantiate this only once
   * @param pn the <python> element before <fileFormats>. Can be NULL
   * @param info only written, when we do something */
  PythonKernel(PtrParamNode pn, PtrParamNode info)
  {
#ifdef USE_EMBEDDED_PYTHON
    Init(pn,info);
#else
    if(pn != nullptr)
      throw Exception("'python' element given but compiled w/o USE_EMBEDDED_PYTHON: " + pn->GetLocation());
#endif
  }

  /** called when cfs dies */
  ~PythonKernel()
  {
#ifdef USE_EMBEDDED_PYTHON
    Release();
#endif
  }

  /** optional python module (script). Can be NULL */
  PyObject* GetKernel() { return kernel_; }

  /** to by called by the pathon mesher to forward the static functions from cfs_modules_, called from python to the class */
  void Register(SimInputPython* mesher, bool remove = false) { this->mesher_ = remove ? NULL : mesher; }

  /** register PythonOptimizer. Note that Optimization is not registered but taken from domain->GetOptimization() */
  void Register(PythonOptimizer* opt, bool remove = false) { this->pyopt_ = remove ? NULL : opt; }

  /** return value of LoadPythonModule. How easy this is with Python :( */
  struct LoadStatus
  {
    /** the python module (script) to call functions from. You need to free it after usage! */
    PyObject* module = NULL;
    /** the actually used path. But module loading might mess this up. Save is to start with empty PYTHONPATH */
    std::string full_file;
    /** the actually used system path. The path of the current script is set first */
    StdVector<std::string> sys_path;
  };

  /** Loads a python module (script) into the running interpreter. Conditionally initializes the interpreter
   * The interpreter got it's cfs_modules (functions to be called from python) when initializing the interpreter.
   * @param file the python module to be loaded. Can include a full path if attribute path is empty
   * @param path optional path (leave empty if not given). "cfs:share:python" is a special key which uses ProgramOption->GetSchemaPath() + "../python" */
  LoadStatus LoadPythonModule(const std::string& file, std::string path = "");

  /** Convenience function which extracts file and possibly path */
  LoadStatus LoadPythonModule(PtrParamNode pn);

  /** Call this in the appropriate locations in cfs. If the hook is registerd, the corresponding
   * python function is called. */
#ifdef USE_EMBEDDED_PYTHON
  static void CallHook(Hook hook);
#else
  static void CallHook(Hook hook) {}
#endif

    /** helper which from all elements "key" and "value" */
  static StdVector<std::pair<std::string, std::string> > ParseOptions(ParamNodeList lst);

  /** helper which creates a python dictionary for a string resource.
   * @see ParseOptions() */
  static PyObject* CreatePythonDict(const StdVector<std::pair<std::string, std::string> > options);

  /** helper which creates a python list of entities which need to be specialized in Create().
   * Any TYPE needs a specialization of Create<>(). */
  template<class TYPE>
  static PyObject* CreatePythonList(const StdVector<TYPE>& list) {
    return CreatePythonList(Container<TYPE>(list));
  }

  template<class TYPE>
  static PyObject* CreatePythonList(const Vector<TYPE>& list) {
    return CreatePythonList(Container<TYPE>(list));
  }

  /** common implementation for Vector and StdVector */
  template<class TYPE>
  static PyObject* CreatePythonList(const Container<TYPE>& list);

  /** convert a single value to PyObject. Works with specialization. */
  template<class TYPE>
  static PyObject* Create(const TYPE& val);

  template<class TYPE>
  static void ConvertPythonList(Vector<TYPE>& vec, PyObject* list) {
    Container<TYPE> cont(vec);
    ConvertPythonList(cont, list);
  }

  template<class TYPE>
  static void ConvertPythonList(StdVector<TYPE>& vec, PyObject* list) {
    Container<TYPE> cont(vec);
    ConvertPythonList(cont, list);
  }

  template<class TYPE>
  static void ConvertPythonList(Container<TYPE>& vec, PyObject* list);


  /** convenience function */
  template<class TYPE>
  static StdVector<TYPE> ConvertPythonList(PyObject* list) {
    StdVector<TYPE> vec;
    ConvertPythonList<TYPE>(vec, list);
    return vec;
  }

  /** counterpart of Create(). Has specialization for double, ... */
  template<class TYPE>
  static TYPE Convert(PyObject* item);

  /** returns the string representation of a PyObject()
   * If obj is NULL an empty string is returned */
  static std::string ToString(PyObject* obj);

  /** convenience function which checks the return value and if it fails calls PyErr_Print() and throws an exception
   * @param pyobject e.g. what you get from PyObject_CallObject
   * @param name if given, used for thrown error message */
  static void CheckPythonReturn(PyObject* pyobject, const char* name = NULL); // { CheckPythonReturn(pyobject == NULL ? 0 : 1, name); }
  static void CheckPythonReturn(int ret, const char* name = NULL);

  /** convenience function which checks if the python function is callable
   * @param pyobject e.g. what you get from PyObject_GetAttrString
   * @param name optional function name for error message */
  static void CheckPythonFunction(PyObject* pyobject, const char* name = NULL);

  /** if an exception had been occurred and PyErr_Print() would print an stacktrace
   * @param call PyErr_Print()
   * @return error message as string (type + value) */
  static std::string PyErr(bool call_PyErr_Print = true);

  /** checks if mesher has been called. It registers a python exception and requires to return with NULL */
  static bool CheckMesher();

  /** check if the PythonOptimizer is registered */
  static bool CheckPyOpt();

  /** gets a the FeaturedDesign design space or throws an error with readable message */
  static FeaturedDesign* GetFeaturedDesign();

  /** check if we have Optimization at all.
   * @see CheckPyOpt() */
  static bool CheckOpt();

  /** make this public such that the c-functions can directly call */
  static PythonOptimizer* GetPyOpt() { return pyopt_; }

  /** DOES NOT WORK! SHALL REPLACE PythonOptimizer::ParseArrays() BUT SEGAULTS ?!
     * Helper which processes a PyTupleObject which needs to consist only of 1dim numpy arrays
     * @param decref if false make sure to decref the objects via the return array
     * @return you must not uses the PyObjects when decref is true */
  static StdVector<PyObject*> ParseNumpyArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref);

  /** Write 2D numpy array data to a Matrix.
   * @param numpy PyArrayObject*
   * @param out is properly set
   * @param decref shall numpy be decrefed? */
  template<class T>
  static Matrix<T> Numpy2DArrayToMatrix(PyObject* numpy, Matrix<T>& out, bool decref);

  template<class T>
  static Matrix<T> Numpy2DArrayToMatrix(PyObject* numpy, bool decref)
  {
    Matrix<T> out;
    return Numpy2DArrayToMatrix(numpy, out, decref);
  }

  /** Write Matrix data to numpy array which needs already to be of proper size */
  template<class T>
  static void MatrixToNumpyArray(const Matrix<T>& in, PyObject* numpy);

private:

  /** actual constructor */
  void Init(PtrParamNode pn, PtrParamNode info);

  /** actual destructor */
  void Release();

  /** initialized the python interpreter which runs in a single instance for all modules (scripts) and
   * loads the c functions for the interpreter before initialization. Register() can be called later
   * but before the corresponding LoadPythonModule() */
  void InitInterpreter();

  struct HookParam
  {
    std::string function; // when empty, not defined */s
    StdVector<std::pair<std::string, std::string> > options;
  };

  std::map<Hook, HookParam> hooks_;

  /** not initialized when not compiled with USE_EMBEDDED_PYTHON or when called with --noPython */
  bool initialized_ = false;

  /** This is the kernel module (python script) if loaded. Not that mesh reader, python optimizer, ... have their
   * own modules (files) */
  PyObject* kernel_ = NULL;

  PtrParamNode info_;

  /** to forward the python mesher calls from python (set_nodes, ...)
   * @see Register() */
  static SimInputPython* mesher_;
  static PythonOptimizer* pyopt_;

// we can have a dummy PyObject* from Vector but not these structures
#ifdef USE_EMBEDDED_PYTHON
  /** all callable functions to be packed to the cfs module to be imported by python scripts.
   * Set in PythonKernelFunctions.cc */
  static PyMethodDef cfs_methods[];

  /** we pack all callable functions (mesher, optimizer, ...) in a single module "cfs" to be imported by python scripts */
  static PyModuleDef cfs_modules;
#endif

  /** give the python interpreter the static functions to be called from python via import cfs*/
  static PyObject* SetModulFunctions(void);


  /** make the static functions class members to allow access to private members (via friend)
   * The other functions are non-class in PythonKernel.cc  */
  static PyObject* mesher_set_nodes(PyObject* self, PyObject* args);
  static PyObject* mesher_set_regions(PyObject* self, PyObject* args);
  static PyObject* mesher_add_elements(PyObject* self, PyObject* args);
  static PyObject* mesher_add_named_nodes(PyObject* self, PyObject* args);
  static PyObject* mesher_add_named_elements(PyObject* self, PyObject* args);

}; // end of class

} // end of namespace

#endif /* PYTHONKERNEL_HH_ */
