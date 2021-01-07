#ifndef OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_
#define OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_

#include <def_use_embedded_python.hh>

#ifdef USE_EMBEDDED_PYTHON
  #define PY_SSIZE_T_CLEAN // https://docs.python.org/3/c-api/intro.html
  #include <Python.h>
#endif

#include <utility>
#include <tuple>
#include <boost/filesystem.hpp>
#include "Optimization/Optimizer/BaseOptimizer.hh"

namespace CoupledField
{
/** Optimizer for shape- and topology optimization! */
class PythonOptimizer : public BaseOptimizer
{
public:
  PythonOptimizer(Optimization* optimization, PtrParamNode pn);
  virtual ~PythonOptimizer();

  /** expects two 1Dim Arrays for design bounds and two 1Dim arrays for constraint bounds. Not normalized, equal is both the same */
  void GetBounds(PyObject* args);

  void GetInitialDesign(PyObject* args);

  /** Evaluate objective and expects 1Dim design as numpy array */
  double EvalObjective(PyObject* args);

  /** automatically makes this a new iteration, could be changed by function to by called by python */
  void EvalGradObjective(PyObject* args);

  /** Evaluate constraints and expects 1Dim design as numpy array and 1Dim array of size m for values */
  void EvalConstraints(PyObject* args);

  void EvalGradConstraints(PyObject* args);

  /** this should be in PythonTools but somehow this gives a segfault which indicates that PyArg_ParseTuple is not defined?! */
  static StdVector<PyObject*> ParseArrays(PyObject* args, int expect, StdVector<Vector<double> >& data, bool decref);

protected:

  /** Implements virtual function from BaseOptimizer */
  void SolveProblem();
private:


  /** the options from the xml file */
  StdVector<std::pair<std::string, std::string> > options;

  /** given filename */
  boost::filesystem::path givenname;

  /** the module givenname opened as embedded python environment */
  PyObject* module = NULL;

  /** number of design variables */
  unsigned int n = 0;

  /** number of constraints */
  int m = -1;

  /** our python info instead of optimizer info_ */
  PtrParamNode pyinf_ = NULL;

};


} // end of namespace

#endif /* OPTIMIZATION_OPTIMIZER_PYTHONOPTIMIZER_HH_ */
