#ifndef COEFFUNCTIONPYTHON_HH_
#define COEFFUNCTIONPYTHON_HH_

#include <boost/shared_ptr.hpp>

#include "CoefFunction.hh"

namespace CoupledField
{


/** this is python based variant of CoefFunctionExpression. Unlike e.g. the python mesh reader it does not
 * hold a own python script (module) but is based on the script (module) of the singleton python kernel.
 * This is defined by the python element in the document root.
 * There are no options to be passed from xml to the function. */
class CoefFunctionPython : public CoefFunction, public boost::enable_shared_from_this<CoefFunctionPython>
{
public:

  /** <python> Element, e.g. for mechanic force */
  CoefFunctionPython(PtrParamNode pn, unsigned int dim);

  virtual ~CoefFunctionPython();

  string GetName() const override { return "CoefFunctionPython"; }

  void GetScalar(double& scal, const LocPointMapped& lpm) override;

  void GetScalar(Complex& scal, const LocPointMapped& lpm)  override{
    EXCEPTION("complex CoefFunctionPython not implemented yet");
  }

  /** the node of interest is transported in lpm.pl.number */
  void GetVector(Vector<double>& vec, const LocPointMapped& lpm)  override;

  void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm)  override {
    EXCEPTION("complex CoefFunctionPython not implemented yet");
  }

  std::string ToString() const override { return function_; }

  /** optional value from xml. If true cfs will add an CoefFunctionCompound */
  bool DoNormalize() const override { return normalize_; }

  unsigned int GetVecSize() const override { return dim_; }

protected: // we have protected instead of private to avoid "private field is not used" errors when compiled w/o

  /** do the actual function call */
  pyObject* CallFunction(const LocPointMapped& lpm);

  unsigned int dim_ = 0;

  /** name of the optional init function, gets the optional options */
  std::string init_;

  /** name of the mandatory evaluation function, Called with a list of coordinates */
  std::string function_;

  /** xml attribute, if true, an CoefFunctionCompound will divide by number of entries */
  bool normalize_;

  /** xml attribute, if we call the eval Python function by coordinate or lpm.lp.number */
  bool by_coord_;

#ifdef USE_EMBEDDED_PYTHON
  /** function object to be repeatedly called */
  pyObject* eval_ = NULL;
#endif
};


} // end of namespace

#endif /* COEFFUNCTIONPYTHON_HH_ */
