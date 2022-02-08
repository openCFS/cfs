#ifndef OPTIMIZATION_SPAGHETTI_PARAM_MAT_HH_
#define OPTIMIZATION_SPAGHETTI_PARAM_MAT_HH_

#include "Optimization/ParamMat.hh"

namespace CoupledField {

/** Feature mapping variant flying spaghetti (SpaghettiDesign).
 * We derive from ParamMat to allow arbitrary angles. The classical case is Spaghetti.
 * We need it for the proper ParamMat::SetElementK() */
class SpaghettiParamMat : public ParamMat
{
public:
  SpaghettiParamMat() : ParamMat()
  {
  }

  virtual ~SpaghettiParamMat()
  {
  }

};

} // end of namespace

#endif /* OPTIMIZATION_SPAGHETTI_PARAM_MAT_HH_ */
