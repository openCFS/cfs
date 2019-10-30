#ifndef SHAPEOPT_HH_
#define SHAPEOPT_HH_

#include <stddef.h>
#include <string>

#include "General/defs.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/ParamMat.hh"

namespace CoupledField {

class Condition;
class Excitation;
class Function;
class Objective;
class SplineBoxDesign;

class SplineBoxOpt : public ParamMat {
public:

  SplineBoxOpt();

  virtual ~SplineBoxOpt();

  virtual void PostInit();

};

}  // end of namespace

#endif // SHAPEOPT_HH_
