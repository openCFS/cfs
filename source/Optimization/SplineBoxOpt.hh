#ifndef SPLINEBOXOPT_HH_
#define SPLINEBOXOPT_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField {

class SplineBoxOpt : public SIMP
{
public:

  SplineBoxOpt();

  virtual ~SplineBoxOpt();

  virtual void PostInit();

};

}  // end of namespace

#endif // SPLINEBOXOPT_HH_
