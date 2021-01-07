#ifndef OPTIMIZATION_SPAGHETTI_HH_
#define OPTIMIZATION_SPAGHETTI_HH_

#include "Optimization/SIMP.hh"

namespace CoupledField {

/** Feature mapping variant flying spaghetti (SpaghettiDesign) */
class Spaghetti : public SIMP
{
public:
  Spaghetti() : SIMP()
  {
  }

  virtual ~Spaghetti()
  {
  }

};

} // end of namespace

#endif /* OPTIMIZATION_SPAGHETTI_HH_ */
