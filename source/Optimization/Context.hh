#ifndef OPTIMIZATION_CONTEXT_HH_
#define OPTIMIZATION_CONTEXT_HH_

#include <cstddef>

namespace CoupledField
{

class Excitation;

/** The context describes where we are within the current optimization process.
 * A basic information is the currently applied Excitation or for a multi state
 * problem the current state of interest. */
class Context
{
  public:

  Context() {
    excitation = NULL;
  }

  Excitation* excitation;
};

} // end of namespace


#endif /* OPTIMIZATION_CONTEXT_HH_ */
