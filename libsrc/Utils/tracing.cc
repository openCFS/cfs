#include "debugger.hh"

namespace CoupledField{

  // Definition and Initialization of static variables of class Debugger
  // has to be outside class scope
#ifdef TRACE
  FcnTrace* Debugger::foo_ = NULL;
  unsigned int Debugger::fcnTraceDepth_ = 0;
#endif
}

