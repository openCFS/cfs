#include "debugger.hh"

namespace CoupledField{

  // Definition and Initialization of static variables of class Debugger
  // has to be outside class scope
#ifdef TRACE
  fcn_trace* Debugger::foo_ = NULL;
  int Debugger::fcn_trace_depth_ = 0;
#endif
}

