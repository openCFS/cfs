#include "tracing.hh"

namespace CoupledField{

  // Definition and Initialization of static variables of class Debugger
  // has to be outside class scope
#ifdef TRACE
  FcnTraceListElem* FcnTraceHandler::foo_ = NULL;
  unsigned int FcnTraceHandler::fcnTraceDepth_ = 0;
#endif
}

