#include "tracing.hh"

namespace OutInfo{

  // Definition and Initialization of static variables of class FcnTraceHandler
  // has to be outside of class scope
#ifdef TRACE
  FcnTraceListElem* FcnTraceHandler::foo_ = NULL;
  unsigned int FcnTraceHandler::fcnTraceDepth_ = 0;
  unsigned int FcnTraceHandler::fcnTraceDepthLimit_ = 0;
#endif
}

