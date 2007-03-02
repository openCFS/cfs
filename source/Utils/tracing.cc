#include "tracing.hh"

namespace OutInfo{

  // Definition and initialization of static variables of
  // class FcnTraceHandler has to be outside of class scope
#ifdef TRACE
  FcnTraceListElem* FcnTraceHandler::foo_ = NULL;
  TraceInt FcnTraceHandler::fcnTraceDepth_ = 0;
  TraceInt FcnTraceHandler::fcnTraceDepthLimit_ = 0;
#endif
}

