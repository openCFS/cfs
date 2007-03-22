// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

