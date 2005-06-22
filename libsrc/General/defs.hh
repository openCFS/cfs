#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

#ifdef TRACE
#include "Utils/tracing.hh"
#endif

#ifdef PROFILING
#include "Utils/profiler.hh"
#endif

#include <typeinfo>

namespace CoupledField{

  // ****************************************************************************
  //   This block deals with function tracing
  // ****************************************************************************

#ifdef TRACE //normal function tracing
#define ENTER_FCN(name) \
OutInfo::FcnTraceObjLocal fcn(name);

#if TRACE>=100 //trace absolutely everything
#define ENTER_IFCN(name)        \
OutInfo::FcnTraceObjLocal fcn(name);
#else //no tracing
#define ENTER_IFCN(name)
#endif//IFCN
#else//no tracing
#define ENTER_FCN(name) 
#define ENTER_IFCN(name)        
#endif

  // ****************************************************************************
  //   This block deals with profiling information
  // ****************************************************************************
#ifdef PROFILING
#define SETPROFILE(name)\
profiler->Trace(name)
#else
#define SETPROFILE(name)
#endif

} // end of namespace
#endif
