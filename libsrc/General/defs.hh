#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

#include "Utils/debugger.hh"
#include <typeinfo>

namespace CoupledField{

// ****************************************************************************
//   This block deals with function tracing
// ****************************************************************************

#ifdef TRACE //normal function tracing
#define ENTER_FCN(name)	\
FcnTraceObjLocal fcn(name);

#if TRACE>=100 //trace absolutely everything
#define ENTER_IFCN(name)	\
FcnTraceObjLocal fcn(name);
#else //no tracing
#define ENTER_IFCN(name)
#endif//IFCN
#else//no tracing
#define ENTER_FCN(name)	
#define ENTER_IFCN(name)	
#endif


} // end of namespace
#endif
