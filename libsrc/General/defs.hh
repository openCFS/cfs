#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

#include "Utils/debugger.hh"
#include <typeinfo>
#include <exception>
#include <new>

namespace CoupledField{

// ****************************************************************************
//   This block deals with function tracing
// ****************************************************************************

#ifdef TRACE //normal function tracing
#define ENTER_FCN(name)	\
FcnObj fcn(name);

#if TRACE>=100 //trace absolutely everything
#define ENTER_IFCN(name)	\
FcnObj fcn(name);
#else //no tracing
#define ENTER_IFCN(name)
#endif//IFCN
#else//no tracing
#define ENTER_FCN(name)	
#define ENTER_IFCN(name)	
#endif



#ifdef CHECK_TYPE_CASTS

#define TRY_CAST try {
#define CATCH_CAST }catch(std::bad_alloc e)\
{Error( WRONG_CAST_MSG, __FILE__, __LINE__ );}

#define PTRCAST(NAME,TYPE,TARGET)		\
TYPE* TARGET = dynamic_cast<TYPE*>(NAME);	\
if (TARGET==NULL) Error( WRONG_CAST_MSG, __FILE__, __LINE__ );

#else
#define TRY_CAST
#define CATCH_CAST

#define PTRCAST(NAME,TYPE,TARGET)\
TYPE* TARGET = dynamic_cast<TYPE*>(NAME);	

#endif

#define REFCAST(NAME,TYPE,TARGET)\
TYPE& TARGET = dynamic_cast<TYPE&>(NAME);	
#define CONSTREFCAST(NAME,TYPE,TARGET)\
const TYPE& TARGET = dynamic_cast<const TYPE&>(NAME);	

//! Error message signaling a dynamic miscast.
#define WRONG_CAST_MSG " Invalid cast attempt! "

} // end of namespace
#endif
