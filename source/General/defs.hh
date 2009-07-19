// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

#include <iostream>
#include <string>
#include <exception>
#include <typeinfo>
#include <complex>
#include <limits>

#include <boost/cstdint.hpp>

// Include build type options header containing the #defines
#include <def_build_type_options.hh>
#include <def_xmlschema.hh>
// Include headers which define if CFS++ should support scripting
#include <def_use_scripting.hh>

#include "OLAS/utils/profiler.hh"

#ifdef CHECK_MEM_ALLOC
#include <new>
#endif

#ifdef PROFILING
#include "Utils/profiler.hh"
#endif

namespace CoupledField{

//! redeclaration of types
typedef boost::int32_t Integer;
typedef boost::uint32_t UInt;
typedef boost::int16_t ShortInt;
typedef float Float;
typedef double Double;
typedef std::complex<Double> Complex;
  // *************************************************************************
  //   This block deals with profiling information
  // *************************************************************************

#ifdef PROFILING
#define SETPROFILE(name)\
profiler->Trace(name)
#else
#define SETPROFILE(name)
#endif

//! Can be used for outputting the value of a macro (expansion gives "a")

//! This macro can be used to obtain the value of the macro at compile time
//! in quotes. This is useful e.g. for passing the value to a stream. It
//! requires a second parameterised macro QUOTEMACRO for its operation.
//! The following piece of code
//! \code
//! #define LAENGE 63
//! std::cout << MACRO2STRING(LAENGE)
//! \endcode
//! looks like this after preprocessing
//! \code
//! std::cout << QUOTEMACRO(63)
//! \endcode
//! which is further expanded to
//! \code
//! std::cout << "63"
//! \endcode
#define MACRO2STRING(a) QUOTEMACRO(a)

//! Auxiliary macro needed by MACRO2STRING
#define QUOTEMACRO(a) #a


//! Value to be used in penalty approach to Dirichlet boundary conditions
#define PENALTY 1e12


// Some compilers, e.g. SGI's CC, seem to have problems with undefined
// parameterised macros. They do not evaluate them to 0 as for
// non-parameterised undefined macros. So to be on the safe side
#ifndef __GNUC__
#define __GNUC_PREREQ(major,minor) 0
#endif

// ****************************************************************************
//   This block deals definition of the CFS++/OLAS "standard" types
// ****************************************************************************

// ****************************************************************************
//   This block deals with function tracing
// ****************************************************************************
//! \def PROFILE(name,ops)
//! This macro enables profiling for the current function
//! name should be a string identifying the function, ops
//! is an estimate on the number of floating-point operations
//! performed in the function.
//! the macro is ignored if PROFILING is not defined

#ifdef PROFILING
#define PROFILE(name,ops)\
FcnProf fcnprof(name,ops);
#else 
#define PROFILE(name,ops)
#endif

// ****************************************************************************
//   This block deals with dynamic memory allocation
// ****************************************************************************

// NOTE: We put the doxygen documentation here in front of the actual defines.
//       The reason is that doxygen also does preprocessing and the
//       documentation might not appear otherwise.

//! \def New
//! \brief Non-throwing new

//! \def New
//! By default the new operator throws an exception, when the memory allocation
//! fails. This macro is a shortcut for std::new(std::nothrow), which tells
//! new to use the old C-style behaviour and return a NULL pointer instead.

//! \def ASSERTMEM(name,size)
//! Assert memory after New and write to memtrace

//! \def NEWARRAY(name,type,size)
//! convenient array allocation (saves the ASSERTMEM, doesn't work for
//! more than 1 template arg)
//! \note ISO C++ forbids initialization in array new, so no need for another
//! fourth parameter here.

//! \def DELETEARRAY(name)
//! Convenience macro for deleting an array.
//! Asserts the pointer, increments it and calls delete.

#ifdef MEMTRACE

#ifdef CHECK_MEM_ALLOC

  // #define New new(std::nothrow)

#define ASSERTMEM(name,size)\
if (!name) EXCEPTION("memory allocation for array failed");\
(*memtrace) << "allocated " << (size)*1e-6 << " MB in " << __FILE__ \
            << ", line " << __LINE__ << std::endl;

#define NEWARRAY(name,type,size) { \
try{\
name = new type [(size)];\
}catch (std::bad_alloc exception){\
EXCEPTION("memory allocation for array failed"); }\
(*memtrace) << "allocated " << (size)*(sizeof(type))*1.0e-6 \
<< " MB in " << __FILE__ << ", line " << __LINE__ << std::endl; }

#else //no check_mem_alloc but memtrace

// #define New new
#define ASSERTMEM(name,size)  \
(*memtrace) << "allocated " << size*1e-6 << " MB in " << __FILE__ << ", line " << __LINE__ << std::endl;
#define NEWARRAY(name,type,size) { \
name = new type[(size)]; \
(*memtrace) << "allocated " << (size)*(sizeof(type))*1e-6 \
<< " MB in " << __FILE__ << ", line " << __LINE__ << std::endl; }
#endif//check_mem_alloc

#else//no memtrace

#ifdef CHECK_MEM_ALLOC

// #define New new(std::nothrow)
#define ASSERTMEM(name,size)  \
if (!name) EXCEPTION("memory allocation failed"); \

#define NEWARRAY(name,type,size) {\
try{\
name = new type[(size)];\
}catch (std::bad_alloc exception){\
EXCEPTION("Memory allocation for array failed\n size = " << (size));\
}\
}
#else //no check_mem_alloc, no memtrace
// #define New new
#define ASSERTMEM(name,size)
#define NEWARRAY(name,type,size) { name = new type[(size)]; }

#endif//check_mem_alloc

#endif //memtrace

// Using delete with a null pointer is assured to be safe, so the assertion is
// superfluous. However, if you increment NULL and then delete it, that is not
// a good idea.
//#define DELETEARRAY(name) if (name) { delete[] (name+1); name = NULL; }
#define DELETEARRAY(name) if (name) { delete[] (name); name = NULL; }


// ****************************************************************************
//   This block deals with dynamic type casting
// ****************************************************************************

#ifdef CHECK_TYPE_CASTS

#define TRY_CAST try {
#define CATCH_CAST } catch( std::bad_cast e )\
{EXCEPTION( WRONG_CAST_MSG );}

#define PTRCAST(NAME,TYPE,TARGET)   \
TYPE* TARGET = dynamic_cast<TYPE*>(NAME); \
if (TARGET==NULL) EXCEPTION( WRONG_CAST_MSG );

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
#define WRONG_CAST_MSG "Invalid cast attempt!"


}
#endif
