// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

#include <exception>
#include <typeinfo>
#include <complex>
#include <limits>
#include <memory>
#include <cstdint>

// Include build type options header containing the #defines
#include <def_build_type_options.hh>
#include <def_xmlschema.hh>

namespace CoupledField{

//! redeclaration of types
typedef int32_t  Integer;
typedef uint32_t UInt;
typedef int16_t  ShortInt;
typedef float Float;
typedef double Double;
typedef std::complex<Double> Complex;

class ParamNode;  
typedef std::shared_ptr<ParamNode> PtrParamNode;
typedef std::weak_ptr<ParamNode> WeakPtrParamNode;
template<typename T> class StdVector; 
typedef StdVector<std::shared_ptr<ParamNode> > ParamNodeList;

// Type definition for shared_ptr<CoefFunction>
class CoefFunction;
typedef std::shared_ptr<CoefFunction> PtrCoefFct;
typedef std::weak_ptr<CoefFunction> WeakPtrCoefFct;

using std::shared_ptr;
using std::make_shared;
using std::weak_ptr;
using std::enable_shared_from_this;
using std::dynamic_pointer_cast;
using std::abs; // note that without code might be wrong! (better do explicitly std::abs)

#define MACRO2STRING(a) QUOTEMACRO(a)
//! Auxiliary macro needed by MACRO2STRING
#define QUOTEMACRO(a) #a

// Some compilers, e.g. SGI's CC, seem to have problems with undefined
// parameterised macros. They do not evaluate them to 0 as for
// non-parameterised undefined macros. So to be on the safe side
#ifndef __GNUC__
#define __GNUC_PREREQ(major,minor) 0
#endif

// ****************************************************************************
//   This block deals with dynamic memory allocation
// ****************************************************************************

// NOTE: We put the doxygen documentation here in front of the actual defines.
//       The reason is that doxygen also does preprocessing and the
//       documentation might not appear otherwise.

//! \def ASSERTMEM(name,size)
//! Assert memory after New and write to memtrace

//! \def NEWARRAY(name,type,size)
//! convenient array allocation (saves the ASSERTMEM, doesn't work for
//! more than 1 template arg)
//! \note ISO C++ forbids initialization in array new, so no need for another
//! fourth parameter here.
#ifdef CHECK_MEM_ALLOC

#define ASSERTMEM(name,size)  \
if (!name) EXCEPTION("memory allocation failed"); \

#define NEWARRAY(name,type,size) {\
try{\
name = new type[(size)];\
}catch (std::bad_alloc& exception){\
EXCEPTION("Memory allocation for array failed\n size = " << (size));\
}\
}
#else //no check_mem_alloc, no memtrace
#define ASSERTMEM(name,size)
#define NEWARRAY(name,type,size) { name = new type[(size)]; }

#endif//check_mem_alloc

//! Error message signaling a dynamic miscast.
#define WRONG_CAST_MSG "Invalid cast attempt!"

}
#endif
