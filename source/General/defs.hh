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

//#define USE_ADOLC 1 // why does this not autmatically work with include/def_use_adolc.hh.in and cmake USE_ADOLC option?
#include "def_use_adolc.hh" // this is a really bad idea. where do i set this to read automatically?
#ifdef USE_ADOLC
#include "adolc/adtl.h"
#endif


namespace CoupledField{

// if AD is to be used, then use adouble data type
#ifdef USE_ADOLC 
typedef boost::int32_t Integer;
typedef boost::uint32_t UInt;
typedef boost::int16_t ShortInt;
typedef adtl::adouble adouble;
typedef adouble Float;
typedef adouble Double;
typedef adouble Adouble; //for good measure
typedef std::complex<Double> Complex;
#else 
// if not, then use standard
  //! redeclaration of types
typedef boost::int32_t Integer;
typedef boost::uint32_t UInt;
typedef boost::int16_t ShortInt;
typedef float Float;
typedef double Double;
typedef std::complex<Double> Complex;
#endif



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
}catch (std::bad_alloc exception){\
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
