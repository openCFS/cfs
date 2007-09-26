// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DEFS_2004
#define FILE_DEFS_2004

// Include build type options header containing the #defines
#include <def_build_type_options.hh>
#include <def_xmlschema.hh>

// Include headers which define if CFS++ should support scripting
#include <def_use_scripting.hh>

#ifdef PROFILING
#include "Utils/profiler.hh"
#endif

#include <typeinfo>

namespace CoupledField{

  // *************************************************************************
  //   This block deals with profiling information
  // *************************************************************************

#ifdef PROFILING
#define SETPROFILE(name)\
profiler->Trace(name)
#else
#define SETPROFILE(name)
#endif

}
#endif
