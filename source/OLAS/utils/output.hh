// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_OUTPUT_HH
#define OLAS_OUTPUT_HH

#include <iostream>


// namespace OLAS 
// {
// #ifdef TRACE
//   extern std::ostream * trace;
// #endif

// #ifdef DEBUG
//   extern std::ostream * debug;
//   extern std::ostream * test;
// #endif

// #ifdef MEMTRACE
//   extern std::ostream * memtrace;

//   extern double sumdmem;
//   extern double sumimem;
// #endif


// } // namespace


namespace OutInfo
{
#ifdef TRACE
  extern std::ostream * trace;
#endif

#ifdef DEBUG
  extern std::ostream * debug;
  extern std::ostream * test;
#endif

#ifdef MEMTRACE
  extern std::ostream * memtrace;

  extern double sumdmem;
  extern double sumimem;
#endif
 } // namespace


namespace OLAS 
{
#ifdef TRACE
  using OutInfo::trace;
#endif

#ifdef DEBUG
  using OutInfo::debug;
  using OutInfo::test;
#endif

#ifdef MEMTRACE
  using OutInfo::memtrace;
  using OutInfo::sumdmem;
  using OutInfo::sumimem;
#endif

} // namespace




#endif // OLAS_OUTPUT_HH
