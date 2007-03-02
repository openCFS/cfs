#include <iostream>
#include <fstream>

#include "utils/output.hh"

// namespace OLAS 
namespace OutInfo
{
  std::ostream *cla = NULL;
  
#ifdef TRACE
  std::ostream * trace    = NULL;
#endif

#ifdef DEBUG
  std::ostream * debug    = NULL;
  std::ostream * test     = NULL;
#endif

#ifdef MEMTRACE
  std::ostream * memtrace = NULL;

  double sumdmem = 0;
  double sumimem = 0;
#endif
}
