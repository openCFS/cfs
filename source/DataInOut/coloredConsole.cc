//------------------------------------------------------------------------------
// coloredConsole.cc: instantiate ColoredConsole::colorise member
//------------------------------------------------------------------------------

#include "coloredConsole.hh"

namespace CoupledField
{
  bool ColoredConsole::colorise = true;

#ifndef SUPPRESS_COLORED_OUTPUT
  const bool ColoredConsole::suppressed = false;
#else
  const bool ColoredConsole::suppressed = true;
#endif
    
}
