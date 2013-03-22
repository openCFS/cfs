//------------------------------------------------------------------------------
// coloredConsole.cc: instantiate ColoredConsole::colorise member
//------------------------------------------------------------------------------

#include "ColoredConsole.hh"

namespace CoupledField
{
  bool ColoredConsole::colorise = true;

#ifndef SUPPRESS_COLORED_OUTPUT
  bool ColoredConsole::suppressed = false;
#else
  bool ColoredConsole::suppressed = true;
#endif

}
