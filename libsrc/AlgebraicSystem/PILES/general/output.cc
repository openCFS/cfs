#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

#include "output.hh"

namespace CoupledField
{

ostream * trace    = NULL;
ostream * debug    = NULL;
ostream * test     = NULL;
ostream * conv     = NULL;
ostream * memtrace = NULL;

  double sumdmem = 0;
  double sumimem = 0;
}
