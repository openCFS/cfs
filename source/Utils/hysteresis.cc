// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <fstream>

#include "Utils/tools.hh"
#include "hysteresis.hh"

namespace CoupledField
{ 

  Hysteresis :: Hysteresis(Integer numElem)
  {

    numElements_ = numElem;
  }

  Hysteresis :: ~Hysteresis()
  {
  }


}
