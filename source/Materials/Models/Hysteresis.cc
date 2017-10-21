#include <iostream>
#include <fstream>


#include "Hysteresis.hh"

#include "Utils/tools.hh"

namespace CoupledField
{ 
  Hysteresis::Hysteresis(Integer numElem)
  {

    numElements_ = numElem;
  }

  Hysteresis::~Hysteresis()
  {
  }
}
