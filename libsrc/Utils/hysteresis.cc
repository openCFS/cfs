#include <iostream>
#include <fstream>

#include "Utils/tools.hh"
#include "hysteresis.hh"

namespace CoupledField
{ 

Hysteresis :: Hysteresis(Integer numElem)
{
  ENTER_FCN("Hysteresis::Hysteresis" );

  numElements_ = numElem;

}

Hysteresis :: ~Hysteresis()
{
}


}
