#include <fstream>
#include <iostream>
#include <string>

#include "BaseEQN.hh"

namespace CoupledField
{

BaseEQN :: BaseEQN(Grid * aptgrid, BCs *aptBCs, std::vector<std::string>& asubdoms, Integer actlevel)
{
#ifdef TRACE
  (*trace) << "entering BaseEQN::BaseEQN" << std::endl;
#endif

  ptgrid_   = aptgrid;
  ptbcs_    = aptBCs;
  subdoms_  = asubdoms;
  actlevel_ = actlevel;

}

BaseEQN :: ~BaseEQN()
{
#ifdef TRACE
  (*trace) << "entering BaseEQN::~BaseEQN" << std::endl;
#endif

}




} // end of namespace
