#include <iostream>
#include <fstream>
#include <string>

#include "timeerror.hh"

namespace CoupledField
{

TimeErrorEstimator::TimeErrorEstimator(BasePDE* aptPDE)
{
#ifdef TRACE
  (*trace) << "entering TimeErrorEstimator::TimeErrorEstimator" << std::endl;
#endif

  if (!aptPDE) Error("Wrong type of pointer in constructor of TimeError");
 
  ptPDE_=aptPDE;
}

} // end of namespace
