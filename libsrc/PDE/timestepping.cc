#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"

namespace CoupledField
{

TimeStepping :: TimeStepping(std::string apdename, BaseSystem * algebraicsystem)
{
#ifdef TRACE
  (*trace) << "entering TimeStepping::TimeStepping" << std::endl;
#endif

  pdename_ = apdename;
  algsys_  = algebraicsystem;

}

TimeStepping :: ~TimeStepping()
{
#ifdef TRACE
  (*trace) << "entering TimeStepping::~TimeStepping" << std::endl;
#endif

}




} // end of namespace
