#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"

namespace CoupledField
{

TimeStepping :: TimeStepping(std::string apdename, BaseSystem * algebraicsystem)
{
  ENTER_FCN( "TimeStepping::TimeStepping" );

  pdename_ = apdename;
  algsys_  = algebraicsystem;

}

TimeStepping :: ~TimeStepping()
{
  ENTER_FCN( "TimeStepping::~TimeStepping" )

}




} // end of namespace
