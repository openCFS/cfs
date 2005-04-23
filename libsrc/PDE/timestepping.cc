#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"
#include "olas.hh"

namespace CoupledField {

TimeStepping::TimeStepping(std::string apdename, BaseSystem * algebraicsystem, 
						   NodeEQN * ptEQN)
{

  ENTER_FCN( "TimeStepping::TimeStepping" );

  pdename_ = apdename;
  algsys_  = algebraicsystem;
  ptEQN_ = ptEQN;
}

TimeStepping::~TimeStepping() {

  ENTER_FCN( "TimeStepping::~TimeStepping" )
}


} // end of namespace
