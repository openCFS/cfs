#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"
#include "olas.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(BaseSystem * algebraicsystem, UInt rhsSize )
  {
    ENTER_FCN( "TimeStepping::TimeStepping" );

    algsys_  = algebraicsystem;
    rhsSize_ = rhsSize;
  }

  TimeStepping::~TimeStepping()
  {
    ENTER_FCN( "TimeStepping::~TimeStepping" );
  }


} // end of namespace
