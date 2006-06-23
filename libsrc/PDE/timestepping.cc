#include <fstream>
#include <iostream>
#include <string>

#include "timestepping.hh"
#include "olas.hh"

namespace CoupledField {

  TimeStepping::TimeStepping(BaseSystem * algebraicsystem )
  {
    ENTER_FCN( "TimeStepping::TimeStepping" );

    algsys_  = algebraicsystem;
    rhsSize_ = 0 ;
  }

  TimeStepping::~TimeStepping()
  {
    ENTER_FCN( "TimeStepping::~TimeStepping" );
  }


} // end of namespace
