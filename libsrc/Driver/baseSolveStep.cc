#include <fstream>
#include <iostream>
#include <string>

#include "baseSolveStep.hh"


namespace CoupledField {

  BaseSolveStep::BaseSolveStep()
  {
    ENTER_FCN( "BaseSolveStep::BaseSolveStep" );

    actStep_ = 0;
    actTime_ = 0.0;
    actFreq_= 0.0;
  }


  
  //   Default Destructor
  // **********************
  BaseSolveStep::~BaseSolveStep() 
  {
    ENTER_FCN( "BaseSolveStep::~BaseSolveStep" );
 
  }


} // end of namespace

