#include <fstream>
#include <iostream>
#include <string>

#include "BaseSolveStep.hh"
#include "Domain/Domain.hh"

namespace CoupledField {

  BaseSolveStep::BaseSolveStep()
  {
    actStep_ = 0;
    actTime_ = 0.0;
    actFreq_ = 0.0;
    couplingIter_ = 0;
  }

  
  //   Default Destructor
  // **********************
  BaseSolveStep::~BaseSolveStep() 
  {
 
  }


} // end of namespace

