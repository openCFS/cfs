// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "baseSolveStep.hh"


namespace CoupledField {

  BaseSolveStep::BaseSolveStep()
  {

    actStep_ = 0;
    actTime_ = 0.0;
    actFreq_= 0.0;
  }


  
  //   Default Destructor
  // **********************
  BaseSolveStep::~BaseSolveStep() 
  {
 
  }


} // end of namespace

