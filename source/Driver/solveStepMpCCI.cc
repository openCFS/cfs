// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "PDE/StdPDE.hh"
#include "DataInOut/ResultCache.hh"

#include "solveStepMpCCI.hh"


namespace CoupledField {

  SolveStepMpCCI::SolveStepMpCCI(StdPDE& apde) : StdSolveStep(apde)
  {

  }


  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMpCCI:: PreStepStatic()
  {
  }
  
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepMpCCI::PreStepTrans()
  {
    ResultCache::SetStepValue( actTime_ );
  }
  
  void SolveStepMpCCI::PostStepTrans()
  { 
  }




  //   Default Destructor
  // **********************
  SolveStepMpCCI::~SolveStepMpCCI() {

 
  }

} // end of namespace

