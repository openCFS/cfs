// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <string>

#include "PDE/StdPDE.hh"
#include "solveStepMpCCI.hh"


namespace CoupledField {

  SolveStepMpCCI::SolveStepMpCCI(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepMpCCI::SolveStepMpCCI" );
  }


  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMpCCI:: PreStepStatic()
  {
    ENTER_FCN( "SolveStepMpCCI::PreStepStatic" );
  }
  
  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepMpCCI::PostStepTrans()
  { 
    ENTER_FCN( "SolveStepMpCCI::PostStepTrans");
  }




  //   Default Destructor
  // **********************
  SolveStepMpCCI::~SolveStepMpCCI() {

    ENTER_FCN( "SolveStepMpCCI::~SolveStepMpCCI" );
 
  }

} // end of namespace

