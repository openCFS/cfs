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

