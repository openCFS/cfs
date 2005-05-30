#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMpCCI.hh"


namespace CoupledField {

  SolveStepMpCCI::SolveStepMpCCI(StdPDE& apde) : StdSolveStep(apde)
  {

    ENTER_FCN( "SolveStepMpCCI::SolveStepMpCCI" );
  }


  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMpCCI:: PreStepStatic(const UInt kstep, const Double asteptime,
                                      const Boolean reset)
  {
    ENTER_FCN( "SolveStepMpCCI::PreStepStatic" );
  }
  
  void SolveStepMpCCI::PostStepStatic(const UInt kstep, const Double asteptime)
  {
    ENTER_FCN( "SolveStepMpCCI::PostStepStatic" );
    
    if (isIterCoupled_)
      iterCoupledCounter_++;
  }
  

  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepMpCCI::PostStepTrans(const UInt kstep, const Double asteptime)
  { 
    ENTER_FCN( "SolveStepMpCCI::PostStepTrans");
  }




  //   Default Destructor
  // **********************
  SolveStepMpCCI::~SolveStepMpCCI() {

    ENTER_FCN( "SolveStepMpCCI::~SolveStepMpCCI" );
 
  }

} // end of namespace

