#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMpCCI.hh"


namespace CoupledField {

  SolveStepMpCCI::SolveStepMpCCI(BasePDE& apde) : BaseSolveStep(apde)
  {

    ENTER_FCN( "SolveStepMpCCI::SolveStepMpCCI" );
  }


  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepMpCCI:: PreStepStatic(const Integer kstep, const Double asteptime,
				const Integer level, const Boolean reset)
  {
    ENTER_FCN( "SolveStepMpCCI::PreStepStatic" );
  }
  
  void SolveStepMpCCI::PostStepStatic(const Integer kstep, const Double asteptime,
				const Integer level)
  {
    ENTER_FCN( "SolveStepMpCCI::PostStepStatic" );
    
    if (pdeIsCoupled_)
      iterCoupledCounter_++;
  }
  

  // ======================================================
  // Solve Step Transient SECTION  
  // ======================================================

  void SolveStepMpCCI::PostStepTrans(const Integer kstep, const Double asteptime,
			       const Integer level)
  { 
    ENTER_FCN( "SolveStepMpCCI::PostStepTrans");
  }




  //   Default Destructor
  // **********************
  SolveStepMpCCI::~SolveStepMpCCI() {

    ENTER_FCN( "SolveStepMpCCI::~SolveStepMpCCI" );
 
  }

} // end of namespace

