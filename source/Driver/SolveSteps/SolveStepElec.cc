#include <fstream>
#include <iostream>
#include <string>

#include "SolveStepElec.hh"
#include "Driver/Assemble.hh"
#include "Materials/Models/Preisach.hh"
#include "PDE/StdPDE.hh"

#include "OLAS/algsys/AlgebraicSys.hh"

namespace CoupledField {

  SolveStepElec::SolveStepElec(StdPDE& apde) : StdSolveStep(apde) {
    doInit_ = true;
  }
  
  
  SolveStepElec::~SolveStepElec() {
  }
  

  // ======================================================
  // Solve Step Static SECTION  
  // ======================================================

  void SolveStepElec::PreStepStatic() {


    // Set right-hand side to zero
    //
    // Note: Though this is PreStepStatic, this is most important
    //       for transient analysis ;-)
    algsys_->InitRHS();

  }

  void SolveStepElec::PreStepTrans() {
    // Update moving ncInterfaces as needed
    ptgrid_->MoveNcInterfaces();
    
    PreStepStatic();
  }


  // time is used for a series of static calculations
  // don't get confused with REAL transient simulations!
  void SolveStepElec::SolveStepStatic() {
    if ( !isHyst_ ) 
      StepStaticLin();
  }

} // end of namespace

