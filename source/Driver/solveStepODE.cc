// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "solveStepODE.hh"

#include "PDE/bubblePDE.hh"

namespace CoupledField {

  SolveStepODE::SolveStepODE(StdPDE& apde) : StdSolveStep(apde) {
  }
  
  SolveStepODE::~SolveStepODE() {
  }
  
  void SolveStepODE::SolveStepTrans() {

    // Perform a cast into a BubblePDE 
    try {
      BubblePDE & bubblePDE = dynamic_cast<BubblePDE&>(PDE_);

    // Call Solve()-method
    bubblePDE.Solve();

    } catch ( ... ) {
      (*error) << "Could not cast pde of type'" 
               << PDE_.GetName()  << "' into a BubblePDE!";
      Error( __FILE__, __LINE__ );
    }



  }


} // end of namespace

