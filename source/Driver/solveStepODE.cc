// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>

#include "Driver/stdSolveStep.hh"
#include "General/exception.hh"
#include "PDE/StdPDE.hh"
#include "solveStepODE.hh"

namespace CoupledField {

  SolveStepODE::SolveStepODE(StdPDE& apde) : StdSolveStep(apde) {
  }
  
  SolveStepODE::~SolveStepODE() {
  }
  
  void SolveStepODE::SolveStepTrans() {
    EXCEPTION( "Don't know how to do a SolveStepTrans for " << PDE_.GetName());
  }


} // end of namespace

