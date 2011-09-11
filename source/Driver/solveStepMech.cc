#include <fstream>
#include <iostream>
#include <string>

#include "solveStepMech.hh"
#include "assemble.hh"
#include "PDE/StdPDE.hh"
#include "OLAS/algsys/algebraicSys.hh"

namespace CoupledField {

  SolveStepMech::SolveStepMech(StdPDE& apde) : StdSolveStep(apde) {
    //as this class was identical to its BaseClass all functions are removed
  }
  
  
  SolveStepMech::~SolveStepMech() {
  }

} // end of namespace

