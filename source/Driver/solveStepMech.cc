#include "Driver/stdSolveStep.hh"
#include "solveStepMech.hh"

namespace CoupledField {

class StdPDE;

  SolveStepMech::SolveStepMech(StdPDE& apde) : StdSolveStep(apde) {
    //as this class was identical to its BaseClass all functions are removed
  }
  
  
  SolveStepMech::~SolveStepMech() {
  }

} // end of namespace

