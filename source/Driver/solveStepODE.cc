#include "solveStepODE.hh"

#include "PDE/bubblePDE.hh"

namespace CoupledField {

  SolveStepODE::SolveStepODE(StdPDE& apde) : StdSolveStep(apde) {
    ENTER_FCN( "SolveStepODE::SolveStepODE" );
  }
  
  SolveStepODE::~SolveStepODE() {
    ENTER_FCN( "SolveStepODE::~SolveStepODE" );
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

