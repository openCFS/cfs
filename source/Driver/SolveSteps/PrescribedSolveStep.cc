#include "PrescribedSolveStep.hh"

#include "FeBasis/FeFunctions.hh"

namespace CoupledField
{

  PrescribedSolveStep::PrescribedSolveStep(StdPDE& apde)
    : StdSolveStep(apde)
  {
  }

  PrescribedSolveStep::~PrescribedSolveStep()
  {
  }

  void PrescribedSolveStep::ApplyPrescribed()
  {
    // Write the prescribed external field(s) straight into the solution vector(s).
    // ApplyExternalData() maps the registered CoefFunction(s) onto the FeSpace and fills
    // the fe-function coefficients directly -- it touches neither the algebraic system nor
    // the matrix graph. If no external source is registered for a function, it is a no-op.
    for (auto& entry : feFunctions_)
    {
      entry.second->ApplyExternalData();
    }
  }

} // end of namespace
