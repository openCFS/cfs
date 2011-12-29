// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPMAGHYST
#define FILE_SOLVESTEPMAGHYST

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "MatVec/vector.hh"
#include "stdSolveStep.hh"

namespace CoupledField {
class StdPDE;
}  // namespace CoupledField

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetic PDE with hysteresis

  class SolveStepMagHyst : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMagHyst(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMagHyst();

    //! solves for one nonlinear transient step 
    //! consideres hystreresis 
    void StepTransNonLinHysteresis(PtrParamNode analysis_id);

    //! solves for one nonlinear transient step 
    //! consideres hystreresis 
    void StepTransNonLinHysteresisDiff(PtrParamNode analysis_id);

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

  private:

    Vector<Double> prevLoadRHS_;
  };

} // end of namespace

#endif

