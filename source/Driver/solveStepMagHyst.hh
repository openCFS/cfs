// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPMAGHYST
#define FILE_SOLVESTEPMAGHYST

#include "stdSolveStep.hh"
#include "Utils/hysteresis.hh"

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
    void StepTransNonLinHysteresis();

    //! solves for one nonlinear transient step 
    //! consideres hystreresis 
    void StepTransNonLinHysteresisDiff();

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

  private:

    Vector<Double> prevLoadRHS_;
  };

} // end of namespace

#endif

