// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPACOUSTIC
#define FILE_SOLVESTEPACOUSTIC

#include <stddef.h>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "stdSolveStep.hh"

namespace CoupledField {
class AdjointParameters;
class StdPDE;
}  // namespace CoupledField

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcoustic : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcoustic(StdPDE& apde, bool justInterpolate = false);

    //! Destructor
    virtual ~SolveStepAcoustic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    void StepTransNonLin(PtrParamNode analysis_id);

    void StepTransLin(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);
    

    //! compute nonlinear part of RHS
    void AddNonLinRHS();

  private:

    bool justInterpolate_;
  };

} // end of namespace

#endif

