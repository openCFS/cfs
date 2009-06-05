// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPFLUIDMECH
#define FILE_SOLVESTEPFLUIDMECH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: fluid mechanics

  class SolveStepFluidMech : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepFluidMech(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepFluidMech();

    //----------------------- STATIC---------------------------------------
    //! solves for one nonlinear static step 
    void StepStaticNonLin(InfoNode* analysis_id)
    {StepTransNonLin(analysis_id);};

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method

    //! solves for one nonlinear transient step 
    void StepTransNonLin(InfoNode* analysis_id);

  private:
    //!< time dependency  of fluidMech PDE (important in case of analysisType=transient)
    // because the fluid can therin be determined staionary or instationary
    bool isInstationary_;
    bool isTrapezoidal_;
    bool isNewton_;
    Vector<Double>  oldVel_, oldSol_;

  };

} // end of namespace

#endif

