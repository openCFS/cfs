// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPACOUSTIC
#define FILE_SOLVESTEPACOUSTIC

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Magnetics

  class SolveStepAcoustic : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepAcoustic(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepAcoustic();


    //----------------------- TRANSIENT---------------------------------------
    //! solves for one nonlinear transient step 
    void StepTransNonLin();

    //! compute nonlinear part of RHS
    void AddNonLinRHS();

  private:


  };

} // end of namespace

#endif

