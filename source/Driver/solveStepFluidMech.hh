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
    void StepStaticNonLin()
    {StepTransNonLin();};

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
//     void PreStepTrans( )
//     {PreStepStatic();};

    //! base method for solving one transient step 
    //! \param reset true: perfrom new assembly, etc
//     void SolveStepTrans()
//     {SolveStepStatic();};

    //! solves for one linear transient step 
//     void StepTransLin( )
//     {StepStaticLin();};

    //! routine for actions after the SolveStep-method
//      void PostStepTrans()
//     {PostStepStatic();};


    //! solves for one nonlinear transient step 
    void StepTransNonLin();
//     {StepStaticNonLin();};

    //! compute nonlinear part of RHS
//     void AddNonLinRHS(){;};

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

