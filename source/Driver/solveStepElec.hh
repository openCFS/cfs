// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPELEC
#define FILE_SOLVESTEPELEC

#include "stdSolveStep.hh"
#include "Utils/hysteresis.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepElec : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepElec(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepElec();


    //----------------------- STATIC---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepStatic( );

    //! base method for solving one static step 
    void SolveStepStatic(PtrParamNode analysis_id);

    //! solves for one nonlinear static step 
    void StepStaticNonLinEpsDiff(PtrParamNode analysis_id);

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans( )
    {PreStepStatic();};

    //! base method for solving one transient step 
    void SolveStepTrans(PtrParamNode analysis_id)
    {SolveStepStatic(analysis_id);};

    //! solves for one linear transient step 
    void StepTransLin(PtrParamNode analysis_id)
    {StepStaticLin(analysis_id);};

    //! routine for actions after the SolveStep-method
     void PostStepTrans()
    {PostStepStatic();};

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

  private:

    bool doInit_;

    Hysteresis * hyst_;

    Double Esat_;
    Double Psat_;
    Double Ec_;
    UInt Dir_;
    bool isVirgin_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

