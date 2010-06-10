// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPMICROPIEZO
#define FILE_SOLVESTEPMICROPIEZO

#include "stdSolveStep.hh"


namespace CoupledField
{

  //! Base class for solution of a single step: Piezotrostatic-PDE

  class SolveStepMicroPiezo : public StdSolveStep
  {

  public:

    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;

    //! Constructor
    SolveStepMicroPiezo(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMicroPiezo();

    void LineSearchPM( Vector<Double>& solIncrement, Vector<Double>& actSol,
                       Vector<Double>& solPrev, Vector<Double>& solPrevderiv2,
                       Double& etaLineSearch, Double& residualL2NormOpt );


    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans();

    //! base method for solving one transient step 
    void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL, const bool reAssembleMatrices = true);

    //! solves for one nonlinear transient step (with hysteresis) 
    void StepTransNonLin(PtrParamNode analysis_id);

    //! save the current values as previous values for the next time step
    void SetPreviousVals();

  private:

    bool doInit_;

  };

} // end of namespace

#endif

