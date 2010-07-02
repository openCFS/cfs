// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_SOLVESTEPPIEZO
#define FILE_SOLVESTEPPIEZO

#include "stdSolveStep.hh"


namespace CoupledField
{

  //! Base class for solution of a single step: Piezotrostatic-PDE

  class SolveStepPiezo : public StdSolveStep
  {

  public:

    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;

    //! Constructor
    SolveStepPiezo(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepPiezo();


    //----------------------- TRANSIENT---------------------------------------

    //! routine for initilizations befor execution the SolveStep-method
    void PreStepTrans();

    //! base method for solving one transient step 
    void SolveStepTrans(PtrParamNode analysis_id, AdjointParameters* adjointParams = NULL);

    //! solves for one nonlinear transient step (with hysteresis) 
    void StepTransNonLinEpsDiff(PtrParamNode analysis_id);

    //! computes differential permittivity
    void SetPreviousVals4Hyst();

    //! computes differential permittivity
    void ComputeDiffEpsilon();

  private:

    bool doInit_;

    //for differential permittivity
    Vector<Double> Eprevious_;
    Vector<Double> Dprevious_;
    Matrix<Double> epsDiff_;
  };

} // end of namespace

#endif

