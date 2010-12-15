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
    void StepStaticNonLin(PtrParamNode analysis_id)
    {StepTransNonLin(analysis_id);};

    //----------------------- TRANSIENT---------------------------------------
    //! routine for initilizations befor execution the SolveStep-method

    //! solves for one nonlinear transient step 
    void StepTransNonLin(PtrParamNode analysis_id);

  private:
    //!< time dependency  of fluidMech PDE (important in case of analysisType=transient)
    // because the fluid can therin be determined staionary or instationary
    bool isInstationary_;
    bool isTrapezoidal_;
    Vector<Double>  oldVel_, oldSol_;

    inline Double NormL2(const Vector<Double>& vec) const
    {
      Double result(0.0);
      const UInt endLoop = vec.GetSize();

#pragma omp parallel for reduction(+:result)
      for(UInt k = 0; k < endLoop; ++k)
      {
        const Double* const dataPtr = &vec[k];
        result += (*dataPtr) * (*dataPtr);
      }

      return std::sqrt(result);
    }
  };

} // end of namespace

#endif

