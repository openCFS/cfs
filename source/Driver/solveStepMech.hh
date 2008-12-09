#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  class SolveStepMech : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMech(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMech();

    //virtual void PredictorStep();

    virtual void StepTransLin();

    //virtual void StepTransNonLin();

  private:
    //Double predRelax_;
    //bool calcPredOutputCoupling_;
    bool firstEi_;
    Double Ei0_;
    Vector<Double> oldInc_;
  };

} // end of namespace

#endif
