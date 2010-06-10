#ifndef FILE_SOLVESTEPSMOOTH
#define FILE_SOLVESTEPSMOOTH

#include "stdSolveStep.hh"

namespace CoupledField
{

  //! Base class for solution of a single step: Electrostatic-PDE

  // as this class was identic to its BaseClass (apart from comments), so all functions are removed
  // if one need old functionality, look into the history

  class SolveStepMech : public StdSolveStep
  {

  public:

    //! Constructor
    SolveStepMech(StdPDE& apde);

    //! Destructor
    virtual ~SolveStepMech();

  private:

  };

} // end of namespace

#endif
